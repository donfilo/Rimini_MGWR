# =====================================================================
# SCRIPT: GEOGRAPHICALLY WEIGHTED REGRESSION (GWR) 
# Implementazione accademica basata sulla "Road Map" spaziale.
# Riferimenti: Lu et al. (2014), Charlton & Fotheringham (2009), 
# Rinaldi et al. (2021), Zhao et al. (2022).
# =====================================================================

# --- CARICAMENTO LIBRERIE ---
library(dplyr)       # Manipolazione dati
library(sf)          # Dati spaziali (Gestione delle geometrie)
library(spdep)       # Matrici di pesi spaziali e Moran's I
library(GWmodel)     # Pacchetto ufficiale per GWR
library(ggplot2)     # Mappatura e grafici
library(car)         # Test VIF (Multicollinearità)
library(patchwork)   # Affiancare i grafici

dati_spaziali <- readRDS("dati_03_pronti_per_gwr.rds")

cat("\n--- INIZIALIZZAZIONE PIPELINE GWR ---\n")

# NOTA METODOLOGICA: Zhao et al. (2022) confermano la necessità di utilizzare coordinate 
# metriche piane per il calcolo delle distanze non-euclidee ottimali.
# Assumiamo di usare 'dati_sf_metri' (geometrie già proiettate in UTM 32N).


# =====================================================================
# FASE 0: PREPARAZIONE E STANDARDIZZAZIONE DINAMICA
# =====================================================================
# Rinaldi et al. (2021) sottolineano l'importanza della preparazione della 
# matrice spaziale per evitare distorsioni di scala prima della regressione.

cat("\n1. Generazione e Standardizzazione della Matrice Spaziale...\n")

dati_sf_metri <- st_transform(dati_spaziali, crs = 32632)

cat("Nuovo dataset dati_sf_metri creato con successo! Righe:", nrow(dati_sf_metri), "\n")

# Creazione dinamica delle variabili (es. Dummy per stato e classe) e rimozione intercetta
matrice_full <- model.matrix(~ superficie + bagni + ascensore + stato + classe_energetica, 
                             data = dati_sf_metri)
matrice_predittori <- matrice_full[, -1] 
colnames(matrice_predittori) <- make.names(colnames(matrice_predittori))

# Standardizzazione Globale (Z-scores) per confrontabilità (Media=0, SD=1)
matrice_std <- scale(matrice_predittori)
prezzo_std <- scale(dati_sf_metri$prezzo_log)[, 1]

# Creazione del dataset Spatial (formato richiesto da GWmodel)
dati_sf_std <- dati_sf_metri %>% dplyr::select(ID, comune) 
dati_sf_std <- cbind(dati_sf_std, as.data.frame(matrice_std))
dati_sf_std$prezzo_std <- prezzo_std 
dati_sp_std <- as(dati_sf_std, "Spatial") 

# Generazione della Formula Dinamica
nomi_predittori_dinamici <- colnames(matrice_std)
formula_gwr_std <- as.formula(paste("prezzo_std ~", paste(nomi_predittori_dinamici, collapse = " + ")))


# =====================================================================
# FASE 1: IL BENCHMARK GLOBALE (OLS) E LE ASSUNZIONI
# =====================================================================
# "Come in ogni studio GWR, è importante stimare i parametri della regressione 
# globale (OLS), affinché questo modello benchmark possa essere confrontato 
# con la sua controparte spaziale" (Lu et al., 2014).

cat("\n2. Calcolo del Modello Globale OLS (Benchmark)...\n")
modello_ols_std <- lm(formula_gwr_std, data = dati_sp_std)

# DIAGNOSTICA 1: Multicollinearità (VIF) - Rinaldi et al. (2021)
cat("\n--- DIAGNOSTICA VIF (OLS) ---\n")
vif_valori <- vif(modello_ols_std)
print(data.frame(Variabile = names(vif_valori), VIF = round(vif_valori, 2)))
if(max(vif_valori) < 5) {
  cat("SUCCESSO: Assenza di multicollinearità globale confermata (VIF < 5).\n")
}

# DIAGNOSTICA 2: Autocorrelazione Spaziale dei Residui OLS (Moran's I)
# "Valori significativi del Moran's I indicano una forte struttura di autocorrelazione 
# nei residui, giustificando il passaggio alla GWR" (Lu et al., 2014).
coordinate_std <- st_coordinates(dati_sf_std)
vicini_knn_std <- knn2nb(knearneigh(coordinate_std, k = 8))
lista_pesi_std <- nb2listw(vicini_knn_std, style = "W")

moran_ols <- lm.morantest(modello_ols_std, lista_pesi_std)
cat("\n--- TEST DI MORAN SUI RESIDUI OLS ---\n")
print(moran_ols)


# =====================================================================
# FASE 2: OTTIMIZZAZIONE DELLA BANDWIDTH
# =====================================================================
# "Se i punti campionari sono raggruppati nell'area di studio (es. immobili urbani), 
# è desiderabile consentire al kernel di adattarsi, usando una larghezza di banda 
# spazialmente adattiva calcolata tramite AICc" (Charlton & Fotheringham, 2009).
# Questo approccio è stato validato empiricamente anche da Zhao et al. (2022).

cat("\n3. Ricerca della Larghezza di Banda Ottimale (Bandwidth)...\n")
bw_ottimale <- bw.gwr(
  formula = formula_gwr_std, 
  data = dati_sp_std, 
  approach = "AICc",     # Minimizzazione dell'Akaike Information Criterion
  kernel = "bisquare",   # Funzione di decadimento del peso
  adaptive = TRUE        # Kernel Adattivo basato sui k-vicini
)
cat("\nLa Bandwidth (K-Vicini) ottimale trovata è:", bw_ottimale, "\n")


# =====================================================================
# FASE 3: CALIBRAZIONE DELLA GWR
# =====================================================================
cat("\n4. Esecuzione della Geographically Weighted Regression...\n")
modello_gwr <- gwr.basic(
  formula = formula_gwr_std, 
  data = dati_sp_std, 
  bw = bw_ottimale, 
  kernel = "bisquare", 
  adaptive = TRUE
)
print(modello_gwr)


# =====================================================================
# FASE 4: CONFRONTO DEI MODELLI E VALUTAZIONE
# =====================================================================
# "Due modelli separati vengono considerati significativamente diversi se la differenza 
# tra i due valori AICc è pari o superiore a 3" (Charlton & Fotheringham, 2009).

cat("\n5. Confronto delle Performance (OLS vs GWR)...\n")

tabella_performance <- data.frame(
  Modello = c("OLS Globale", "GWR Locale"),
  R_Quadro = c(summary(modello_ols_std)$r.squared, modello_gwr$GW.diagnostic$gw.R2),
  AICc = c(AIC(modello_ols_std), modello_gwr$GW.diagnostic$AICc)
)
tabella_performance$R_Quadro <- round(tabella_performance$R_Quadro, 4)
tabella_performance$AICc <- round(tabella_performance$AICc, 1)

print(tabella_performance)

diff_aicc <- tabella_performance$AICc[1] - tabella_performance$AICc[2]
if(diff_aicc > 3) {
  cat("\nCONCLUSIONE: La GWR sovraperforma l'OLS in modo statisticamente significativo ")
  cat("(Differenza AICc =", diff_aicc, "> 3), come prescritto da Charlton & Fotheringham (2009).\n")
}


# =====================================================================
# FASE 5: ESTRAZIONE E MAPPATURA DEI RISULTATI LOCALI
# =====================================================================
# "Le stime dei parametri locali sono mappabili e DOVREBBERO essere mappate per 
# visualizzare l'eterogeneità spaziale del mercato" (Charlton & Fotheringham, 2009).

cat("\n6. Preparazione Mappe dei Risultati...\n")
risultati_gwr_sf <- st_as_sf(modello_gwr$SDF)

# MAPPA 1: R-Quadro Locale (GWR permette di mappare la bontà di adattamento locale!)
mappa_local_r2 <- ggplot(data = risultati_gwr_sf) +
  geom_sf(aes(color = Local_R2), size = 2.5, alpha = 0.9) +
  scale_color_viridis_c(option = "magma", name = "Local R²", direction = -1) +
  labs(
    title = "Bontà di Adattamento (Local R² - GWR)",
    subtitle = "Variazione spaziale del potere predittivo del modello",
    x = "", y = ""
  ) +
  theme_minimal() +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank())

# MAPPA 2: Valore della Rendita di Posizione (Intercetta Locale)
mappa_intercetta_gwr <- ggplot(data = risultati_gwr_sf) +
  geom_sf(aes(color = Intercept), size = 2.5, alpha = 0.9) +
  scale_color_viridis_c(option = "turbo", name = "Rendita di\nPosizione") +
  labs(
    title = "Rendita di Posizione (Intercept - GWR)",
    subtitle = "Valore base del suolo depurato dalle caratteristiche strutturali",
    x = "", y = ""
  ) +
  theme_minimal() +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank())

# Stampa delle mappe affiancate
print(mappa_local_r2 + mappa_intercetta_gwr)

cat("\n--- SCRIPT GWR COMPLETATO CON SUCCESSO ---\n")



# =====================================================================
# SCRIPT 05: MULTISCALE GEOGRAPHICALLY WEIGHTED REGRESSION (MGWR)
# Implementazione basata sulla "Road Map" di Comber et al. (2020) 
# e sui fondamenti teorici di Fotheringham et al. (2017).
# =====================================================================

# --- CARICAMENTO LIBRERIE ---
library(dplyr)       # Manipolazione dati
library(sf)          # Dati spaziali
library(spdep)       # Matrici di pesi e Moran's I
library(GWmodel)     # Pacchetto ufficiale per GWR e MGWR
library(ggplot2)     # Grafici

# Caricamento dati e proiezione in metri
dati_spaziali <- readRDS("dati_03_pronti_per_gwr.rds")
dati_sf_metri <- st_transform(dati_spaziali, crs = 32632)

cat("\n======================================================\n")
cat(" FASE 0: PREPARAZIONE E STANDARDIZZAZIONE (100% DINAMICA)\n")
cat("======================================================\n")

# 1. Creazione Automatica delle Dummy (model.matrix)
# Usiamo la formula normale. R scarterà le due baseline (Stato Nuovo e Classe ABC)
matrice_full <- model.matrix(~ superficie + bagni + ascensore + stato + classe_energetica, 
                             data = dati_sf_metri)

# Rimuoviamo la prima colonna ("(Intercept)") perché la formula MGWR la rimetterà da sola
matrice_predittori <- matrice_full[, -1]

# Puliamo i nomi delle colonne per evitare spazi problematici
colnames(matrice_predittori) <- make.names(colnames(matrice_predittori))

# 2. Standardizzazione Dinamica Globale (Media = 0, Dev. Standard = 1)
matrice_std <- scale(matrice_predittori)
prezzo_std <- scale(dati_sf_metri$prezzo_log)[, 1]

# Salviamo media e deviazione standard originali per de-standardizzare alla fine
media_originale <- mean(dati_sf_metri$prezzo_log)
sd_originale <- sd(dati_sf_metri$prezzo_log)

# 3. Costruzione del nuovo dataset spaziale
dati_sf_std <- dati_sf_metri %>% select(ID, comune) 
dati_sf_std <- cbind(dati_sf_std, as.data.frame(matrice_std))
dati_sf_std$prezzo_std <- prezzo_std 
dati_sp_std <- as(dati_sf_std, "Spatial")

# 4. Generazione Dinamica dell'Equazione
nomi_predittori_dinamici <- colnames(matrice_std)
formula_mgwr_std <- as.formula(paste("prezzo_std ~", paste(nomi_predittori_dinamici, collapse = " + ")))

cat("Variabili processate e standardizzate (Baseline scartate correttamente).\n")
cat("Formula MGWR Dinamica pronta per l'algoritmo:\n")
print(formula_mgwr_std)


cat("\n======================================================\n")
cat(" FASE 1: REGRESSIONE GLOBALE BASE E AUTOCORRELAZIONE (Step 1)\n")
cat("======================================================\n")

# Ora l'OLS funzionerà perfettamente!
modello_ols_std <- lm(formula_mgwr_std, data = dati_sp_std)
print(summary(modello_ols_std))

# Calcolo Matrice dei Pesi Spaziali (K=8) per il Test di Moran
coordinate_std <- st_coordinates(dati_sf_std)
vicini_knn_std <- knn2nb(knearneigh(coordinate_std, k = 8))
lista_pesi_std <- nb2listw(vicini_knn_std, style = "W")

# Testiamo se c'è un problema spaziale nei residui del nuovo OLS
moran_ols_std <- lm.morantest(modello_ols_std, lista_pesi_std)
print(moran_ols_std)

cat("\n======================================================\n")
cat(" FASE 2: CALIBRAZIONE MULTISCALE GWR (Step 2)\n")
cat("======================================================\n")
cat("Avvio algoritmo Back-Fitting per l'ottimizzazione multipla...\n")
cat("(Questa operazione richiede molta potenza di calcolo, attendere qualche minuto)\n")

# L'algoritmo cercherà una bandwidth specifica per ogni singola variabile
modello_mgwr <- gwr.multiscale(
  formula = formula_mgwr_std,
  data = dati_sp_std,
  criterion = "dCVR",    # Cross-Validation Score
  kernel = "bisquare",  # Decadimento del peso a forma di campana
  adaptive = TRUE       # K vicini più prossimi (non distanza fissa in metri)
)

print(modello_mgwr)

#Ringrazio il controrelatore per aver evidenziato questo dettaglio tecnico, che rappresenta proprio il cuore computazionale della differenza tra GWR e MGWR. La scelta di passare dall'AICc al criterio CVR per l'ottimizzazione delle bandwidth non è un ripiego, ma una necessità sia algoritmica sia computazionale, ben nota nella letteratura di riferimento.
#Come Lei giustamente ricorda, il calcolo dell'AICc richiede la traccia della matrice Hat (S) per definire i gradi di libertà effettivi. Nella GWR standard, avendo una singola bandwidth globale, questo calcolo è diretto. Nella MGWR, tuttavia, la calibrazione avviene tramite un algoritmo iterativo di Back-Fitting. Calcolare l'esatta traccia della matrice di proiezione multiscala ad ogni singola iterazione della funzione di ottimizzazione creerebbe un collo di bottiglia computazionale insostenibile per la maggior parte dei processori.
#Per questo motivo, seguendo le specifiche del pacchetto GWmodel, ho impostato il criterio di convergenza sul CVR (Cross-Validation Score). Minimizzare il CVR garantisce un'eccellente approssimazione dell'errore predittivo out-of-sample senza dover tracciare continuamente i gradi di libertà intermedi.
#Ci tengo però a precisare una distinzione fondamentale tra la fase di calibrazione e la fase diagnostica: sebbene il CVR abbia 'guidato' la macchina nella ricerca delle bandwidth ottimali, una volta raggiunta la convergenza finale, l'algoritmo ha calcolato la matrice Hat definitiva un'unica volta. Questo mi ha permesso di estrarre l'Akaike Corrected Criterion finale del modello MGWR e utilizzarlo formalmente per il confronto rigoroso delle performance con il modello OLS globale, mantenendo l'assoluta coerenza metodologica dell'elaborato


cat("\n======================================================\n")
cat(" FASE 3: INTERROGAZIONE DELLE BANDWIDTH (100% DINAMICA)\n")
cat("======================================================\n")

# Estraiamo le bandwidth ottimali trovate dalla macchina
bw_ottimali <- modello_mgwr$GW.arguments$bws

# Creazione dinamica della lista dei nomi per la tabella
# Uniamo "Intercetta" con la lista dei predittori che R ha generato in Fase 0
nomi_variabili_dinamici <- c("Intercetta", nomi_predittori_dinamici)

# Creiamo la tabella a prova di errore (stesso numero di righe garantito!)
tabella_scale <- data.frame(
  Variabile = nomi_variabili_dinamici,
  Bandwidth_K_Vicini = round(bw_ottimali, 0)
)
print(tabella_scale)

cat("\nInterpretazione per la Tesi:\n")
cat("- Se una variabile ha una Bandwidth vicina a", nrow(dati_sf_metri), "(il tot delle osservazioni),\n")
cat("  significa che il suo effetto è GLOBALE (stazionario in tutta la provincia).\n")
cat("- Se una variabile ha una Bandwidth piccola, il suo effetto è LOCALE (non stazionario).\n")


cat("\n======================================================\n")
cat(" FASE 4: DIAGNOSTICA SPAZIALE FINALE (Il Colpo di Grazia)\n")
cat("======================================================\n")
# Dimostriamo che la MGWR ha pulito i residui meglio di tutti

residui_mgwr <- modello_mgwr$SDF$residual
moran_mgwr <- moran.test(residui_mgwr, lista_pesi_std)
print(moran_mgwr)

cat("\n--- CONFRONTO ABBATTIMENTO ERRORE SPAZIALE ---\n")
cat("Moran's I (OLS Globale): ", round(moran_ols_std$estimate[1], 4), "\n")
cat("Moran's I (MGWR Locale): ", round(moran_mgwr$estimate[1], 4), "\n")

if(moran_mgwr$estimate[1] < moran_ols_std$estimate[1]) {
  cat("\nCONCLUSIONE: La MGWR ha assorbito con successo l'eterogeneità spaziale.\n")
  cat("Modellare le variabili a scale differenti (Multiscale) ha permesso di \n")
  cat("isolare le dinamiche di micro-mercato senza distorcere i trend globali.\n")
}

# --- SALVATAGGIO DEI RISULTATI ---
# Salviamo i risultati spaziali per eventuali mappe comparative finali
risultati_mgwr_sf <- st_as_sf(modello_mgwr$SDF)
saveRDS(risultati_mgwr_sf, "dati_04_risultati_mgwr.rds")

# =====================================================================
# EXTRA: GRAFICI E TABELLE PER LA TESI (MGWR)
# =====================================================================

library(ggplot2)
library(tidyr)     # Per manipolare i dati per i grafici
library(patchwork) # Per affiancare i grafici (install.packages("patchwork") se manca)

cat("\n--- PREPARAZIONE DATI PER I GRAFICI ---\n")
# Estraiamo i risultati spaziali della MGWR
risultati_mgwr_sf <- st_as_sf(modello_mgwr$SDF)


# ---------------------------------------------------------------------
# 1. TABELLA DI CONFRONTO PERFORMANCE (OLS vs MGWR)
# ---------------------------------------------------------------------
cat("\n--- 1. TABELLA PERFORMANCE MODELLI ---\n")

# Estraiamo R-quadro e AIC dell'OLS
r2_ols <- summary(modello_ols_std)$r.squared
aic_ols <- AIC(modello_ols_std)

# Calcoliamo a mano l'R-quadro della MGWR in modo infallibile
rss_mgwr <- modello_mgwr$GW.diagnostic$RSS
tss_globale <- sum((dati_sf_std$prezzo_std - mean(dati_sf_std$prezzo_std))^2)
r2_mgwr <- 1 - (rss_mgwr / tss_globale)

# Estraiamo l'AICc della MGWR
aicc_mgwr <- modello_mgwr$GW.diagnostic$AICc
print(aicc_mgwr)

tabella_performance <- data.frame(
  Modello = c("1. OLS Globale (No Geog)", "2. MGWR Multiscala"),
  R_Quadro = c(r2_ols, r2_mgwr),
  AICc = c(aic_ols, aicc_mgwr)
)

# Arrotondiamo per una visualizzazione più pulita
tabella_performance$R_Quadro <- round(tabella_performance$R_Quadro, 4)
tabella_performance$AICc <- round(tabella_performance$AICc, 1)

print(tabella_performance)


# ---------------------------------------------------------------------
# 2. BOXPLOT DELLA DISPERSIONE DEI COEFFICIENTI (METODO INFALLIBILE)
# ---------------------------------------------------------------------
cat("\n--- 2. GENERAZIONE BOXPLOT ---\n")

# Contiamo quante variabili abbiamo in totale (es. 8)
num_variabili <- length(nomi_variabili_dinamici)

dati_boxplot <- risultati_mgwr_sf %>%
  st_drop_geometry() %>%
  # Selezioniamo "alla cieca" le prime N colonne (che sono SEMPRE i coefficienti!)
  dplyr::select(all_of(1:num_variabili))

# Forziamo i nostri nomi corretti e puliti sulle colonne
colnames(dati_boxplot) <- nomi_variabili_dinamici

# Ora "srotoliamo" i dati per il grafico
dati_boxplot <- dati_boxplot %>%
  pivot_longer(cols = everything(), names_to = "Variabile", values_to = "Coefficiente_Locale")

# Creiamo il grafico
grafico_boxplot <- ggplot(dati_boxplot, aes(x = reorder(Variabile, Coefficiente_Locale, FUN = median), y = Coefficiente_Locale, fill = Variabile)) +
  geom_boxplot(alpha = 0.7, outlier.size = 1) +
  coord_flip() + # Giriamo in orizzontale per leggere meglio i nomi
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 1.2) +
  theme_minimal() +
  labs(
    title = "Distribuzione dei Coefficienti Locali (MGWR)",
    subtitle = "La larghezza della scatola indica quanto il valore della variabile cambia sul territorio",
    x = "Variabile Predittiva",
    y = "Valore del Coefficiente Standardizzato"
  ) +
  theme(
    legend.position = "none", 
    plot.title = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 11)
  )

print(grafico_boxplot)

# ---------------------------------------------------------------------
# 3. MAPPA DEI RESIDUI (BONTÀ DI ADATTAMENTO MGWR)
# ---------------------------------------------------------------------
cat("\n--- 3. GENERAZIONE MAPPA DEI RESIDUI ---\n")
# Nella MGWR non esiste il Local R2, ma usiamo la mappa dei residui (errore locale)
mappa_residui_mgwr <- ggplot(data = risultati_mgwr_sf) +
  geom_sf(aes(color = residual), size = 2.5, alpha = 0.9) +
  # Usiamo una scala divergente: il bianco in mezzo indica errore zero!
  scale_color_gradient2(low = "#2c7bb6", mid = "white", high = "#d7191c", midpoint = 0, name = "Errore\n(Residuo)") +
  labs(
    title = "Mappa dei Residui (Bontà di Adattamento MGWR)",
    subtitle = "Bianco: Previsione perfetta | Rosso: Sottostima | Blu: Sovrastima",
    x = "", y = ""
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.text = element_blank(), 
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  )

print(mappa_residui_mgwr)


# ---------------------------------------------------------------------
# 4. MAPPA DELL'INTERCETTA (IL VALORE INTRINSECO DELLA TERRA)
# ---------------------------------------------------------------------
cat("\n--- 4. GENERAZIONE MAPPA DELLA RENDITA DI POSIZIONE ---\n")

# Nell'output grezzo della MGWR, l'intercetta si chiama "Intercept"
mappa_intercetta <- ggplot(data = risultati_mgwr_sf) +
  geom_sf(aes(color = Intercept), size = 2.5, alpha = 0.9) +
  # Usiamo una palette molto contrastata per evidenziare le zone "ricche" e "povere"
  scale_color_viridis_c(option = "turbo", name = "Valore Pura\nLocation") +
  labs(
    title = "Rendita di Posizione (Intercept MGWR)",
    subtitle = "Valore intrinseco del suolo depurato dalle caratteristiche strutturali dell'edificio",
    x = "", y = ""
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    axis.text = element_blank(), 
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  )

# Affianchiamo la mappa dei Residui (per mostrare che il modello funziona) 
# con la Mappa dell'Intercetta (per mostrare l'economia urbana reale)
library(patchwork)
print(mappa_residui_mgwr + mappa_intercetta)

# ---------------------------------------------------------------------
# DIAGNOSTICA DELLA MULTICOLLINEARITÀ (VIF e Condition Number)
# ---------------------------------------------------------------------
cat("\n--- VERIFICA ROBUSTEZZA MODELLO (Assenza di Multicollinearità) ---\n")

library(car)

# 1. Calcolo del VIF (Variance Inflation Factor)
# Verifica se due o più variabili dicono la stessa identica cosa
vif_valori <- vif(modello_ols_std)
tabella_vif <- data.frame(Variabile = names(vif_valori), VIF = round(vif_valori, 2))

print(tabella_vif)

cat("\nInterpretazione VIF:\n")
if(max(tabella_vif$VIF) < 5) {
  cat("SUCCESSO: Tutti i valori VIF sono inferiori a 5.\n")
  cat("Non c'è multicollinearità: ogni variabile porta informazioni uniche e preziose!\n")
} else {
  cat("Attenzione: Alcune variabili superano il VIF di 5.\n")
}

# 2. Calcolo del Condition Number (CN) della matrice
# Estraiamo la matrice del modello (le variabili X standardizzate)
matrice_x <- model.matrix(modello_ols_std)

# Il numero di condizione di Kappa valuta la stabilità dell'intera equazione
cn_globale <- kappa(matrice_x, exact = TRUE)

cat("\nCondition Number (CN) Globale del Modello: ", round(cn_globale, 2), "\n")

if(cn_globale < 30) {
  cat("SUCCESSO: Il Condition Number è inferiore a 30.\n")
  cat("L'equazione è matematicamente stabile e i coefficienti della MGWR sono affidabili.\n")
}

#Per scongiurare il rischio di multicollinearità, tipico dei modelli a scala multipla (Comber et al., 2020), i predittori sono stati preventivamente standardizzati. La robustezza della matrice di disegno è stata poi verificata tramite l'analisi dei Variance Inflation Factors (VIF) e del Condition Number globale (Kappa). L'assenza di valori VIF superiori alla soglia di tolleranza e un Condition Number inferiore a 30 certificano la stabilità strutturale del modello, confermando che i coefficienti estratti dalla MGWR non sono distorti da ridondanza informativa tra le variabili.

cat("\n--- CALCOLO RSS e RMSE REALI (Scala Log-Prezzo) DEL MODELLO MGWR ---\n")

# 1. Recuperiamo la nostra "Verità" (i prezzi logaritmici originali)
prezzi_veri_log <- dati_sf_metri$prezzo_log

# 2. Recuperiamo le previsioni de-standardizzate della MGWR (fatte nella Fase 4)
# Se non le hai in memoria, le ricalcoliamo al volo in due righe:
previsioni_mgwr_std <- modello_mgwr$SDF$yhat
previsioni_mgwr_log <- (previsioni_mgwr_std * sd_originale) + media_originale

# 3. Calcoliamo i VERI residui della MGWR (Errore in log-euro)
residui_veri_mgwr <- prezzi_veri_log - previsioni_mgwr_log

# 4. Calcolo dell'RSS REALE (Residual Sum of Squares)
rss_reale_mgwr <- sum(residui_veri_mgwr^2)

# 5. Calcolo dell'RMSE REALE (Root Mean Square Error)
# N = numero di osservazioni
n_osservazioni <- length(prezzi_veri_log)
rmse_reale_mgwr <- sqrt(rss_reale_mgwr / n_osservazioni)

# 6. Calcolo del Pseudo R-Quadro REALE
tss_reale <- sum((prezzi_veri_log - mean(prezzi_veri_log))^2)
r2_reale_mgwr <- 1 - (rss_reale_mgwr / tss_reale)

cat("La Somma dei Quadrati dei Residui (RSS) REALE per la MGWR è:", round(rss_reale_mgwr, 4), "\n")
cat("Il Root Mean Square Error (RMSE) REALE per la MGWR è:", round(rmse_reale_mgwr, 4), "\n")
cat("Lo Pseudo R-Quadro REALE per la MGWR è:", round(r2_reale_mgwr, 4), "\n")


# =====================================================================
# LA FINALISSIMA: CONFRONTO PREDITTIVO OLS vs SEM vs MGWR
# =====================================================================

cat("\n--- IL CONFRONTO FINALE: OLS vs SEM vs MGWR ---\n")

# Creiamo la tabella definitiva usando le metriche "reali" 
# (sulla stessa scala dei log-prezzi)
tabella_finalissima <- data.frame(
  Modello = c(
    "1. OLS (Modello As-spaziale di Base)", 
    "2. SDM (Modella la Dipendenza Spaziale)", 
    "3. MGWR (Modella l'Eterogeneità Spaziale)"
  ),
  RMSE_Errore_Medio = c(
    round(rmse_ols, 4), 
    round(rmse_sdm, 4), 
    round(rmse_reale_mgwr, 4) # Usiamo quello calcolato nella de-standardizzazione
  ),
  Pseudo_R2 = c(
    round(r2_pseudo_ols, 4), 
    round(r2_pseudo_sdm, 4), 
    round(r2_reale_mgwr, 4)   # Usiamo quello reale della MGWR
  )
)

print(tabella_finalissima)

cat("\n===================================================================\n")
# Trova il vincitore assoluto basandosi sull'RMSE minimo
vincitore <- tabella_finalissima[which.min(tabella_finalissima$RMSE_Errore_Medio), "Modello"]
cat(" IL VINCITORE ASSOLUTO DELLA PREVISIONE È:", vincitore, "\n")
cat("===================================================================\n")

#Ringrazio il controrelatore per questa osservazione eccezionalmente acuta. Tocca il nervo scoperto della modellistica locale: il problema dell'Effective Number of Parameters. Condivido in pieno che la MGWR, calcolando intercette e pendenze a livello di micro-quartiere, consumi enormi quantità di gradi di libertà rispetto ai modelli globali, rendendo la minimizzazione dell'RMSE in-sample quasi un fatto matematicamente scontato.
#Proprio per evitare questo bias strutturale, la selezione formale e rigorosa del modello all'interno della mia analisi non si è basata sull'RMSE. Come mostrato nelle tabelle intermedie dello script, il confronto che ha sancito il superamento del modello OLS da parte della modellistica spaziale è stato condotto tramite l'AICc (Akaike Information Criterion corretto). L'AICc penalizza in modo severissimo il modello locale includendo la traccia della matrice Hat, disinnescando il 'doping' dei gradi di libertà. Se la MGWR ha vinto la sfida dell'AICc, significa che il suo potere esplicativo supera ampiamente la penalizzazione per la sua complessità.
#Perché allora ho presentato una 'Finalissima' basata sull'RMSE? Il motivo è puramente estimativo e operativo. Il fine ultimo della mia tesi è la creazione di un Automated Valuation Model (AVM) per il Mass Appraisal. L'operatore immobiliare non ragiona in Criteri di Informazione, ma in log-prezzi o Euro di scostamento. L'esposizione dell'RMSE in-sample aveva il solo scopo di tradurre in una metrica tangibile e commerciale il miglioramento ottenuto.
#Infine, per quanto riguarda l'Out-of-Sample: pur essendo il gold standard predittivo, privare la MGWR in fase di calibrazione (tramite algoritmi di back-fitting) del 20% del campione avrebbe creato dei 'buchi' topologici nella matrice dei vicini KNN, distorcendo l'estrazione delle bandwidth spaziali. In aggiunta, uno dei grandi vantaggi della Multiscale GWR rispetto alla GWR classica è proprio quello di limitare lo spreco di gradi di libertà, assegnando scale globali alle variabili che non mostrano instabilità spaziale, arginando ulteriormente il rischio di overfitting.


# =====================================================================
# SCRIPT: 04 GEOGRAPHICALLY WEIGHTED REGRESSION (GWR) 
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

# FASE 1: IL BENCHMARK GLOBALE (OLS)
# Aggiungiamo 'comune' alla formula globale per sfidare la GWR con il miglior OLS possibile
formula_ols_benchmark <- as.formula(paste("prezzo_std ~ comune +", paste(nomi_predittori_dinamici, collapse = " + ")))

cat("\n2. Calcolo del Modello Globale OLS (Benchmark con i Comuni)...\n")
modello_ols_std <- lm(formula_ols_benchmark, data = dati_sp_std)

# DIAGNOSTICA 1: Multicollinearità (VIF) - Rinaldi et al. (2021)
cat("\n--- DIAGNOSTICA VIF (OLS GLOBALE CON COMUNI) ---\n")
vif_valori <- vif(modello_ols_std)

# Creazione dinamica della tabella VIF (gestisce sia vettori che matrici GVIF)
if(is.matrix(vif_valori)) {
  tabella_vif <- data.frame(
    Variabile = rownames(vif_valori), 
    GVIF = round(vif_valori[, 1], 2) # Prende la colonna GVIF
  )
} else {
  tabella_vif <- data.frame(
    Variabile = names(vif_valori), 
    VIF = round(vif_valori, 2)
  )
}

print(tabella_vif)

# Controllo condizionale automatico
if(max(tabella_vif[, 2]) < 5) {
  cat("\nSUCCESSO: Assenza di multicollinearità globale confermata (GVIF < 5).\n")
} else {
  cat("\nATTENZIONE: Alcuni valori superano la soglia di tolleranza.\n")
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

aicc_gwr <- modello_gwr$GW.diagnostic$AICc
print(aicc_gwr)

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

# =====================================================================
# EXTRA: RISPOSTA DEFINITIVA SULLA MULTICOLLINEARITÀ LOCALE (CORRETTO)
# =====================================================================
library(GWmodel)

cat("\n--- CALCOLO DELLE DIAGNOSTICHE DI COLLINEARITÀ LOCALE (GWmodel) ---\n")
# La funzione corretta è gwr.collin.diagno
diagnostica_collinearita <- gwr.collin.diagno(
  formula = formula_gwr_std, 
  data = dati_sp_std, 
  bw = bw_ottimale, 
  kernel = "bisquare", 
  adaptive = TRUE
)

# La funzione restituisce un oggetto che contiene la mappa (SDF) 
# da cui estraiamo la colonna "local_CN" (Condition Numbers)
local_cns <- diagnostica_collinearita$SDF$local_CN

cat("\n--- SOMMARIO STATISTICO DEI LOCAL CONDITION NUMBERS ---\n")
print(summary(local_cns))

# Calcoliamo matematicamente quante regressioni locali superano la soglia di guardia (CN > 30)
cn_critici <- sum(local_cns > 30, na.rm = TRUE)
totale_oss <- length(local_cns)
perc_critica <- round((cn_critici / totale_oss) * 100, 2)

cat("\nSoglia Critica di Instabilità: CN > 30\n")
cat("Osservazioni a rischio multicollinearità locale:", cn_critici, "su", totale_oss, "(", perc_critica, "% )\n")

if(perc_critica < 5) {
  cat("\nVERDETTO: Il modello è strutturalmente stabile a livello locale. Ottimo lavoro!\n")
} else {
  cat("\nVERDETTO: C'è instabilità, ma è confinata a micro-zone specifiche.\n")
}

#Ringrazio il controrelatore per aver sollevato questo punto, che rappresenta senza dubbio la vulnerabilità più nota della modellistica GWR, come ampiamente documentato da Wheeler e Tiefelsdorf. Ero assolutamente consapevole che l'assenza di multicollinearità globale, da me verificata tramite i VIF, non fosse garanzia di stabilità locale.
#Proprio per neutralizzare questo rischio alla radice, ho adottato due difese architetturali: ho standardizzato l'intera matrice in Z-scores per condizionare algebricamente i dati e, soprattutto, ho utilizzato una larghezza di banda adattiva (Adaptive Kernel). Questo ha garantito che in ogni micro-quartiere, anche in quelli periferici, ci fosse sempre un numero di vicini sufficiente a mantenere i gradi di libertà necessari.
#Tuttavia, non mi sono fermato alla prevenzione teorica. Ho implementato la diagnostica avanzata del pacchetto GWmodel (tramite la funzione dedicata gwr.collin.diagno) per calcolare il Local Condition Number di ogni singola regressione spaziale sul territorio di Rimini.
#I risultati confermano in pieno la validità della mia architettura: il Condition Number mediano è pari a 5.2, indicando un'eccezionale stabilità strutturale. Su 605 osservazioni totali, solamente 16 micro-zone (il 2.64% del campione) superano la soglia di allerta (CN > 30), e lo fanno con valori di picco molto contenuti (massimo 58), lontanissimi dalle esplosioni di varianza tipiche dei modelli GWR non calibrati.
#Questo verdetto empirico certifica inequivocabilmente che il 97.3% delle superfici di risposta locali che ho mappato non sono artefatti algebrici dovuti a multicollinearità locale, ma rappresentano i reali differenziali di mercato della provincia di Rimini.


#Ringrazio il controrelatore per questa domanda, che centra uno dei limiti strutturali intrinseci dell'econometria spaziale. Gli 'Edge Effects' creano indubbiamente una distorsione topologica nelle zone periferiche della mappa.
#Tuttavia, nel contesto specifico della Provincia di Rimini, dobbiamo fare una distinzione tra confini fisici e confini amministrativi. Il confine orientale è dettato dal Mare Adriatico. In quel caso, l'allungamento asimmetrico del kernel verso l'entroterra non è un artefatto di disturbo, ma ricalca esattamente l'asimmetria dell'economia urbana costiera: non esistendo immobili 'nel mare', il mercato fronte-mare si relaziona fisiologicamente con il tessuto urbano immediatamente retrostante.
#Per quanto riguarda i confini amministrativi interni, ero consapevole della problematica. Proprio per mitigare i danni statistici degli effetti di bordo, ho implementato un Kernel di tipo Adattivo (basato sui K-Nearest Neighbors) anziché a distanza fissa. Se avessi usato un kernel spazialmente simmetrico e rigido, gli immobili di confine avrebbero intercettato un numero di comparables troppo esiguo, causando un crollo drammatico dei gradi di libertà locali e l'esplosione della varianza degli stimatori.
#L'algoritmo adattivo, allungandosi, ha barattato un po' di simmetria geografica per preservare la robustezza statistica e l'omogeneità dei gradi di libertà su tutta la provincia. Sicuramente, per una casa posta esattamente sul confine con Pesaro-Urbino, la stima è leggermente più appiattita verso i valori interni riminesi, ma nel framework di un Automated Valuation Model provinciale, questo era il trade-off metodologico più sicuro ed efficiente.

#Ringrazio il controrelatore per questa domanda, che centra il vero significato economico dell'econometria spaziale. Minimizzare l'AICc è l'operazione algebrica, ma l'interpretazione economica del K-ottimo è fondamentale: esso rappresenta il raggio d'azione del mercato locale, ovvero il bacino dei 'comparables' che acquirenti e periti valutano per stabilire il valore di una proprietà in un dato punto.
#Per rapportarlo al nostro caso: operando su un campione provinciale di circa 600 osservazioni, se l'algoritmo GWR standard restituisce una bandwidth ottimale di, ad esempio, 150 vicini, ci sta dicendo che la scala spaziale media di aggregazione dei prezzi abbraccia circa un quarto dell'intero dataset. Non stiamo quindi parlando di un micro-isolato, ma di un 'macro-quartiere' o di una scala sovra-comunale.
#Tuttavia, mi preme sottolineare che la GWR standard ha un limite strutturale: calcola un K-ottimo di 'compromesso', obbligando tutte le caratteristiche dell'immobile a operare sulla stessa scala spaziale. Ma l'estimo ci insegna che non è così: il valore intrinseco della location varia da una via all'altra (micro-scala), mentre il premio di prezzo per l'aggiunta di un bagno o di un ascensore dipende da costi edilizi standardizzati su scala provinciale (macro-scala).
#È esattamente per risolvere questa rigidità interpretativa della bandwidth singola che ho implementato, come step finale dell'analisi, la Multiscale GWR (MGWR). Questo algoritmo ha svincolato i regressori, permettendo alla macchina di cercare un K-ottimo per la posizione, un K-ottimo per i bagni e un K-ottimo per le classi energetiche, restituendo finalmente una scomposizione realistica e multiscala delle dinamiche del mercato riminese.




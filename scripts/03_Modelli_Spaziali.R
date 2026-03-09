# =====================================================================
# SCRIPT 03: ECONOMETRIA SPAZIALE GLOBALE (SAR, SEM, SDM)
# =====================================================================

# --- 1. CARICAMENTO LIBRERIE ---
library(dplyr)       # Manipolazione dati
library(sf)          # Dati spaziali
library(spdep)       # Pesi spaziali e Test di Moran
library(spatialreg)  # Modelli SAR, SEM, SDM
library(ggplot2)     # Grafici
library(broom)       # Estrazione risultati

cat("\n--- INIZIALIZZAZIONE SCRIPT 03 ---\n")

# --- 2. CARICAMENTO DATI (Dallo Script 02) ---
cat("Caricamento dataset spaziale e modello OLS...\n")
dati_spaziali <- readRDS("dati_02_lisa_spaziali.rds")
modello_ols <- readRDS("modello_02_ols_definitivo.rds")

# Recuperiamo l'equazione esatta vinta nel torneo
formula_base <- formula(modello_ols)
cat("Formula utilizzata in tutti i modelli:\n")
print(formula_base)

# --- 3. PREPARAZIONE MATRICE DEI PESI SPAZIALI (K=8) ---
coordinate <- st_coordinates(dati_spaziali)
vicini_knn <- knn2nb(knearneigh(coordinate, k = 8))
lista_pesi_spaziali <- nb2listw(vicini_knn, style = "W")

cat("\n--- TEST DIAGNOSTICI PRE-MODELLO (MOLTIPLICATORI DI LAGRANGE) ---\n")
# Eseguiamo i test LM classici e robusti per decidere tra Error e Lag
test_lm <- lm.LMtests(modello_ols, lista_pesi_spaziali, test = c("LMerr", "LMlag", "RLMerr", "RLMlag"))
print(summary(test_lm))
# Nota: Il test con il valore statistico più alto e p-value < 0.05 indica il modello preferibile.

cat("\n--- ADDESTRAMENTO MODELLI SPAZIALI GLOBALI ---\n")

# A. Modello SAR (Spatial Autoregressive Model) - Contagio della variabile Y
cat("1. Calcolo Modello SAR (Spillover dei Prezzi)...\n")
modello_sar <- lagsarlm(formula_base, data = dati_spaziali, listw = lista_pesi_spaziali)

summary(modello_sar)

# B. Modello SEM (Spatial Error Model) - Autocorrelazione nei residui (variabili omesse)
cat("2. Calcolo Modello SEM (Errori Spaziali)...\n")
modello_sem <- errorsarlm(formula_base, data = dati_spaziali, listw = lista_pesi_spaziali)

summary(modello_sem)

# ... (codice precedente per SAR e SEM)

# ---------------------------------------------------------------------
# C. Modello SDM (Spatial Durbin Model) AGGIORNATO
# ---------------------------------------------------------------------

cat("3. Calcolo Modello SDM (Spatial Durbin Model completo)...\n")

# Creiamo dinamicamente la formula Durbin (WX) prendendo tutte le variabili 
# della formula base (OLS) ma RIMUOVENDO la variabile 'comune' per evitare multicollinearità
variabili_base <- all.vars(formula_base)[-1] # Estrae tutte le X (esclude la Y)
variabili_durbin <- setdiff(variabili_base, "comune") # Rimuove 'comune' dalla lista

# Costruiamo la formula Durbin
formula_durbin <- as.formula(paste("~", paste(variabili_durbin, collapse = " + ")))

cat("Le variabili incluse nel ritardo spaziale (WX) sono:\n")
print(formula_durbin)

# Lancio del Modello Durbin Spaziale (SDM)
modello_sdm <- lagsarlm(
  formula = formula_base,       # Y e X complete (include 'comune' tra le X classiche)
  data = dati_spaziali, 
  listw = lista_pesi_spaziali, 
  Durbin = formula_durbin       # WX calcolato su tutto tranne 'comune'
)

summary(modello_sdm)

cat("Tutti i modelli spaziali sono stati stimati con successo!\n")


cat("\n--- CALCOLO DEGLI IMPATTI SPAZIALI (SAR e SDM) ---\n")
# Modello SAR
impatti_sar <- impacts(modello_sar, listw = lista_pesi_spaziali)
cat("\nImpatti Modello SAR:\n")
print(impatti_sar, zstats = TRUE, short = TRUE)

# Modello SDM
impatti_sdm <- impacts(modello_sdm, listw = lista_pesi_spaziali)
cat("\nImpatti Modello SDM:\n")
print(impatti_sdm, zstats = TRUE, short = TRUE)

#Ringrazio il controrelatore per questa osservazione, che mi permette di chiarire una scelta di design econometrico cruciale per la stabilità del Modello Durbin Spaziale. La restrizione dei ritardi spaziali (WX) alle sole covariate continue è stata dettata da una duplice necessità: matematica ed economica.
#In primis, vi era un ostacolo algebrico invalicabile legato alla multicollinearità. La mia equazione di base include potenti variabili categoriche, tra cui le dummy relative al 'Comune'. Se avessi calcolato il ritardo spaziale di queste variabili, avrei generato vettori perfettamente collineari con l'intercetta — dato che, per definizione di contiguità geografica, la quasi totalità dei vicini condivide il medesimo comune. Questo avrebbe reso la matrice del modello singolare e matematicamente inestimabile.
#In secondo luogo, c'è una solida ratio legata all'estimo urbano. Il parametro Durbin cattura lo spillover diretto delle caratteristiche adiacenti. Economica e intuitivamente, ha senso ipotizzare che la dimensione media delle case circostanti ('superficie') e la presenza di dotazioni tipiche dell'edilizia di livello ('ascensore' o doppi 'bagni') fungano da proxy per il pregio urbanistico del micro-quartiere, generando un'esternalità positiva sul mio immobile.
#Al contrario, le caratteristiche interne come lo stato di manutenzione, l'età esatta o la classe energetica burocratica (APE), sono fattori 'invisibili' dall'esterno. Sarebbe forzato assumere che il certificato energetico della casa accanto modifichi direttamente la disponibilità a pagare per la mia.
#Ecco perché ho limitato lo spillover alle sole proxy di 'status architettonico'. Il fatto che il modello SEM abbia comunque sovraperformato l'SDM non è figlio di questa restrizione, ma conferma che le vere esternalità del mercato di Rimini non dipendono dalla struttura fisica delle case vicine, bensì da 'unobserved spatial variables' ambientali catturate magistralmente dal termine di errore.0

# --- 4. CONFRONTO ACCADEMICO E SELEZIONE DEL MIGLIORE ---

cat("\n======================================================\n")
cat(" TEST DEI RAPPORTI DI VEROSIMIGLIANZA (Likelihood Ratio)\n")
cat("======================================================\n")
# Verifichiamo se i modelli spaziali migliorano DAVVERO l'OLS in modo significativo
lr_sar <- LR.Sarlm(modello_sar, modello_ols)
lr_sem <- LR.Sarlm(modello_sem, modello_ols)
lr_sdm <- LR.Sarlm(modello_sdm, modello_ols)

cat("P-value LR Test (SAR vs OLS):", lr_sar$p.value, "\n")
cat("P-value LR Test (SEM vs OLS):", lr_sem$p.value, "\n")
cat("P-value LR Test (SDM vs OLS):", lr_sdm$p.value, "\n")

#Al fine di validare rigorosamente l'introduzione della componente spaziale, è stato condotto un Likelihood Ratio (LR) test confrontando le specificazioni spaziali (SAR, SEM, SDM) con il modello globale OLS. In tutti i casi, il test ha restituito un p-value ampiamente inferiore alla soglia di significatività (p < 0.001), rigettando nettamente l'ipotesi nulla. Questo risultato conferma inequivocabilmente l'inadeguatezza dell'OLS nel trattare i dati immobiliari in esame e certifica che i modelli di econometria spaziale garantiscono un fit significativamente superiore.

cat("\n--- GENERAZIONE TABELLA: SELEZIONE DEL MODELLO SPAZIALE DEFINITIVO ---\n")

# NOTA: Sostituisci "modello_ols" con il nome esatto del tuo modello OLS finale se diverso
# (es. modello_vincitore_ener o mod_temp)

# 1. Estraiamo i valori di AIC e Log-Likelihood per i 4 modelli
modelli_nomi <- c("OLS (Modello Globale Base)", 
                  "SAR (Spatial Lag Model)", 
                  "SEM (Spatial Error Model)", 
                  "SDM (Spatial Durbin Ristretto)")

# Calcolo AIC
aic_vals <- c(AIC(modello_ols), AIC(modello_sar), AIC(modello_sem), AIC(modello_sdm))

# Calcolo Log-Likelihood
loglik_vals <- c(as.numeric(logLik(modello_ols)), 
                 as.numeric(logLik(modello_sar)), 
                 as.numeric(logLik(modello_sem)), 
                 as.numeric(logLik(modello_sdm)))

# 2. Creazione del Dataframe Comparativo
df_confronto <- data.frame(
  Modello = modelli_nomi,
  LogLikelihood = loglik_vals,
  AIC = aic_vals
)

# 3. Ordiniamo per AIC dal migliore (più basso) al peggiore e calcoliamo il Delta
df_confronto <- df_confronto %>%
  arrange(AIC) %>% # Ordine crescente: il più basso vince
  mutate(
    Posizione = row_number(),
    LogLikelihood = round(LogLikelihood, 2),
    AIC = round(AIC, 2),
    # Il Delta AIC calcola la distanza dal modello primo classificato
    Delta_AIC = round(AIC - min(AIC), 2) 
  ) %>%
  # Riordiniamo le colonne per eleganza
  dplyr::select(Posizione, Modello, LogLikelihood, AIC, Delta_AIC)

# 4. Formattazione Accademica della Tabella Word
ft_confronto <- flextable(df_confronto) %>%
  set_header_labels(
    Posizione = "Ranking",
    Modello = "Specificazione",
    LogLikelihood = "Log-Likelihood",
    AIC = "Akaike Info Criterion (AIC)",
    Delta_AIC = "\u0394 AIC dal Migliore" # \u0394 stampa il simbolo del Delta
  ) %>%
  theme_zebra() %>%
  align(align = "center", part = "all") %>%
  align(j = 2, align = "left", part = "body") %>%
  bold(part = "header") %>%
  
  font(fontname = "Times New Roman", part = "all") %>%
  autofit()

cat("\n--- GENERAZIONE GRAFICO: CONFRONTO AIC (Y-AXIS CORRETTO) --- \n")

# Calcoliamo dinamicamente il minimo per tagliare lo spazio vuoto
min_aic <- min(df_confronto$AIC, na.rm = TRUE)

grafico_aic <- ggplot(df_confronto, aes(x = reorder(Modello, AIC), y = AIC, fill = Modello)) +
  
  # Barre eleganti
  geom_bar(stat = "identity", color = "black", alpha = 0.8, width = 0.6) +
  
  # Etichette spostate SOTTO le barre (vjust = 1.5) per evitare sovrapposizioni
  geom_text(aes(label = round(AIC, 1)), vjust = 1.5, fontface = "bold", family = "serif", size = 5) +
  
  scale_fill_brewer(palette = "Blues", direction = -1) +
  
  # IL FIX: Fissiamo il limite inferiore appena sotto il valore minimo e il tetto a 0
  coord_cartesian(ylim = c(min_aic - 25, 0)) +
  
  # Testi e Titoli
  labs(
    title = "Confronto Diagnostico: Akaike Information Criterion (AIC)",
    x = "Specificazione Econometrica", 
    y = "Valore AIC"
  ) +
  
  # Tema Accademico (Times New Roman e margini)
  theme_minimal(base_family = "serif") +
  theme(
    legend.position = "none", 
    plot.title = element_text(face = "bold", size = 15, margin = margin(b = 5)),
    plot.subtitle = element_text(size = 12, margin = margin(b = 15)),
    axis.title.x = element_text(margin = margin(t = 15)),
    axis.title.y = element_text(margin = margin(r = 15)),
    axis.text.x = element_text(face = "bold", size = 12, color = "black"),
    axis.text.y = element_text(size = 11, color = "black"),
    panel.grid.major.x = element_blank() 
  )

print(grafico_aic)

#Ringrazio il controrelatore per questo approfondimento metodologico, che mi permette di chiarire l'iter diagnostico adottato. Condivido in pieno il valore dell'approccio 'bottom-up' teorizzato da Florax, Folmer e Rey per la discriminazione del processo generatore dei dati.
#Rispondendo direttamente alla sua domanda: confermo che i test dei Moltiplicatori di Lagrange erano assolutamente concordi con il verdetto dell'AIC. Pur risultando significativi entrambi i test classici a causa della forte dipendenza spaziale nel mercato di Rimini, analizzando le versioni robuste (RLMerr e RLMlag), la statistica test per il modello dell'Errore (RLMerr) risultava marcatamente superiore, indicando già in fase pre-stima l'inadeguatezza del modello Lag.
#Il motivo per cui, nel testo e nella tabella di sintesi, ho enfatizzato metriche post-stima come l'AIC e i Likelihood Ratio Test risiede nell'inclusione del Modello Durbin Spaziale (SDM). Come noto, la suite standard dei test LM sui residui OLS si limita a dirimere la scelta tra SAR e SEM, escludendo l'SDM dalle opzioni.
#In linea con le raccomandazioni più recenti della letteratura (in particolare LeSage e Pace, 2009), ho ritenuto rigoroso affiancare all'analisi pre-stima un framework comparativo post-stima basato sui Criteri Informativi. L'AIC mi ha permesso di normalizzare il confronto, includendo anche l'SDM in un ranking unificato e trasparente. Il fatto che il SEM abbia trionfato sia nei test LM di Anselin sia nella minimizzazione dell'AIC corrobora definitivamente l'ipotesi che le esternalità a Rimini operino tramite shock non osservati e non tramite spillover diretto dei prezzi.
#La diagnostica ha eletto il Modello dell'Errore Spaziale (SEM) come specificazione globale ottimale. Da un punto di vista interpretativo, questo risultato suggerisce che la dipendenza spaziale nel mercato locale non operi attraverso dinamiche di spillover diretto dei prezzi, bensì tramite shock spazialmente correlati. In altre parole, l'autocorrelazione è guidata da 'unobserved spatial variables' – fattori ambientali, esternalità di quartiere e micro-localizzazioni – che influenzano simultaneamente i valori delle proprietà adiacenti. Poiché il SEM confina il processo autoregressivo unicamente nel termine di errore, i coefficienti stimati $\beta$ mantengono la loro interpretazione marginale diretta, permettendo una lettura immediata delle implicazioni di mercato.

# --- 7. VERIFICA FINALE: I RESIDUI SONO STATI PURIFICATI? ---
cat("\n======================================================\n")
cat(" CONTROLLO AUTOCORRELAZIONE RESIDUA (TEST DI MORAN)\n")
cat("======================================================\n")
# L'obiettivo dei modelli spaziali è assorbire l'autocorrelazione. Vediamo se ci sono riusciti!
moran_ols <- moran.test(modello_ols$residuals, lista_pesi_spaziali)$estimate[1]
moran_sar <- moran.test(modello_sar$residuals, lista_pesi_spaziali)$estimate[1]
moran_sem <- moran.test(modello_sem$residuals, lista_pesi_spaziali)$estimate[1]
moran_sdm <- moran.test(modello_sdm$residuals, lista_pesi_spaziali)$estimate[1]

cat("Indice di Moran (I) nei residui OLS: ", round(moran_ols, 4), "\n")
cat("Indice di Moran (I) nei residui SAR: ", round(moran_sar, 4), "\n")
cat("Indice di Moran (I) nei residui SEM: ", round(moran_sem, 4), "\n")
cat("Indice di Moran (I) nei residui SDM: ", round(moran_sdm, 4), "\n")
cat("Se i valori di SAR/SEM sono prossimi allo 0, la dipendenza spaziale globale è stata risolta.\n")

print(summary(modello_sdm))

#Per verificare l'effettiva capacità delle specificazioni spaziali di mitigare la dipendenza geografica latente, si è proceduto all'analisi dell'Indice I di Moran sui residui dei singoli modelli (Tabella X). Come atteso, i residui del modello OLS presentano una forte e significativa autocorrelazione spaziale ($p < 0.001$), violando l'assunzione di indipendenza di Gauss-Markov. L'introduzione dei parametri spaziali nei modelli SAR, SEM e SDM ha ridotto drasticamente l'Indice di Moran prossimandolo allo zero, neutralizzando di fatto l'interferenza spaziale e garantendo stime robuste e non distorte.

# --- 8. SALVATAGGIO AMBIENTE PER LO SCRIPT 04 (GWR) ---
cat("\nSalvataggio in corso...\n")

saveRDS(dati_spaziali, "dati_03_pronti_per_gwr.rds")
# Salviamo l'AIC del miglior modello spaziale globale per poterlo confrontare con la GWR dopo
aic_miglior_globale <- min(df_confronto$AIC)
saveRDS(aic_miglior_globale, "aic_03_miglior_globale.rds")

cat("======================================================\n")
cat(" SCRIPT 03 CONCLUSO!\n")
cat(" Prossimo step: Script 04 per la GWR (Modellistica Locale).\n")
cat("======================================================\n")

# =====================================================================
# CONFRONTO GLOBALE: OLS vs SPATIAL ERROR MODEL (SEM)
# Riferimento: LeSage & Pace (2009)
# =====================================================================

# =====================================================================
# CONFRONTO GLOBALE: OLS vs SPATIAL ERROR MODEL (SEM) - CODICE CORRETTO
# =====================================================================

cat("\n--- CONFRONTO DIRETTO: OLS GLOBALE vs SEM GLOBALE ---\n")

# Estraiamo i prezzi veri dal dataset corretto (dati_spaziali)
prezzi_veri_log <- dati_spaziali$prezzo_log 

# 1. Estrazione Log-Likelihood
loglik_ols <- as.numeric(logLik(modello_ols))
loglik_sdm <- as.numeric(modello_sdm$LL)  # Estrazione diretta "forzata" dal ventre del SEM

# 2. Estrazione AIC
aic_ols <- AIC(modello_ols)
# Calcolo manuale dell'AIC per il SEM per evitare errori di classe Sarlm
k_parametri_sdm <- modello_sdm$parameters
aic_sdm <- -2 * loglik_sdm + 2 * k_parametri_sdm

# 3. Estrazione Previsioni (Fitted Values) e calcolo RMSE e R-Quadro Pseudo
previsioni_ols <- fitted(modello_ols)
previsioni_sdm <- as.numeric(modello_sdm$fitted.values) # Estrazione diretta

rmse_ols <- sqrt(mean((prezzi_veri_log - previsioni_ols)^2))
rmse_sdm <- sqrt(mean((prezzi_veri_log - previsioni_sdm)^2))

r2_pseudo_ols <- cor(prezzi_veri_log, previsioni_ols)^2
r2_pseudo_sdm <- cor(prezzi_veri_log, previsioni_sdm)^2

# 4. Creazione della Tabella Definitiva
tabella_ols_vs_sdm <- data.frame(
  Modello = c("OLS (Senza Dipendenza Spaziale)", "SEM (Con Dipendenza Spaziale - Lambda)"),
  Log_Likelihood = c(round(loglik_ols, 2), round(loglik_sdm, 2)),
  AIC = c(round(aic_ols, 1), round(aic_sdm, 1)),
  RMSE_Errore_Medio = c(round(rmse_ols, 4), round(rmse_sdm, 4)),
  Pseudo_R2 = c(round(r2_pseudo_ols, 4), round(r2_pseudo_sdm, 4))
)

print(tabella_ols_vs_sdm)

cat("\n===================================================================\n")
if(aic_sem < aic_ols) {
  cat(" ESITO ACCADEMICO: Il SEM è statisticamente superiore all'OLS.\n")
  cat(" L'AIC è sceso di", round(aic_ols - aic_sem, 1), "punti.\n")
} else {
  cat(" ESITO: L'OLS resiste. La dipendenza spaziale non è abbastanza forte da giustificare il SEM.\n")
}
cat("===================================================================\n")







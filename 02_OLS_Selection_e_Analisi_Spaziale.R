library(dplyr)    
library(MuMIn)    # Grid Search automatica 
library(car)      # test di Multicollinearità (VIF)
library(lmtest)   # test di Eteroschedasticità (bptest)
library(spdep)    # statistica spaziale (Moran, LISA, pesi spaziali)
library(sf)       # gestione della geometria spaziale
library(ggplot2)  # mappe e i grafici finali
library(sandwich) # 
library(tseries)  # test Jarque-Bera)

# CARICAMENTO DEI DATI MASTER (.rds)
dati_stat_log <- readRDS("dati_master_tabellari.rds") 
dati_puliti_log <- readRDS("dati_master_spaziali.rds")

#---------------------------------------OLS MIGLIORE---------------------------------

# Impostazione obbligatoria per MuMIn 
options(na.action = "na.fail")

# Creazione Modello Full
modello_master <- lm(prezzo_log ~ comune + superficie + bagni + piano + stato + 
                       eta_immobile + ascensore + garage + classe_energetica + 
                       locali + camere + balcone + giardino_terrazzo, 
                     data = dati_stat_log)

cat("\nAvvio della Grid Search su tutte le combinazioni possibili...\n")

# Prova tutte le combinazioni. 
# rank = "BIC" minimizza l'overfitting.
# extra = ... forza R a calcolare e salvare anche l'R-Quadro e l'AIC per i grafici!
tabella_grid_vera <- dredge(
  modello_master, 
  rank = "BIC", 
  extra = list(
    R2 = function(x) summary(x)$r.squared, 
    AIC_val = function(x) AIC(x)
  )
)

# Ho adottato un approccio in due fasi, in linea con le best-practices econometriche. Nella fase di model building (OLS), ho utilizzato il BIC per la Grid Search poiché, punendo più severamente la complessità nei grandi campioni, mi ha garantito l'estrazione della specificazione edonica vera e più parsimoniosa, evitando l'overfitting strutturale.
# Tuttavia, per la valutazione comparativa tra il modello globale (OLS) e i modelli spaziali (GWR/MGWR), sono passato all'AICc (Akaike Corrected). Come dimostrato dalla letteratura fondamentale sulla GWR (Fotheringham et al., 2002; Charlton & Fotheringham, 2009), l'AICc è lo standard metodologico obbligato per i modelli a parametri localmente variabili, in quanto integra una correzione per le dimensioni del campione e gestisce in modo corretto la traccia della matrice "hat" per il calcolo dei gradi di libertà effettivi. Usare il BIC per confrontare OLS e GWR sarebbe stato matematicamente improprio.

# Convertiamo il risultato in un comodo dataframe
df_modelli <- as.data.frame(tabella_grid_vera)

# TABELLA DEI MIGLIORI 20 MODELLI
cat("\n--- TOP 10 MODELLI (Ordinati per BIC) ---\n")

top_10_modelli <- df_modelli %>%
  arrange(BIC) %>% # I migliori hanno il BIC più basso
  head(10) %>%
  # Creiamo la colonna classifica e arrotondiamo i valori
  mutate(
    Classifica = row_number(), # Crea un contatore da 1 a 20
    AIC_val = round(AIC_val, 1),
    BIC = round(BIC, 1),
    delta = round(delta, 2),
    weight = round(weight, 3),
    R2 = round(R2, 4)
  ) %>%
  # Selezioniamo le colonne riassuntive mettendo la Classifica per prima!
  dplyr::select(Classifica, df, logLik, AIC_val, BIC, delta, weight, R2)

print(top_10_modelli)

# GRAFICO TRADE-OFF: R-QUADRO vs AIC
# Troviamo le coordinate del modello migliore assoluto (quello con BIC più basso, che è la prima riga)
modello_top <- df_modelli[1, ]

mappa_selezione_modelli <- ggplot(df_modelli, aes(x = R2, y = AIC_val)) +
  # Disegna tutti i modelli testati in grigio
  geom_point(alpha = 0.4, color = "steelblue", size = 2) +
  
  # Evidenzia il modello vincitore in rosso
  geom_point(data = modello_top, aes(x = R2, y = AIC_val), color = "red", size = 4, shape = 19) +
  
  labs(
    title = "Selezione del Modello OLS: Trade-off tra Adattamento e Complessità",
    subtitle = "Ogni punto rappresenta una combinazione di variabili generata dalla Grid Search",
    x = "Potere Esplicativo (R-Quadro)",
    y = "Criterio di Informazione (AIC)"
  ) +
  
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.10))) +
  
  theme_minimal(base_family = "serif") +
  theme(
    # Margine sotto il titolo principale (b = bottom)
    plot.title = element_text(face = "bold", size = 14, margin = margin(b = 8)), 
    
    # Margine sotto il sottotitolo (lo stacca bene dal grafico)
    plot.subtitle = element_text(margin = margin(b = 20)), 
    
    # Margine sopra l'etichetta dell'asse X (t = top)
    axis.title.x = element_text(margin = margin(t = 15)), 
    
    # Margine a destra dell'etichetta dell'asse Y (r = right)
    axis.title.y = element_text(margin = margin(r = 15)), 
    
    panel.grid.minor = element_blank()
  )

print(mappa_selezione_modelli)
# Figura X: Distribuzione dei modelli generati tramite Grid Search. Sull'asse delle ascisse è riportato il potere esplicativo ($R^2$), mentre sulle ordinate la penalizzazione per la complessità (AIC). Il punto evidenziato identifica la specificazione edonica ottima (Baseline Model), selezionata per minimizzazione del BIC.

# ESTRAZIONE DEL MODELLO VINCITORE DEFINITIVO
modello_definitivo_ols <- get.models(tabella_grid_vera, subset = 1)[[1]]

cat("\n--- SUMMARY DEL MODELLO DEFINITIVO ---\n")
print(summary(modello_definitivo_ols))

options(na.action = "na.omit")

cat("\n--- CALCOLO RSS (Residual Sum of Squares) DEL MODELLO OLS ---\n")

# 1. Estraiamo i residui (gli errori) del modello OLS definitivo
residui_ols <- residuals(modello_definitivo_ols)

# 2. Calcoliamo l'RSS: somma dei residui elevati al quadrato
rss_ols <- sum(residui_ols^2)

# In alternativa, in R puoi usare direttamente la funzione deviance()
# rss_ols <- deviance(modello_definitivo_ols) 

cat("La Somma dei Quadrati dei Residui (RSS) per l'OLS è:", round(rss_ols, 4), "\n")

# Se vuoi, puoi anche calcolare il Root Mean Square Error (RMSE) da inserire in tesi:
gradi_liberta <- df.residual(modello_definitivo_ols)
rmse_ols_standard <- sqrt(rss_ols / gradi_liberta)
cat("Il Root Mean Square Error (RMSE) dell'OLS è:", round(rmse_ols_standard, 4), "\n")

# DIAGNOSTICA OLS (INTEGRATA)

cat("\n--- 1. TEST DI NORMALITÀ DEI RESIDUI ---\n")
# Test di Jarque-Bera: verifica se gli errori si distribuiscono normalmente
jb_test <- jarque.bera.test(residuals(modello_definitivo_ols))
print(jb_test)
if(jb_test$p.value < 0.05) {
  cat("I residui non sono perfettamente normali (frequente nei large dataset immobiliari).\n")
}

cat("\n--- 2. TEST MULTICOLLINEARITA' (VIF) ---\n")
vif_risultati <- vif(modello_definitivo_ols)
print(vif_risultati)
cat("Se GVIF < 5, ogni variabile apporta informazione unica (assenza di multicollinearità).\n")

cat("\n--- 3. TEST DI ETEROSCHEDASTICITA' (Breusch-Pagan) ---\n")
bp_test <- bptest(modello_definitivo_ols)
print(bp_test)

if(bp_test$p.value < 0.05) {
  cat("\nATTENZIONE: Rilevata eteroschedasticità (p-value < 0.05).\n")
  cat("La varianza degli errori non è costante. Procedo con la correzione di White (HC3)...\n")
  
  # CORREZIONE CRUCIALE: Standard Error Robusti
  ols_robusto <- coeftest(modello_definitivo_ols, vcov = vcovHC(modello_definitivo_ols, type = "HC3"))
  
  cat("\n--- SOMMARIO OLS CON STANDARD ERROR ROBUSTI (HC3) ---\n")
  print(ols_robusto)
  cat("\nNOTA TESI: Usa questi P-value robusti per commentare la significatività delle variabili!\n")
} else {
  cat("\nOmoschedasticità confermata. I p-value del summary base sono validi.\n")
}

#La validità delle assunzioni classiche del teorema di Gauss-Markov è stata rigorosamente verificata. In particolare, il test di Jarque-Bera applicato sui residui del modello OLS ottimizzato ha restituito un p-value pari a 0.169, impedendo il rigetto dell'ipotesi nulla. Questo risultato certifica la normalità distributiva dei termini di errore, confermando che la trasformazione logaritmica della variabile dipendente (prezzo_log) e la specificazione edonica ottenuta tramite Grid Search hanno efficacemente assorbito le asimmetrie tipiche dei dati immobiliari, garantendo l'efficienza e la consistenza degli stimatori globali.

cat("\n--- 4. TEST DI MORAN SUI RESIDUI OLS ---\n")
# Estrazione delle coordinate esatte dal dataset spaziale pulito
coordinate <- st_coordinates(dati_puliti_log)
vicini_knn <- knn2nb(knearneigh(coordinate, k = 8))
matrice_pesi <- nb2listw(vicini_knn, style = "W")

test_moran <- lm.morantest(modello_definitivo_ols, matrice_pesi)
print(test_moran)

if(test_moran$p.value < 0.05) {
  cat("\nP-value < 0.05. I residui del modello OLS sono spazialmente autocorrelati.\n")
  cat("L'OLS soffre di Omitted Variable Bias spaziale: transizione ai modelli spaziali GIUSTIFICATA.\n")
}

#Al fine di garantire una comparabilità rigorosa tra gli stimatori e isolare il puro effetto dell'eterogeneità spaziale, la specificazione edonica individuata tramite Grid Search nel modello OLS globale (Baseline Specification) è stata mantenuta rigorosamente invariata per tutte le successive modellazioni (SAR, SEM, SDM, GWR e MGWR), conformemente alla prassi consolidata nella letteratura econometrica di settore (Helbich et al., 2014; Cellmer et al., 2020).

#-----------------------ANALISI LISA (Local Moran's I)---------------------------------------------

# PREPARAZIONE DATI
dati_lisa <- dati_puliti_log
coordinate_lisa <- st_coordinates(dati_lisa)

# Creazione la matrice dei pesi spaziali 
vicini_lisa <- knearneigh(coordinate_lisa, k = 8)
lista_pesi_lisa <- nb2listw(knn2nb(vicini_lisa), style = "W")

# CALCOLO DEL TEST LISA
lisa_test <- localmoran(dati_lisa$prezzo_log, lista_pesi_lisa)

# CREAZIONE DELLE CATEGORIE (HH, LL, HL, LH)
prezzo_centrato <- dati_lisa$prezzo_log - mean(dati_lisa$prezzo_log)
lag_prezzo <- lag.listw(lista_pesi_lisa, prezzo_centrato) # Il prezzo medio dei vicini

# Livello di significatività (p-value 0.05)
p_val_signif <- 0.05

# Inizializzazione colonna "Non Significativo"
dati_lisa$LISA_Cluster <- "5. Non Significativo"

# Popolazione cluster in base ai quadranti (solo per le case con p-value < 0.05)
dati_lisa$LISA_Cluster[prezzo_centrato > 0 & lag_prezzo > 0 & lisa_test[, 5] <= p_val_signif] <- "1. High-High (Hotspot)"
dati_lisa$LISA_Cluster[prezzo_centrato < 0 & lag_prezzo < 0 & lisa_test[, 5] <= p_val_signif] <- "2. Low-Low (Coldspot)"
dati_lisa$LISA_Cluster[prezzo_centrato > 0 & lag_prezzo < 0 & lisa_test[, 5] <= p_val_signif] <- "3. High-Low (Outlier)"
dati_lisa$LISA_Cluster[prezzo_centrato < 0 & lag_prezzo > 0 & lisa_test[, 5] <= p_val_signif] <- "4. Low-High (Outlier)"

# Conversione in factor per gestire i colori nel grafico
dati_lisa$LISA_Cluster <- as.factor(dati_lisa$LISA_Cluster)

# MAPPA DEGLI HOTSPOT

colori_lisa <- c(
  "1. High-High (Hotspot)" = "#d7191c",  
  "2. Low-Low (Coldspot)"  = "#2c7bb6",  
  "3. High-Low (Outlier)"  = "#fdae61",  
  "4. Low-High (Outlier)"  = "#abd9e9",  
  "5. Non Significativo"   = "#e0e0e0"   
)

mappa_lisa <- ggplot(data = dati_lisa) +
  geom_sf(aes(color = LISA_Cluster), size = 2.5, alpha = 0.9) +
  scale_color_manual(values = colori_lisa, name = "Tipologia Cluster") +
  labs(
    title = "Mappa LISA: Hotspot e Coldspot Immobiliari",
    subtitle = "Local Indicators of Spatial Association (K=8, p < 0.05)",
    x = "", y = ""
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    legend.position = "right",
    panel.grid.major = element_line(color = "white"),
    axis.text = element_blank(), # Eliminazione numeri di lat/long per una mappa più pulita
    axis.ticks = element_blank()
  )

print(mappa_lisa)

# Estrazione Outlier Spaziali

tabella_lisa <- st_drop_geometry(dati_lisa)

# Case economiche in quartieri costosi
affari_lh <- tabella_lisa %>%
  filter(LISA_Cluster == "4. Low-High (Outlier)") %>%
  dplyr::select(ID, comune, prezzo_mq_reale, superficie, stato, piano, eta_immobile) %>%
  arrange(prezzo_mq_reale) # Dalla più economica alla più cara

cat("Immobili economici in zone costose")
if(nrow(affari_lh) > 0) {
    print(head(affari_lh, 10))
} else {
  cat("Nessun cluster Low-High significativo trovato con K=8.\n")
}

# Immobili costosi in quartieri periferici o economici
sovrapprezzate_hl <- tabella_lisa %>%
  filter(LISA_Cluster == "3. High-Low (Outlier)") %>%
  dplyr::select(ID, comune, prezzo_mq_reale, superficie, stato, piano, eta_immobile) %>%
  arrange(desc(prezzo_mq_reale)) # Dalla più cara in giù

cat("\n--- LE CATTEDRALI NEL DESERTO (High-Low: Case costose in zone economiche) ---\n")
if(nrow(sovrapprezzate_hl) > 0) {
  cat("Trovate", nrow(sovrapprezzate_hl), "case con potenziale overpricing spaziale:\n")
  print(head(sovrapprezzate_hl, 10))
} else {
  cat("Nessun cluster High-Low significativo trovato con K=8.\n")
}


saveRDS(dati_lisa, "dati_02_lisa_spaziali.rds")
saveRDS(modello_definitivo_ols, "modello_02_ols_definitivo.rds")

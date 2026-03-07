# =====================================================================
# SCRIPT 05: AVM SPAZIALE CON GEOLOCALIZZAZIONE AUTOMATICA
# =====================================================================

# 1. CARICAMENTO LIBRERIE
library(sf)
library(dplyr)
library(tidygeocoder)

cat("\n--- INIZIALIZZAZIONE SIMULATORE PREZZI (AVM) ---\n")

# 2. RECUPERO PARAMETRI GLOBALI (Servono all'algoritmo per standardizzare i nuovi dati)
medie_x <- attr(matrice_std, "scaled:center")
sd_x <- attr(matrice_std, "scaled:scale")
media_y <- mean(dati_sf_metri$prezzo_log)
sd_y <- sd(dati_sf_metri$prezzo_log)

# 3. LA "RICETTA" DELL'ALGORITMO (Con Forchetta Commerciale al 68%)
simula_prezzo_mgwr <- function(indirizzo, mq, num_bagni, ha_ascensore, stato_immobile, classe_en) {
  
  cat("\nSto cercando le coordinate per:", indirizzo, "...\n")
  
  # A) Geolocalizzazione
  geo_data <- data.frame(addr = indirizzo) %>%
    geocode(addr, method = 'osm', lat = lat, long = lon, quiet = TRUE)
  
  if(is.na(geo_data$lat[1])) {
    cat("ERRORE: Indirizzo non trovato. Prova ad aggiungere la città (es. 'Rimini, Italia').\n")
    return(invisible(NULL))
  }
  
  # B) Trasformazione da GPS a Metri UTM
  nuova_casa_gps <- st_as_sf(geo_data, coords = c("lon", "lat"), crs = 4326)
  nuova_casa_metri <- st_transform(nuova_casa_gps, 32632)
  
  # C) Individuazione Micro-Quartiere
  indice_vicino <- st_nearest_feature(nuova_casa_metri, risultati_mgwr_sf)
  coef_locali <- st_drop_geometry(risultati_mgwr_sf[indice_vicino, ])
  colnames(coef_locali)[1:length(nomi_variabili_dinamici)] <- nomi_variabili_dinamici
  
  # D) Costruzione del Vettore Utente
  nomi_veri <- names(medie_x)
  valori_grezzi <- setNames(numeric(length(nomi_veri)), nomi_veri)
  
  valori_grezzi[grep("superficie", nomi_veri, ignore.case = TRUE)] <- mq
  valori_grezzi[grep("bagni", nomi_veri, ignore.case = TRUE)] <- num_bagni
  valori_grezzi[grep("ascensore", nomi_veri, ignore.case = TRUE)] <- ha_ascensore
  valori_grezzi[grep("Abitabile", nomi_veri, ignore.case = TRUE)] <- ifelse(stato_immobile == "Abitabile", 1, 0)
  valori_grezzi[grep("Ristrutturare", nomi_veri, ignore.case = TRUE)] <- ifelse(stato_immobile == "Da Ristrutturare", 1, 0)
  valori_grezzi[grep("DE", nomi_veri, ignore.case = TRUE)] <- ifelse(classe_en == "2_DE", 1, 0)
  valori_grezzi[grep("FG", nomi_veri, ignore.case = TRUE)] <- ifelse(classe_en == "3_FG", 1, 0)
  
  # E) Standardizzazione e Matematica MGWR
  valori_std <- (valori_grezzi - medie_x[nomi_veri]) / sd_x[nomi_veri]
  intercetta_locale <- as.numeric(coef_locali[["Intercetta"]])
  coef_variabili <- as.numeric(coef_locali[nomi_veri])
  
  prezzo_std_stimato <- intercetta_locale + sum(valori_std * coef_variabili, na.rm = TRUE)
  
  # F) Conversione in Euro Reali
  prezzo_log_stimato <- (prezzo_std_stimato * sd_y) + media_y
  prezzo_euro_mq <- exp(prezzo_log_stimato)
  valore_totale <- prezzo_euro_mq * mq
  
  # G) CALCOLO ACCURATEZZA (Intervallo commerciale al 68% -> 1 Deviazione Standard)
  rmse_std <- sqrt(sum(risultati_mgwr_sf$residual^2) / nrow(risultati_mgwr_sf))
  rmse_log <- rmse_std * sd_y 
  
  # ABBIAMO TOLTO IL 1.96! Ora usiamo 1 singola deviazione standard.
  log_min <- prezzo_log_stimato - (1 * rmse_log)
  log_max <- prezzo_log_stimato + (1 * rmse_log)
  
  valore_minimo <- exp(log_min) * mq
  valore_massimo <- exp(log_max) * mq
  
  # Margine di errore % 
  margine_perc <- round(((valore_massimo - valore_totale) / valore_totale) * 100, 1)
  
  # H) Stampa dello Scontrino
  cat("\n====================================================\n")
  cat("     REPORT VALUTAZIONE IMMOBILIARE (Motore MGWR)\n")
  cat("====================================================\n")
  cat("Indirizzo:", indirizzo, "\n")
  cat("Immobile: ", mq, "mq | Bagni:", num_bagni, "| Ascensore:", ha_ascensore, "\n")
  cat("Stato:    ", stato_immobile, "| Energia:", classe_en, "\n")
  cat("----------------------------------------------------\n")
  cat("PREZZO STIMATO:     €", formatC(valore_totale, format="f", digits=0, big.mark=".", decimal.mark=","), "\n")
  cat("Stima al mq:        €", formatC(prezzo_euro_mq, format="f", digits=0, big.mark=".", decimal.mark=","), "\n")
  cat("----------------------------------------------------\n")
  cat(" >> ACCURATEZZA DEL MODELLO <<\n")
  cat("Range Commerciale (68% dei casi simili):\n")
  cat("Da: €", formatC(valore_minimo, format="f", digits=0, big.mark=".", decimal.mark=","), 
      " a  €", formatC(valore_massimo, format="f", digits=0, big.mark=".", decimal.mark=","), "\n")
  cat("Margine di Errore:  +/-", margine_perc, "%\n")
  cat("====================================================\n")
  
  return(invisible(valore_totale))
}

# INPUT CARATTERISTICHE IMMOBILE

stima_test <- simula_prezzo_mgwr(
  indirizzo = "Viale Giuseppe Mazzini 16, Santarcangelo, Italia", 
  mq = 100, 
  num_bagni = 2, 
  ha_ascensore = 1, 
  stato_immobile = "Da Ristrutturare", 
  classe_en = "3_FG"
)    


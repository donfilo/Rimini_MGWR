library(readxl)    # read_excel()
library(dplyr)     # manipolazione: %>%, mutate(), filter(), case_when(), ecc.
library(sf)        # gestione spaziale: st_as_sf(), st_transform(), st_drop_geometry()
library(leaflet)   # interattiva (controllo del jittering)
library(MASS)      # trasformazione di Box-Cox (funzione boxcox())
library(ggplot2)   # grafici analitici e diagnostici
library(patchwork) # unire più grafici insieme

#------------------------------------------DATA CLEANING E SETUP INIZIALE-------------------------------------

# IMPORTAZIONE DATASET
percorso_file <- "C:/Users/flaviano/Desktop/Filippo/Dataset.xlsx"
df_raw <- read_excel(percorso_file)

# PULIZIA, TRASFORMAZIONE E CALCOLO PREZZI
df_ready <- df_raw %>%
  rename(
    prezzo_mq_asking = Prezzo_mq,
    superficie = Superficie_mq,
    indirizzo_input = Indirizzo,
    giardino_terrazzo = `Giardino/terrazzo` 
  ) %>%
  mutate(
    # Gestione Variabili Categoriche
    comune = as.factor(Comune),
    stato = factor(Stato, levels = c(1, 2, 3), labels = c("Nuovo", "Abitabile", "Da Ristrutturare")),
    
    # PULIZIA CLASSE ENERGETICA 
    Classe_ener = tolower(trimws(Classe_ener)),
    Classe_ener = ifelse(grepl("^a", Classe_ener), "a", Classe_ener), # Raggruppamento A1, A2, A3 in 'a'
    Classe_ener = as.factor(Classe_ener),
    
    locali = as.numeric(Locali),
    camere = as.numeric(Camere),
    balcone = as.numeric(Balcone),
    giardino_terrazzo = as.numeric(giardino_terrazzo),
    garage = as.numeric(Garage),
    ascensore = as.numeric(Ascensore),
    eta_immobile = 2026 - anno_di_costruzione,
    piano = as.numeric(Piano),
    bagni = as.numeric(Bagni),
    
    # Costruzione Indirizzo Completo per Geocoding
    indirizzo_completo = paste(indirizzo_input, comune, "RN, Italia", sep = ", "),
    
    # APPLICAZIONE SCONTI TRATTATIVA (Market Value Correction)
    perc_sconto = case_when(
      stato == "Nuovo" ~ 0.05,             # -5%
      stato == "Abitabile" ~ 0.07,         # -7%
      stato == "Da Ristrutturare" ~ 0.078,  # -7.8%
      TRUE ~ 0.09
    ),
    
    # Calcolo Prezzo Reale Stimato 
    prezzo_mq_reale = prezzo_mq_asking * (1 - perc_sconto),
    prezzo_totale_reale = Prezzo * (1 - perc_sconto)
  )
#------------------------------STATISTICHE DESCRITTIVE (PRE-PULIZIA: CAMPIONE GREZZO)----------------------------------------------

cat("Numero totale di immobili nel dataset grezzo:", nrow(df_ready), "\n\n")

# VARIABILI CONTINUE

var_continue_pre <- c("prezzo_totale_reale", "prezzo_mq_reale", 
                      "superficie", "eta_immobile", "locali", "camere", "bagni", "piano")

tabella_continue_pre <- data.frame(
  Variabile = var_continue_pre,
  N_Oss = sapply(df_ready[var_continue_pre], function(x) sum(!is.na(x))),
  Minimo = sapply(df_ready[var_continue_pre], function(x) round(min(x, na.rm=TRUE), 1)),
  Massimo = sapply(df_ready[var_continue_pre], function(x) round(max(x, na.rm=TRUE), 1)),
  Media = sapply(df_ready[var_continue_pre], function(x) round(mean(x, na.rm=TRUE), 1)),
  Mediana = sapply(df_ready[var_continue_pre], function(x) round(median(x, na.rm=TRUE), 1)),
  Dev_Standard = sapply(df_ready[var_continue_pre], function(x) round(sd(x, na.rm=TRUE), 1))
)

rownames(tabella_continue_pre) <- NULL
cat("--- VARIABILI CONTINUE (PRE-PULIZIA) ---\n")
print(tabella_continue_pre)
cat("\n")

# VARIABILI CATEGORIALI / BINARIE
var_categoriali_pre <- c("comune", "Quartiere", "stato", "Classe_ener", 
                         "ascensore", "balcone", "garage", "giardino_terrazzo")

crea_tabella_cat_pre <- function(var_name) {
  colonna <- df_ready[[var_name]]
  
  if (is.null(colonna)) {
    stop(paste("\n❌ ERRORE: La variabile '", var_name, "' NON ESISTE nel dataset df_ready!"))
  }
  
  t <- table(colonna, useNA = "ifany")
  p <- prop.table(t) * 100
  data.frame(
    Variabile = var_name,
    Categoria = names(t),
    Frequenza_N = as.numeric(t),
    Percentuale = paste0(round(as.numeric(p), 1), "%")
  )
}

tabella_categoriali_pre <- do.call(rbind, lapply(var_categoriali_pre, crea_tabella_cat_pre))
rownames(tabella_categoriali_pre) <- NULL

cat("--- VARIABILI CATEGORIALI E BINARIE (PRE-PULIZIA) ---\n")
print(tabella_categoriali_pre)

# ---------------------------------GEOCODING--------------------------------

#df_coords <- df_ready %>%
#geocode(address = indirizzo_completo, 
#          method = 'arcgis', # Metodo gratuito e preciso per l'Italia
#         lat = lat, 
#          long = long)

# CARICAMENTO IL CSV 
percorso_csv <- "C:/Users/flaviano/Desktop/Filippo/università/tesi/PROGETTO_R/DatasetGeocoded"
df_coords_storico <- read.csv(percorso_csv)

df_ready_con_coords <- df_ready %>%
  left_join(df_coords_storico %>% dplyr::select(ID, lat, long), by = "ID")

# Pulizia righe non geocodificate
df_clean_coords <- df_ready_con_coords %>%
  filter(!is.na(lat) & !is.na(long))

cat("Immobili localizzati:", nrow(df_clean_coords), "\n")

# JITTERING 
# Spostamento delle coordinate di ~5-10 metri per evitare sovrapposizioni esatte
set.seed(999) # Seme casuale per riproducibilità

df_gwr_final <- df_clean_coords %>%
  mutate(
    lat_jit = jitter(lat, amount = 0.0001),
    long_jit = jitter(long, amount = 0.0001)
  )

# CONTROLLO SU MAPPA
leaflet(df_gwr_final) %>%
  addTiles() %>%
  addCircleMarkers(
    lng = ~long_jit, 
    lat = ~lat_jit, 
    radius = 3, 
    color = "red", 
    stroke = FALSE, 
    fillOpacity = 0.6,
    popup = ~paste("<b>ID:</b>", ID, "<br>",
                   "<b>Indirizzo:</b>", indirizzo_input, "<br>", 
                   "<b>Prezzo Reale MQ:</b>", round(prezzo_mq_reale, 0), "€")
  )

# ---------------TRASFORMAZIONE SPAZIALE-----------------------------------------------

# Creazione oggetto spaziale (WGS84)

dati_sf <- st_as_sf(df_gwr_final, coords = c("long_jit", "lat_jit"), crs = 4326)
# In UTM 32N (Metri) per le distanze GWR
dati_sf_metri <- st_transform(dati_sf, crs = 32632)

#-----------------------------------TRASFORMAZIONE VARIABILE DIPENDENTE (BOX-COX)----------------------

# RICERCA LAMBDA OTTIMALE
# Creazione modello lineare temporaneo (sul dataset geocodificato) per stimare il lambda
modello_temp <- lm(prezzo_mq_reale ~ comune + superficie + bagni + piano + stato + 
                     eta_immobile + ascensore + garage + Classe_ener + 
                     locali + camere + balcone + giardino_terrazzo, 
                   data = dati_sf_metri)

bc <- boxcox(modello_temp, plotit = TRUE)

# Estrazione dell valore esatto di Lambda (punto più alto della curva)
lambda_ottimale <- bc$x[which.max(bc$y)]

cat("LAMBDA OTTIMALE TROVATO:", round(lambda_ottimale, 4), "\n")

# Applicazione della trasformazione al dataset
df_boxcox <- dati_sf_metri %>%
  mutate(
    prezzo_boxcox = (prezzo_mq_reale^lambda_ottimale - 1) / lambda_ottimale
  )

# GRAFICI DIAGNOSTICI DELLA TRASFORMAZIONE
# GRAFICO NORMALITA' PREZZI REALI
grafico_prima <- ggplot(df_boxcox, aes(x = prezzo_mq_reale)) +
  geom_histogram(aes(y = after_stat(density)), bins = 35, fill = "salmon", color = "white") +
  geom_density(alpha = 0.5, fill = "#FF9999", color = "darkred") +
  
  # Normale teorica
  stat_function(fun = dnorm, 
                args = list(mean = mean(df_boxcox$prezzo_mq_reale, na.rm=TRUE), 
                            sd = sd(df_boxcox$prezzo_mq_reale, na.rm=TRUE)), 
                color = "black", linetype = "dashed", linewidth = 1) +
  
  theme_minimal() +
  labs(title = "Prezzi Reali (€/mq)",
       subtitle = "Forte asimmetria positiva (Right-skewness)",
       x = "Prezzo al metro quadro (€)", y = "Densità") +
  theme(plot.title = element_text(face = "bold"))

# GRAFICO NORMALITA' PREZZI TRASFORMATI (DOPO)
grafico_dopo <- ggplot(df_boxcox, aes(x = prezzo_boxcox)) +
  geom_histogram(aes(y = after_stat(density)), bins = 35, fill = "lightblue", color = "white") +
  geom_density(alpha = 0.5, fill = "#99CCFF", color = "darkblue") +
  
  stat_function(fun = dnorm, 
                args = list(mean = mean(df_boxcox$prezzo_boxcox, na.rm=TRUE), 
                            sd = sd(df_boxcox$prezzo_boxcox, na.rm=TRUE)), 
                color = "black", linetype = "dashed", linewidth = 1) +
  
  theme_minimal() +
  labs(title = "Trasformazione Box-Cox",
       subtitle = paste("Normalizzazione raggiunta (Lambda =", round(lambda_ottimale, 3), ")"),
       x = "Prezzo Trasformato", y = "") +
  theme(plot.title = element_text(face = "bold"))

# GRAFICI UNITI COMPARATI (usando patchwork)
confronto_densita <- grafico_prima + grafico_dopo + 
  plot_annotation(
    title = "Analisi della Normalità: Effetto della Trasformazione di Box-Cox",
    theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
  )

print(confronto_densita)

# --- Q-Q PLOT ---
qqnorm(df_boxcox$prezzo_boxcox, main = "Q-Q Plot: I residui seguono la linea rossa?")
qqline(df_boxcox$prezzo_boxcox, col = "red", lwd = 2)

# Trasformazione logaritmica
df_trasformato <- dati_sf_metri %>%
  mutate(
    prezzo_boxcox = (prezzo_mq_reale^lambda_ottimale - 1) / lambda_ottimale,
    prezzo_log = log(prezzo_mq_reale)
  )

# Correlazione tra trasformazione Box-Cox e Logaritmica
correlazione_metodi <- cor(df_trasformato$prezzo_boxcox, df_trasformato$prezzo_log)
cat("\nCorrelazione lineare tra Box-Cox e Logaritmo:", round(correlazione_metodi, 4), "\n")
cat("(Un valore vicino a 1 giustifica in pieno l'uso del Logaritmo!)\n")

# GRAFICI COMPARATIVI 

# La Box-Cox esatta
plot_bc <- ggplot(df_trasformato, aes(x = prezzo_boxcox)) +
  geom_histogram(aes(y = after_stat(density)), bins = 35, fill = "lightblue", color = "white") +
  geom_density(alpha = 0.5, fill = "#99CCFF", color = "darkblue") +
  stat_function(fun = dnorm, 
                args = list(mean = mean(df_trasformato$prezzo_boxcox, na.rm=TRUE), 
                            sd = sd(df_trasformato$prezzo_boxcox, na.rm=TRUE)), 
                color = "black", linetype = "dashed", linewidth = 1) +
  theme_minimal() +
  labs(title = "Trasformazione Box-Cox esatta",
       subtitle = paste("Lambda =", round(lambda_ottimale, 3)),
       x = "Valore Box-Cox", y = "Densità")

# Il Logaritmo Naturale
plot_log <- ggplot(df_trasformato, aes(x = prezzo_log)) +
  geom_histogram(aes(y = after_stat(density)), bins = 35, fill = "#a1d99b", color = "white") +
  geom_density(alpha = 0.5, fill = "#e5f5e0", color = "darkgreen") +
  stat_function(fun = dnorm, 
                args = list(mean = mean(df_trasformato$prezzo_log, na.rm=TRUE), 
                            sd = sd(df_trasformato$prezzo_log, na.rm=TRUE)), 
                color = "black", linetype = "dashed", linewidth = 1) +
  theme_minimal() +
  labs(title = "Trasformazione Logaritmica (ln)",
       subtitle = "Forma funzionale adottata per il modello",
       x = "Logaritmo Naturale del Prezzo", y = "")

# Scatter Plot
plot_scatter <- ggplot(df_trasformato, aes(x = prezzo_boxcox, y = prezzo_log)) +
  geom_point(alpha = 0.6, color = "purple") +
  geom_smooth(method = "lm", color = "red", linetype = "dashed") +
  theme_minimal() +
  labs(title = "Confronto Diretto: Box-Cox vs Logaritmo",
       subtitle = paste("Le due trasformazioni sono quasi perfettamente correlate (R =", round(correlazione_metodi, 4), ")"),
       x = "Prezzo Trasformato (Box-Cox)",
       y = "Prezzo Trasformato (Logaritmo)")

# Uniione grafici
layout_finale <- (plot_bc | plot_log) / plot_scatter + 
  plot_annotation(
    title = "Giustificazione Metodologica: L'equivalenza delle Trasformazioni",
    theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
  )

print(layout_finale)

qqnorm(df_trasformato$prezzo_log, main = "Q-Q Plot: I residui seguono la linea rossa?")
qqline(df_trasformato$prezzo_log, col = "red", lwd = 2)

#Per correggere la forte asimmetria positiva (right-skewness) tipica delle distribuzioni dei prezzi immobiliari, è stata esplorata la trasformazione di Box-Cox. La stima di massima verosimiglianza ha individuato un parametro Lambda ottimale pari a -0.10. Essendo tale valore asintoticamente prossimo allo zero, la teoria statistica suggerisce la convergenza verso la trasformazione logaritmica naturale. A conferma di ciò, l'analisi di correlazione lineare tra la trasformazione esatta di Box-Cox (lambda = -0.10$) e il Logaritmo Naturale ha restituito un coefficiente di Pearson quasi perfetto (R = 0.9993). Alla luce di questa equivalenza matematica, si è optato per l'adozione del Logaritmo Naturale come forma funzionale definitiva. Questa scelta, standard nella letteratura edonica (Helbich et al., 2014), garantisce non solo la normalizzazione dei residui, ma permette la fondamentale interpretazione dei coefficienti di regressione in termini di semi-elasticità (variazioni percentuali del prezzo).

#---------------------------PULIZIA OUTLIER MULTIVARIATA (DISTANZA DI COOK)--------------------------------------

n_osservazioni <- nrow(df_trasformato)

# 1. Addestrameno modello FULL sul prezzo LOGARITMICO
modello_diagnostico_log <- lm(prezzo_log ~ comune + superficie + bagni + piano + stato + 
                                eta_immobile + ascensore + garage + Classe_ener + 
                                locali + camere + balcone + giardino_terrazzo, 
                              data = df_trasformato)

# 2. Calcolo Distanza di Cook per ogni immobile
cook_d <- cooks.distance(modello_diagnostico_log)

# Calcolo dinamico della soglia canonica di Bollen e Jackman
soglia_cook <- 4/n_osservazioni
cat("\nSoglia di Cook calcolata (4/n):", round(soglia_cook, 5), "\n")

outliers_index <- which(cook_d > soglia_cook)

# Pulizia outliers
immobili_eliminati <- df_trasformato[outliers_index, ]
dati_puliti_log <- df_trasformato[-outliers_index, ]

cat("Osservazioni iniziali:", nrow(df_trasformato), "\n")
cat("Numero di outlier influenti rimossi (D > 0.01):", length(outliers_index), "\n")
cat("Osservazioni nel dataset pulito definitivo:", nrow(dati_puliti_log), "\n")

# Estrazione dati senza geometria 
dati_stat_log <- st_drop_geometry(dati_puliti_log)

# Creazione tabella degli immobili più distorcenti
immobili_eliminati_df <- st_drop_geometry(immobili_eliminati)
tabella_outliers <- data.frame(
  ID_Immobile = immobili_eliminati_df$ID,            
  Comune = immobili_eliminati_df$comune,
  Superficie = immobili_eliminati_df$superficie,
  Prezzo_Reale = immobili_eliminati_df$prezzo_totale_reale,
  Distanza_Cook = round(cook_d[outliers_index], 4)
)
tabella_outliers <- tabella_outliers[order(-tabella_outliers$Distanza_Cook), ]

cat("\n--- TOP 10 IMMOBILI PIÙ DISTORCENTI ELIMINATI ---\n")
print(head(tabella_outliers, 10))

cat("\n--- GENERAZIONE GRAFICI SEPARATI: DIAGNOSTICA OUTLIER ---\n")

# Preparazione del dataframe per ggplot (rimuoviamo la geometria spaziale se presente)
dati_grafico_cook <- st_drop_geometry(df_trasformato)
dati_grafico_cook$Categoria <- ifelse(cook_d > soglia_cook, "Outlier (Eliminato)", "Dato Valido")

# GRAFICO 1: SCATTERPLOT CON LE DUE RETTE (DISTORTA VS VERA)
grafico_distorsione_log <- ggplot(dati_grafico_cook, aes(x = superficie, y = prezzo_log)) +
  geom_point(aes(color = Categoria, size = Categoria), alpha = 0.6) +
  scale_color_manual(values = c("Dato Valido" = "darkgray", "Outlier (Eliminato)" = "red")) +
  scale_size_manual(values = c("Dato Valido" = 1.5, "Outlier (Eliminato)" = 3.5)) +
  
  # Retta distorta (con gli outlier)
  geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed", linewidth = 1) +
  
  # Retta pulita (senza outlier)
  geom_smooth(data = subset(dati_grafico_cook, Categoria == "Dato Valido"), 
              method = "lm", se = FALSE, color = "blue", linewidth = 1.2) +
  
  labs(
    title = "Distorsione del Mercato Causata dagli Outlier",
    subtitle = "La linea blu mostra il vero trend senza le influenze estreme (linea rossa)",
    x = "Superficie (mq)",
    y = "Prezzo (Logaritmo)"
  ) +
  theme_minimal(base_family = "serif") +
  theme(
    plot.title = element_text(face = "bold", size = 15, margin = margin(b = 5)),
    plot.subtitle = element_text(size = 12, margin = margin(b = 15)),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
    legend.position = "bottom", # Legenda indipendente in basso
    legend.title = element_blank(),
    legend.text = element_text(size = 12, family = "serif")
  )

print(grafico_distorsione_log)

# GRAFICO 2: INDEX PLOT DELLA DISTANZA DI COOK

df_index_cook <- data.frame(
  Indice = 1:length(cook_d),
  Cook = as.numeric(cook_d),
  Categoria = dati_grafico_cook$Categoria
)

grafico_index_cook <- ggplot(df_index_cook, aes(x = Indice, y = Cook, color = Categoria)) +
  
  # Effetto Lollypop visivo
  geom_segment(aes(x = Indice, xend = Indice, y = 0, yend = Cook), alpha = 0.5, color = "gray80") +
  geom_point(aes(size = Categoria), alpha = 0.7) +
  
  # Linea di soglia rossa tratteggiata
  geom_hline(yintercept = soglia_cook, color = "red", linetype = "dashed", linewidth = 1) +
  scale_color_manual(values = c("Dato Valido" = "darkgray", "Outlier (Eliminato)" = "red")) +
  scale_size_manual(values = c("Dato Valido" = 1.5, "Outlier (Eliminato)" = 3.5)) +
  
  labs(
    title = "Distanza di Cook (Index Plot)",
    subtitle = "Identificazione delle osservazioni ad alta leva (Leverage)",
    x = "Indice dell'Immobile nel Dataset",
    y = "Valore della Distanza di Cook (D)"
  ) +
  theme_minimal(base_family = "serif") +
  theme(
    plot.title = element_text(face = "bold", size = 15, margin = margin(b = 5)),
    plot.subtitle = element_text(size = 12, margin = margin(b = 15)),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
    legend.position = "bottom", # Legenda indipendente in basso
    legend.title = element_blank(),
    legend.text = element_text(size = 12, family = "serif")
  )

print(grafico_index_cook)

#---------CATEGORIZZAZIONE CLASSI ENERGETICHE---------------

#La selezione della combinazione ottimale di regressori, con particolare attenzione alla classificazione dell'Attestato di Prestazione Energetica (APE), è stata condotta seguendo un approccio 'Information-Theoretic' (Burnham & Anderson, 2002). Anziché ricorrere a metodi algoritmici di regolarizzazione (es. regressione LASSO) o ad approcci euristici come lo Stepwise, si è optato per una Exhaustive Grid Search su 63 combinazioni di classi energetiche. > Questa scelta metodologica trova forte giustificazione nella letteratura dell'economia immobiliare (Brounen & Kok, 2011; Fuerst et al., 2015), la quale evidenzia come l'impiego diretto delle sette classi energetiche (A-G) nei modelli edonici generi spesso distorsioni statistiche dovute all'esiguità campionaria nelle code della distribuzione (Sparse Data Bias). L'algoritmo di selezione ha massimizzato l'Adjusted R-squared e minimizzato l'AIC, portando all'identificazione empirica di una polarizzazione binaria del mercato (Classi A-B-C vs D-E-F-G), isolando un 'Brown Discount' in perfetta aderenza con i risultati empirici internazionali sui mercati europei.

# Pulizia sul dataset STATISTICO (Per il Torneo OLS)
dati_stat_log <- dati_stat_log %>%
  mutate(
    Classe_ener = tolower(trimws(Classe_ener)),
    Classe_ener = ifelse(grepl("^a", Classe_ener), "a", Classe_ener),
    Classe_ener = as.factor(Classe_ener),
    classe_singola = factor(Classe_ener, levels = c("a", "b", "c", "d", "e", "f", "g"))
  )

# Pulizia identica sul dataset SPAZIALE (Per la futura GWR)
dati_puliti_log <- dati_puliti_log %>%
  mutate(
    Classe_ener = tolower(trimws(Classe_ener)),
    Classe_ener = ifelse(grepl("^a", Classe_ener), "a", Classe_ener),
    Classe_ener = as.factor(Classe_ener),
    classe_singola = factor(Classe_ener, levels = c("a", "b", "c", "d", "e", "f", "g"))
  )

cat("\n--- CONTEGGIO IMMOBILI PER CLASSE ORIGINALE (POST-COOK) ---\n")
print(table(dati_stat_log$Classe_ener))

# ANALISI ESPLORATIVA

colori_ape_completa <- c("a"="#00a651", "b"="#50b848", "c"="#c4d400", "d"="#fff200", "e"="#fdb913", "f"="#f37021", "g"="#ed1c24")

plot_singole <- ggplot(dati_stat_log %>% filter(!is.na(classe_singola)), aes(x = toupper(classe_singola), y = prezzo_mq_reale, fill = classe_singola)) +
  geom_boxplot(alpha = 0.8, outlier.colour = "black", outlier.size = 1.5, outlier.alpha = 0.5) +
  scale_fill_manual(values = colori_ape_completa) +
  coord_cartesian(ylim = quantile(dati_stat_log$prezzo_mq_reale, c(0.01, 0.95), na.rm=TRUE)) +
  labs(title = "Prezzi per Singola Classe Energetica (Dataset Pulito)", x = "Classe Originale (A-G)", y = "Prezzo Reale (€/mq)") +
  theme_minimal() + theme(legend.position = "none", plot.title = element_text(face="bold"))

print(plot_singole)

# MODELLO BASE (Tutte le classi, CONTROLLANDO PER IL COMUNE)
modello_base <- lm(prezzo_log ~ comune + superficie + bagni + piano + stato + eta_immobile + Classe_ener, data = dati_stat_log)

coeff_classi <- as.data.frame(summary(modello_base)$coefficients) %>%
  mutate(Variabile = rownames(.)) %>%
  filter(grepl("Classe_ener", Variabile)) %>%
  mutate(Classe = toupper(gsub("Classe_ener", "", Variabile)))

plot_salti <- ggplot(coeff_classi, aes(x = Classe, y = Estimate)) +
  
  # La linea dello zero (Baseline = Classe A) messa SOTTO per non coprire i punti
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 1) +
  
  geom_errorbar(aes(ymin = Estimate - (1.96 * `Std. Error`), 
                    ymax = Estimate + (1.96 * `Std. Error`)), 
                width = 0.2, color = "gray40", linewidth = 0.8) +
  
  # I punti (Le stime dei coefficienti)
  geom_point(size = 4, color = "darkorange") +
  
  # Testi e Titoli (Leggermente resi più formali per la tesi)
  labs(
    title = "Impatto Marginale delle Classi Energetiche sul Valore", 
    subtitle = "Benchmark (Linea Rossa) = Classe A. I punti indicano la penalizzazione sul prezzo (Log).", 
    y = "Stima del Coefficiente (Penalità)", 
    x = "Classe Energetica (Dalla B alla G)"
  ) +
  
  # Tema Accademico (Times New Roman e Margini)
  theme_minimal(base_family = "serif") + 
  theme(
    plot.title = element_text(face = "bold", size = 15, margin = margin(b = 5)),
    plot.subtitle = element_text(size = 12, margin = margin(b = 15)),
    axis.title.x = element_text(margin = margin(t = 15)),
    axis.title.y = element_text(margin = margin(r = 15)),
    # Nomi delle classi sull'asse X ben in evidenza
    axis.text.x = element_text(face = "bold", size = 12, color = "black"), 
    axis.text.y = element_text(size = 11, color = "black"),
    panel.grid.minor = element_blank()
  )

print(plot_salti)
#Come si evince dall'analisi dei coefficienti (Coefficient Plot con intervalli di confidenza al 95%), le classi energetiche intermedie (B e C) non presentano una differenza statisticamente significativa rispetto al benchmark di massima efficienza (Classe A). Questo risultato empirico dimostra che il mercato immobiliare locale non sconta le singole variazioni micro-burocratiche dell'attestato APE, ma prezza l'efficienza procedendo per 'gradoni' o macro-fasce. Pertanto, al fine di evitare la frammentazione campionaria, mitigare il rischio di overfitting e preservare i gradi di libertà del modello globale, si è ritenuto metodologicamente stringente procedere all'aggregazione (feature engineering) della variabile originaria in tre macro-categorie (Alta, Media e Bassa Efficienza). Tale approccio parsimonioso ha garantito stimatori OLS significativamente più robusti e stabili.

# GRID SEARCH (63 COMBINAZIONI) OTTIMIZZATO

classi <- c("a", "b", "c", "d", "e", "f", "g")
dividers <- expand.grid(d1=0:1, d2=0:1, d3=0:1, d4=0:1, d5=0:1, d6=0:1)
risultati_gs_completi <- list()

for (i in 1:nrow(dividers)) {
  div <- as.numeric(dividers[i, ])
  nome_gruppi <- c(toupper(classi[1]))
  mappa <- numeric(7)
  mappa[1] <- 1
  idx_gruppo <- 1
  
  for (j in 1:6) {
    if (div[j] == 1) {
      nome_gruppi <- c(nome_gruppi, toupper(classi[j+1])) 
      idx_gruppo <- idx_gruppo + 1
    } else {
      nome_gruppi[length(nome_gruppi)] <- paste0(nome_gruppi[length(nome_gruppi)], toupper(classi[j+1]))
    }
    mappa[j+1] <- idx_gruppo
  }
  
  nome_finale <- paste(nome_gruppi, collapse = "_")
  idx_match <- match(as.character(dati_stat_log$classe_singola), classi)
  dati_stat_log$cl_temp <- as.factor(mappa[idx_match])
  
  if(length(levels(droplevels(dati_stat_log$cl_temp))) < 2) next 
  
  # MODELLO AGGIORNATO CON 'comune' e 'prezzo_log'
  mod_temp <- lm(prezzo_log ~ comune + superficie + bagni + piano + stato + eta_immobile + cl_temp, data = dati_stat_log)
  coefs <- summary(mod_temp)$coefficients
  
  # Estraiamo i p-value ESCLUDENDO l'intercetta
  p_values <- coefs[-1, "Pr(>|t|)"]
  p_value_peggiore <- max(p_values)
  num_non_sig <- sum(p_values > 0.05)
  
  # Estrazione p-value specifico di statoNuovo
  pval_nuovo <- ifelse("statoNuovo" %in% rownames(coefs), coefs["statoNuovo", "Pr(>|t|)"], NA)
  
  risultati_gs_completi[[i]] <- data.frame(
    Modello = nome_finale, 
    AIC = AIC(mod_temp), 
    R2_Adj = summary(mod_temp)$adj.r.squared,
    P_Value_Max = p_value_peggiore,
    Num_Non_Sig = num_non_sig,
    P_Value_Nuovo = pval_nuovo
  )
}

df_globale <- bind_rows(risultati_gs_completi)
dati_stat_log$cl_temp <- NULL # Pulizia colonna temporanea

# FILTRAGGIO DEI MODELLI "VERDI" (con tutte variabili significative) 
df_verdi <- df_globale %>% filter(P_Value_Max <= 0.05)

cat(" TROVATI", nrow(df_verdi), "MODELLI CON TUTTE VARIABILI SIGNIFICATIVE SU 63\n")

top_verdi_aic <- df_verdi %>% arrange(AIC) %>% head(10)

cat("TOP 10 PER AIC")
print(top_verdi_aic %>% dplyr::select(Modello, AIC, R2_Adj, P_Value_Max), row.names = FALSE)

# TOP 10 PER R-QUADRO (Maggiore è meglio)
top_verdi_r2 <- df_verdi %>%
  arrange(desc(R2_Adj)) %>%
  head(10)

cat("TOP 10 PER R-QUADRO")
print(top_verdi_r2 %>% dplyr::select(Modello, R2_Adj, AIC, P_Value_Max), row.names = FALSE)


# AGGIORNAMENTO DELLA VARIABILE DEFINITIVA NEL DATASET
migliore_assoluto <- top_verdi_aic$Modello[1]
cat("\n---> IL MODELLO CHE MINIMIZZA L'AIC È:", migliore_assoluto, "<---\n")

# GENERAZIONE GRAFICO: GRID SEARCH CON FILTRO SIGNIFICATIVITÀ

# 1. Creiamo una colonna per il colore nel dataset globale
df_globale <- df_globale %>%
  mutate(Validita = ifelse(P_Value_Max <= 0.05, 
                           "Valido (Tutti i P-Value \u2264 0.05)", 
                           "Scartato (Almeno un P-Value > 0.05)"))

# 2. Isuriamo il vincitore assoluto (Il primo della tua lista top_verdi_aic)
modello_vincitore_assoluto <- top_verdi_aic[1, ]

# 3. Creiamo il Grafico Accademico
plot_gs_significativita <- ggplot(df_globale, aes(x = R2_Adj, y = AIC)) +
  
  # Disegna tutti i 63 modelli con i colori assegnati
  geom_point(aes(color = Validita), size = 3, alpha = 0.7) +
  
  # Definiamo i colori: Verde foresta per i validi, Rosso pomodoro per gli scartati
  scale_color_manual(values = c(
    "Valido (Tutti i P-Value \u2264 0.05)" = "forestgreen", 
    "Scartato (Almeno un P-Value > 0.05)" = "tomato"
  )) +
  
  # Evidenzia il modello vincitore (tra i verdi) con un pallino più grande, nero e verde acceso
  geom_point(data = modello_vincitore_assoluto, aes(x = R2_Adj, y = AIC), 
             color = "black", fill = "chartreuse", size = 5, shape = 21, stroke = 1.2) +
  
  # Etichetta di testo attaccata al vincitore
  annotate("text", x = modello_vincitore_assoluto$R2_Adj, y = modello_vincitore_assoluto$AIC, 
           label = paste("Modello Selezionato:\n", modello_vincitore_assoluto$Modello), 
           hjust = 1.15, vjust = -0.5, 
           color = "darkgreen", fontface = "bold", family = "serif", size = 4.5) +
  
  # Testi e Titoli
  labs(
    title = "Ottimizzazione Classi Energetiche: Trade-off e Significatività",
    subtitle = "Filtro a due step: selezione del minimo AIC esclusivamente tra i modelli statisticamente validi",
    x = "Potere Esplicativo (R-Quadro Adattato)",
    y = "Criterio di Informazione (AIC)",
    color = "Stato del Modello:"
  ) +
  
  # Margini espansi per non tagliare il testo a sinistra
  scale_x_continuous(expand = expansion(mult = c(0.15, 0.10))) +
  
  # Stile Tesi: Times New Roman e legenda in basso
  theme_minimal(base_family = "serif") +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 15, margin = margin(b = 5)),
    plot.subtitle = element_text(size = 12, margin = margin(b = 15)),
    axis.title.x = element_text(margin = margin(t = 15)),
    axis.title.y = element_text(margin = margin(r = 15)),
    panel.grid.minor = element_blank()
  )

print(plot_gs_significativita)

# Creazione della mappa delle classi dal vincitore
classi <- c("a", "b", "c", "d", "e", "f", "g")
gruppi_migliori <- strsplit(migliore_assoluto, "_")[[1]]
mappa_migliore <- rep(NA, 7)

for (k in 1:length(gruppi_migliori)) {
  lettere <- strsplit(tolower(gruppi_migliori[k]), "")[[1]]
  mappa_migliore[classi %in% lettere] <- paste0(k, "_", toupper(gruppi_migliori[k]))
}

# Applicazione della mappa dinamicamente al dataset STATISTICO
idx_match_stat <- match(as.character(dati_stat_log$classe_singola), classi)
dati_stat_log$classe_energetica <- as.factor(mappa_migliore[idx_match_stat])

# Applicazione della mappa dinamicamente al dataset SPAZIALE 
idx_match_spaz <- match(as.character(dati_puliti_log$classe_singola), classi)
dati_puliti_log$classe_energetica <- as.factor(mappa_migliore[idx_match_spaz])

cat("Nuovi livelli applicati a entrambi i dataset:", levels(dati_stat_log$classe_energetica), "\n")

# Addestramento del modello OLS con il nuovo raggruppamento per estrarre i dati
modello_vincitore_ener <- lm(prezzo_log ~ comune + superficie + bagni + piano + stato + 
                               eta_immobile + ascensore + garage + classe_energetica + 
                               locali + camere + balcone + giardino_terrazzo, 
                             data = dati_stat_log)

# Estrazione dinamica dei livelli della classe energetica migliore
livelli_ener <- levels(dati_stat_log$classe_energetica)
baseline_ener <- livelli_ener[1] # Il primo livello è sempre la Baseline (Zero)

colori_dinamici <- colorRampPalette(c("#31a354", "#fec44f", "#e6550d", "#de2d26"))(length(livelli_ener))
names(colori_dinamici) <- livelli_ener

# =====================================================================
# GRAFICO 1: IL DIVARIO DI PREZZO (BOXPLOT SUI DATI OSSERVATI) 
# =====================================================================
plot_box_vincitore <- ggplot(dati_stat_log, aes(x = classe_energetica, y = prezzo_log, fill = classe_energetica)) +
  geom_boxplot(alpha = 0.8, outlier.colour = "red", outlier.shape = 4) +
  scale_fill_manual(values = colori_dinamici) +
  labs(
    title = "Distribuzione dei Prezzi",
    subtitle = "Prezzi logaritmici osservati sul mercato",
    x = "Macro-Classe Energetica",
    y = "Prezzo Logaritmico (€/mq)"
  ) +
  theme_minimal(base_family = "serif") + 
  theme(
    legend.position = "none", 
    plot.title = element_text(face = "bold", margin = margin(b = 5)),
    plot.subtitle = element_text(margin = margin(b = 15)),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.text = element_text(size = 11, color = "black")
  )

# =====================================================================
# ESTRAZIONE DATI PER GRAFICO 2: IL PREMIO ENERGETICO (OLS)
# =====================================================================
coef_ener <- summary(modello_vincitore_ener)$coefficients

# Costruzione del dataframe degli effetti partendo dalla Baseline 
df_effetti <- data.frame(
  Classe = baseline_ener,
  Effetto = 0,
  Min = 0,
  Max = 0
)

# Ciclo for dinamico per estrarre i coefficienti di TUTTI gli altri gruppi 
for (i in 2:length(livelli_ener)) {
  nome_livello <- livelli_ener[i]
  nome_coef <- paste0("classe_energetica", nome_livello) 
  
  # Estrazione delle stime
  stima <- coef_ener[nome_coef, "Estimate"]
  err <- coef_ener[nome_coef, "Std. Error"]
  
  df_effetti <- rbind(df_effetti, data.frame(
    Classe = nome_livello,
    Effetto = stima,
    Min = stima - (1.96 * err),   # Intervallo di confidenza al 95%
    Max = stima + (1.96 * err)
  ))
}

# =====================================================================
# GRAFICO 2: L'EFFETTO PURO (COEFFICIENT PLOT)
# =====================================================================
plot_effetto_vincitore <- ggplot(df_effetti, aes(x = Classe, y = Effetto, color = Classe)) +
  
  # Linea di base a zero (messa per prima così sta sotto)
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray30", linewidth = 1) +
  
  # Barre di errore e Punti
  geom_errorbar(aes(ymin = Min, ymax = Max), width = 0.2, linewidth = 0.8) +
  geom_point(size = 4) +
  scale_color_manual(values = colori_dinamici) +
  labs(
    title = "L'Effetto Marginale Puro (OLS)",
    subtitle = paste("Penalizzazione % rispetto al benchmark:", baseline_ener),
    x = "Macro-Classe Energetica",
    y = "Variazione Percentuale Stimata (Log)"
  ) +
  theme_minimal(base_family = "serif") + 
  theme(
    legend.position = "none", 
    plot.title = element_text(face = "bold", margin = margin(b = 5)),
    plot.subtitle = element_text(margin = margin(b = 15)),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.text = element_text(size = 11, color = "black", face = "bold")
  )

# =====================================================================
# UNIONE DEI GRAFICI E ESPORTAZIONE
# =====================================================================
grafico_finale_energia <- plot_box_vincitore + plot_effetto_vincitore + 
  plot_annotation(
    title = "Dinamica del Mercato: Prezzi Osservati vs Effetto Ceteris Paribus",
    theme = theme(
      plot.title = element_text(family = "serif", size = 16, face = "bold", hjust = 0.5, margin = margin(b = 20))
    )
  )

print(grafico_finale_energia)

# SALVATAGGIO FILE .RDS

# Salva il dataset statistico (senza geometria)
saveRDS(dati_stat_log, "dati_master_tabellari.rds")

# Salva il dataset spaziale (con geometria e proiezione UTM)
saveRDS(dati_puliti_log, "dati_master_spaziali.rds")

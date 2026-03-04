# 🏠 Analisi Econometrica Spaziale e AVM del Mercato Immobiliare (Rimini)

Questo repository contiene l'intera pipeline analitica in R sviluppata per la mia **Tesi di Laurea Magistrale**. Il progetto riguarda la modellazione edonica e la stima dei valori immobiliari (Automated Valuation Model - AVM) applicata al mercato della provincia di Rimini. 

Il progetto confronta rigorosamente modelli as-spaziali (OLS), modelli di econometria spaziale globale (SAR, SEM, SDM) e modelli a coefficienti localmente variabili (GWR e MGWR) per risolvere l'eterogeneità e la dipendenza spaziale tipica dei dati immobiliari.

## 📊 Panoramica del Progetto
Il codice è strutturato come una "Road Map" econometrica:
1. **Data Prep & Feature Engineering:** Geocoding, trasformazione Box-Cox (convergenza al logaritmo), pulizia degli outlier tramite Distanza di Cook e ottimizzazione delle Classi Energetiche (APE) tramite Exhaustive Grid Search.
2. **Modellazione Globale:** Selezione del modello Baseline OLS minimizzando il BIC e diagnostica spaziale dei residui (Test di Moran e mappe LISA).
3. **Econometria Spaziale:** Addestramento e confronto tramite Likelihood Ratio Test dei modelli SAR, SEM e SDM.
4. **Modellistica Locale (GWR & MGWR):** Implementazione della *Geographically Weighted Regression* e della sua evoluzione multiscala (*MGWR*), permettendo ai predittori di agire su scale spaziali (bandwidth) differenti.
5. **AVM (Automated Valuation Model):** Un simulatore finale che accetta in input un indirizzo e le caratteristiche di un immobile, lo geolocalizza in tempo reale e restituisce una stima commerciale basata sui coefficienti spaziali del micro-quartiere.

## 📂 Struttura del Codice
- `01_Data_Prep.R`: Pulizia, trasformazioni spaziali (UTM 32N), Box-Cox, e selezione ottima delle variabili.
- `02_OLS_Selection_e_Analisi_Spaziale.R`: Torneo OLS (Grid Search), diagnostica classica (VIF, Breusch-Pagan) e Local Moran's I (Hotspot/Coldspot).
- `03_Modelli_Spaziali.R`: Modelli globali (SAR, SEM, SDM) e calcolo degli impatti.
- `04_GWR.R`: Calibrazione della GWR standard con bandwidth adattiva.
- `05_MGWR.R`: Calibrazione Multiscale GWR, estrazione delle bandwidth ottimali e confronto finale delle performance predittive.
- `06_AVM_Spaziale.R`: Algoritmo predittivo per la stima in real-time di nuovi immobili.

## ⚠️ Nota sui Dati (Data Privacy)
**Il dataset originale utilizzato per questo studio non è presente nel repository.** I dati sono stati raccolti, puliti e geocodificati manualmente attraverso un lungo lavoro di ricerca indipendente, costituendo proprietà intellettuale. 

Tuttavia, il codice è progettato per essere **completamente riproducibile** con qualsiasi altro dataset immobiliare. Affinché gli script funzionino correttamente, è sufficiente fornire un file Excel (`Dataset.xlsx`) strutturato con le variabili descritte nell'analisi (es. Prezzo_mq, Superficie, Indirizzo e caratteristiche strutturali dell'immobile).

## 🛠 Requisiti e Librerie
Il codice è scritto in `R`. Per eseguirlo, assicurati di avere installato i seguenti pacchetti principali:
* **Spaziali:** `sf`, `spdep`, `GWmodel`, `spatialreg`, `tidygeocoder`
* **Data Manipulation & Stats:** `dplyr`, `MASS`, `car`, `lmtest`, `MuMIn`
* **Visualizzazione:** `ggplot2`, `leaflet`, `patchwork`

## 🙏 Ringraziamenti
Un ringraziamento speciale alla mia relatrice, la **Prof.ssa Costanza Torricelli**, per il supporto, i preziosi consigli e la guida fondamentale durante la stesura di questa Tesi di Laurea Magistrale e la realizzazione di questo progetto.

## ✒️ Autore
Filippo Vincenzi https://www.linkedin.com/in/filippo-vincenzi-3320bb1a9/

# Quale delle 4 idee ORBIS: analisi comparativa empirica

> Screening quantitativo delle quattro applicazioni candidate su dati ORBIS reali, per
> scegliere quella più fattibile e più informativa per il paper panel e2tree. Prodotto
> da `code/02_compare_ideas.R`.
> Tabelle: `output/idea_comparison.csv`, `output/idea_between_within_drivers.csv`.
> Grafici: `output/plots/idea_icc.png`, `output/plots/idea_between_within_R2.png`.

## Disegno dell'analisi

Universo: imprese europee (EU-27 + EFTA/EEA + UK), conti non consolidati, EUR,
2013-2022. Ogni idea è valutata sul proprio universo naturale: panel bilanciato
decennale, complete-case sul suo specifico target più il set di 11 predittori, senza
pre-filtri presi da un'altra idea. Per ciascuna idea misuro due dimensioni.

1. Fattibilità e necessità del metodo. Il panel e2tree serve solo quando il target è
fortemente between-dominated (ICC > 0.8). Se l'ICC è basso la varianza è già within e un
modello pooled basta, quindi il metodo non è motivato. Riporto anche la numerosità
(imprese complete-case) e l'ICC mediano dei predittori.

2. Quanto informa la decomposizione. Decomposti i dati à la Mundlak, fitto due random
forest, uno sulle medie d'impresa (between, il livello) e uno sui dati demeanati (within,
le variazioni anno-su-anno), e confronto:
- `R²_between`, `R²_within`: quanto segnale c'è nelle due dimensioni;
- divergenza dei driver: `drivers_rankcorr` (correlazione di rango tra le importance
  between e within) e `top3_overlap` (Jaccard dei primi 3). Valori bassi vogliono dire che
  livello e cambiamento sono guidati da fattori diversi, il caso in cui la doppia
  spiegazione between/within non è ridondante.

## Risultati

| Idea | Target | Imprese | **ICC target** | R²pooled | R²between | R²within | rankcorr B/W | overlap top3 |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| 1 — Redditività | ROA (P/L before tax) | 241.893 | 0.396 | 0.595 | 0.393 | 0.472 | 0.65 | 0.50 |
| 2 — Leva | Gearing (debt/equity) | 183.526 | 0.617 | 0.462 | 0.243 | 0.164 | **−0.10** | **0.20** |
| 3 — Produttività | Valore aggiunto/addetto | 243.653 | 0.574 | 0.755 | 0.516 | 0.296 | 0.57 | 0.50 |
| **4 — Solvibilità** | **Solvency ratio (equity/assets)** | **214.144** | **0.811** | 0.638 | 0.525 | 0.377 | 0.82 | 0.50 |

Gli ICC di questa tabella sono allo stadio di selezione, calcolati con lo stimatore
random-effects sugli universi grezzi (oggi non più su disco). Sul panel finale, con la
convenzione unica between-share, la solvibilità è a 0.826 (pool) / 0.839 (campione).

Driver principali (between → livello; within → cambiamento):

| Idea | Between top-3 | Within top-3 |
|---|---|---|
| 1 ROA | va/addetto, rotazione magazzino, current ratio | va/addetto, costo lavoro/ricavi, rotazione magazzino |
| 2 Leva | interest cover, liquidity, current ratio | **ROA, interest cover, dimensione** |
| 3 Produttività | dimensione, costo lavoro/ricavi, EBITDA margin | **cash-flow margin, EBITDA margin, dimensione** |
| 4 Solvibilità | current ratio, liquidity, interest cover | current ratio, liquidity, **ROA** |

## Lettura

La disponibilità di dati non è un vincolo per nessuna idea: tutte hanno 180k–244k imprese
con panel bilanciato decennale. La discriminante è la necessità del metodo.

- Solo l'Idea 4 (Solvibilità) supera la soglia ICC > 0.8 (0.811). Le altre sono troppo
  within per giustificare il metodo: ROA 0.40 (la redditività oscilla troppo
  anno-su-anno), produttività 0.57, leva 0.62. Il confronto equo a quattro vie, con
  l'Idea 3 valutata sul suo vero target (valore aggiunto per addetto), conferma la scelta
  già operata: tra gli outcome di bilancio, la solvibilità è l'unica abbastanza
  persistente da richiedere il panel e2tree.

- Per l'Idea 4 c'è il **mix giusto**: ICC alto (≈81% della varianza è between → un
  surrogato/SHAP pooled spiegherebbe quasi solo i livelli) **ma con un segnale within
  reale e sostanziale** (R²within = 0.377) che la sola spiegazione di livello
  *nasconderebbe*. È esattamente la situazione in cui la decomposizione recupera
  informazione che il modello pooled confonde. I driver lo mostrano: il *livello* di
  solvibilità è guidato da liquidità e copertura degli interessi (struttura), mentre
  le *variazioni* tirano dentro la **redditività (ROA)** — gli utili trattenuti che
  muovono il capitale anno per anno.

- Nota su un caso "interessante ma non idoneo": l'**Idea 2 (Leva)** ha la **massima
  divergenza tra driver** (rankcorr −0.10, overlap 0.20): ciò che determina il
  *livello* di leva (liquidità, copertura) è diverso da ciò che la *muove* (ROA,
  dimensione). Concettualmente è la dissociazione between/within più netta, ma con
  ICC 0.62 (< 0.8) e R² deboli (0.16–0.24): il metodo non è *necessario* e il segnale
  è scarso. Resta un buon **caso di robustezza/contrasto**, non l'applicazione
  principale.

- Idea 1 e 3 sono ben spiegate già in pooled (R²pooled 0.60 e **0.76**): aggiungono
  poco alla tesi del paper.

## Verdetto

| | Più fattibile | Risultati più interessanti per il metodo |
|---|---|---|
| **Vincitrice** | **Idea 4 — Solvibilità** | **Idea 4 — Solvibilità** |

L'Idea 4 (early-warning / solvibilità, target Solvency ratio = equity/assets) è insieme
la più fattibile e la più informativa: è l'unica che soddisfa la condizione di necessità
del metodo (ICC 0.811 > 0.8) e ha al contempo un segnale within robusto (R²within 0.377)
che un'analisi di solo livello nasconderebbe, con una narrativa decisionale doppia
(banca/regolatore sul livello strutturale, CFO/CRO sulle dinamiche guidate dalla
redditività).

Questa stessa tabella può entrare nel paper: mostra empiricamente che, tra i principali
outcome di bilancio, la solvibilità è l'unico regime ad alto ICC, così la scelta
dell'applicazione è motivata dai dati e non per assunzione. L'Idea 2 (leva) può figurare
come caso di contrasto, con alta divergenza dei driver ma ICC sotto soglia.

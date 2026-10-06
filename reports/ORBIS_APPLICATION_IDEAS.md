# Panel e2tree — Applicazione empirica su dati ORBIS (Bureau van Dijk)

> Documento di scoping per la sostituzione dell'applicazione empirica del paper
> *"An explainable decision-support framework for longitudinal machine-learning models"*
> (target journal: *Decision Support Systems*).
>
> Sostituiamo il caso "aspettativa di vita dei paesi" (WDI/WHO GHO) con un caso
> **corporate / creditizio / ESG** costruito sul panel ORBIS.
>
> Snapshot ORBIS: Dicembre 2025, copertura globale, ~53,8 mln di entità,
> panel bilanci-anno ~2000–2024 (≈20 anni). Le tabelle finanziarie esistono in
> valuta locale, EUR e USD (confronti cross-country immediati).

---

## 0. Perché ORBIS si presta a panel e2tree

Il cuore metodologico del paper è la scomposizione à la Mundlak delle feature in una
componente between (medie di lungo periodo della singola impresa, $\bar{x}_i$) e una
componente within (deviazioni anno-su-anno dalla media d'impresa, $x_{it}-\bar{x}_i$), da
cui due alberi surrogati distinti. Serve solo quando l'ICC è alto (> 0,8): quando cioè un
modello globale confonde *"chi è l'impresa"* con *"cosa le sta succedendo quest'anno"*.

I bilanci d'impresa sono un caso tipico di questa patologia:

- Inerzia strutturale (between). Dimensione (log attivo), tangibilità degli asset,
  intensità di capitale, struttura settoriale, paese e forma di finanziamento sono
  quasi-fissi nel tempo: un'impresa grande e capital-intensive lo resta per anni →
  varianza between dominante.
- Fluttuazioni annuali genuine (within). Marginalità, liquidità, leva, working capital e
  cash flow oscillano con il ciclo, con gli shock (2008-09, COVID-2020, energia-2022) e con
  le scelte gestionali → varianza within non trascurabile, ma schiacciata da quella between
  in un modello pooled.

Questo produce ICC tipicamente 0,75–0,95 sulle variabili candidate-target qui sotto:
esattamente il regime in cui SHAP/alberi surrogati pooled attribuiscono erroneamente alla
"dinamica" ciò che è "livello strutturale". È la condizione necessaria che giustifica il
metodo. (Per ogni idea va riportato l'ICC stimato come *necessity check*, in linea con
`16_panel_e2tree_necessity.R` dell'applicazione precedente.)

> **Nota di costruzione comune a tutte le idee** (riusabile in `01_build_harmonized.R`):
> - Chiave panel: `BvD ID number` × anno (da `Closing date`).
> - Filtrare **consolidamento**: usare un solo regime per impresa (preferenza `U1/U2`
>   non consolidato per le PMI; `C1/C2` per i gruppi quotati) per evitare doppi conteggi.
> - Filtrare `Number of months = 12` ed escludere closing date duplicate.
> - Panel bilanciato/quasi-bilanciato: imprese con ≥ 6–8 anni consecutivi (la within
>   richiede serie storiche per-impresa decenti).
> - Stratificare per `NACE Rev.2 main section` × `paese (ISO)` per controllare composizione.
> - Winsorizzare i ratio all'1°/99° percentile (i ratio ORBIS hanno code estreme).
> - Lavorare sulla versione **EUR** delle tabelle finanziarie per comparabilità.
> - Tabella finanziaria di riferimento: **`Industry_Global_financials_and_ratios`**
>   (~90 voci + 40 ratio già calcolati). Integrazioni da `Detailed_format_industries`
>   (breakdown debito/PP&E), `Cash_flow_Global` (rendiconto), `Legal_info` /
>   `Industry_classifications` (anagrafica, settore, età, quotazione).

---

## Idea 1 — ROA / redditività operativa (corporate finance & equity screening)

**Tema:** previsione e spiegazione della **redditività** d'impresa.

### Target (regressione, alta ICC)
**`ROA using P/L before tax (%)`** (in alternativa robustness: `EBIT margin (%)` o
`ROCE`). La redditività ha forte persistenza cross-firm (esistono imprese
strutturalmente più redditizie — qualità del business, posizionamento, settore) ma
anche un'evidente componente ciclica/idiosincratica anno-su-anno. ICC atteso ≈ 0,80–0,90.
La letteratura sulla *profitability persistence* (mean reversion lenta) documenta proprio
questa struttura.

### Predictors (8–12, tutti ORBIS `Industry_G_F&R` + `Detailed_format`)
1. `ln(Total assets)` — dimensione (proxy scala/economie).
2. `Operating revenue per employee (th)` — produttività del lavoro.
3. `Net assets turnover (x)` — efficienza nell'uso del capitale.
4. `Gross margin (%)` — potere di prezzo / struttura di costo.
5. `Tangible fixed assets / Total assets` — tangibilità / intensità di capitale.
6. `Intangible fixed assets / Total assets` — asset immateriali, R&D capitalizzata.
7. `Gearing (%)` o `Long term debt / Total assets` — leva finanziaria.
8. `Current ratio (x)` — liquidità di breve.
9. `Stock turnover (x)` e `Collection period (days)` — efficienza del working capital.
10. `Costs of employees / Operating revenue (%)` — incidenza del costo del lavoro.
11. `Operating revenue` growth (within-year) — momento di crescita.
12. `NACE Rev.2 main section` + `paese` — controlli strutturali (categoriali).

### Doppio scenario decisionale (DSS)
- **Scenario A — Between (investitore equity / analista buy-side).**
  *"Quali caratteristiche strutturali separano le imprese cronicamente più redditizie da
  quelle marginali?"* L'albero **between** spiega il livello di ROA in funzione di
  dimensione, tangibilità, produttività media, settore → **stock screening e
  benchmarking** di portafoglio: l'investitore identifica i tratti permanenti di
  "qualità" su cui costruire un universo investibile.
- **Scenario B — Within (CFO / controllo di gestione).**
  *"Dato il profilo della mia impresa, quali leve fanno muovere il ROA anno-su-anno?"*
  L'albero **within** spiega le deviazioni di ROA dalla media dell'impresa in funzione di
  variazioni di margine lordo, working capital, leva, costo del lavoro → **piano d'azione
  gestionale**: il CFO interviene sulle leve modificabili nel breve, senza essere
  fuorviato dai tratti strutturali immutabili.

---

## Idea 2 — Leva finanziaria / struttura del capitale (capital structure & rischio)

**Tema:** determinanti e spiegazione della **struttura del capitale** (trade-off vs
pecking order), con taglio di **rischio finanziario**.

### Target (regressione, alta ICC)
**`Gearing (%)`** oppure **`Long term debt / Total assets`** (leverage ratio).
La struttura del capitale è notoriamente *sticky*: le imprese hanno leverage-target
relativamente stabili (politica finanziaria, settore, tangibilità del collaterale), con
aggiustamenti parziali e lenti verso il target — esattamente alto between + within non
nullo. ICC atteso ≈ 0,82–0,92.

### Predictors (8–12)
1. `ln(Total assets)` — dimensione (accesso al credito).
2. `Tangible fixed assets / Total assets` — collaterale disponibile (trade-off theory).
3. `ROA` o `Profit margin (%)` — redditività (pecking order: più profitti → meno debito).
4. `Operating revenue` growth — opportunità di crescita.
5. `Intangible fixed assets / Total assets` — asset non collateralizzabili.
6. `Interest cover (x)` — sostenibilità del servizio del debito.
7. `Current ratio (x)` / `Liquidity ratio (x)` — buffer di liquidità.
8. `EBITDA margin (%)` — capacità di generare cassa.
9. `Net assets turnover (x)` — intensità di capitale.
10. `Depreciation & Amortization / Total assets` — scudi fiscali non da debito.
11. `firm age` (da `Date of incorporation`) — maturità / reputazione creditizia.
12. `NACE main section` + `paese` — norma settoriale di leva, sistema bancocentrico vs mercato.

### Doppio scenario decisionale (DSS)
- **Scenario A — Between (banca / regolatore creditizio / risk manager di sistema).**
  *"Che tipo di impresa è strutturalmente più indebitata, e perché?"* L'albero
  **between** mappa il livello di leva su tangibilità, dimensione, settore, paese →
  **politiche di credito e calibrazione del rischio per segmento**: la banca segmenta il
  portafoglio prestiti e il regolatore identifica cluster strutturalmente fragili.
- **Scenario B — Within (CFO / tesoreria).**
  *"Cosa spinge la mia leva sopra o sotto il mio target storico, anno per anno?"*
  L'albero **within** spiega gli scostamenti di gearing dalla media d'impresa (shock di
  redditività, crescita, copertura degli interessi) → **gestione attiva della struttura
  finanziaria** e timing del rifinanziamento.

---

## Idea 3 — Produttività del lavoro / valore aggiunto per addetto (real economy & policy)

**Tema:** **produttività** e sua dispersione, con rilevanza per policy industriale e
*misallocation* (à la Hsieh–Klenow).

### Target (regressione, alta ICC)
**`Added value / Number of employees`** (valore aggiunto per addetto; in alternativa
`Operating revenue per employee (th)` o `Profit per employee (th)`, già pronti in ORBIS).
La produttività ha enorme eterogeneità *between* persistente (gap tra imprese frontiera e
laggard che dura anni) e una componente *within* legata a ciclo, ristrutturazioni,
investimenti. ICC atteso ≈ 0,85–0,93 — il caso con la between più dominante, ottimo per
mostrare il *fallimento* di SHAP pooled.

### Predictors (8–12)
1. `ln(Total assets)` — scala.
2. `Tangible fixed assets / Number of employees` — capitale per addetto (capital deepening).
3. `Intangible fixed assets / Total assets` — capitale immateriale / tecnologia.
4. `R&D expenses / Operating revenue (%)` — intensità di innovazione (dove disponibile).
5. `Costs of employees / Operating revenue (%)` — costo del lavoro relativo.
6. `Net assets turnover (x)` — efficienza degli asset.
7. `EBITDA margin (%)` — marginalità operativa.
8. `Export revenue / Operating revenue (%)` — esposizione internazionale.
9. `firm age` — apprendimento / vintage.
10. `ln(Number of employees)` — dimensione occupazionale.
11. `Additions to Fixed Assets / Total assets` (da `Cash_flow`) — investimento corrente.
12. `NACE main section` + `regione NUTS / paese` — settore e contesto geografico.

### Doppio scenario decisionale (DSS)
- **Scenario A — Between (policy maker / agenzia di sviluppo regionale).**
  *"Quali tratti strutturali distinguono le imprese ad alta produttività da quelle a
  bassa, per settore e regione?"* L'albero **between** identifica i driver permanenti del
  divario di produttività → **disegno di policy industriale e allocazione di
  incentivi** verso i cluster e i fattori giusti (capitale immateriale, dimensione,
  apertura ai mercati).
- **Scenario B — Within (direttore operations / management aziendale).**
  *"Come faccio crescere la produttività della mia impresa rispetto al suo storico?"*
  L'albero **within** isola le leve che muovono il VA per addetto anno-su-anno
  (investimenti, mix di costo del lavoro, capacità) → **roadmap di miglioramento
  operativo** depurata dai vantaggi strutturali di partenza.

---

## Idea 4 — Early-warning di distress: Interest Coverage / liquidità (rischio & ESG-G)

**Tema:** **monitoraggio del rischio di insolvenza** in chiave continua (non binaria),
ponte verso credit risk ed ESG-governance (resilienza).

### Target (regressione, alta ICC)
**`Interest cover (x)`** (EBIT / interessi passivi) — misura continua di sostenibilità del
debito; in alternativa `Solvency ratio (%)` o un **Altman Z-score** ricalcolato dalle voci
di bilancio. La capacità di copertura ha un livello strutturale per impresa (modello di
business, leva-target) ma reagisce nettamente agli shock → alto between + within ciclico.
ICC atteso ≈ 0,78–0,88.

### Predictors (8–12)
1. `Gearing (%)` / `Long term debt / Total assets` — leva.
2. `EBITDA margin (%)` — generazione di cassa operativa.
3. `Current ratio (x)` e `Liquidity ratio (x)` — liquidità.
4. `Cash flow / Operating revenue (%)` — qualità della cassa.
5. `Net Cash from Operating Activities / Total assets` (da `Cash_flow`) — cash flow operativo.
6. `Collection period (days)` − `Credit period (days)` — squilibrio del working capital.
7. `Stock turnover (x)` — efficienza delle scorte.
8. `ROA` — redditività.
9. `ln(Total assets)` e `firm age` — dimensione e maturità (proxy di resilienza).
10. `Solvency ratio (Asset based) (%)` — patrimonializzazione.
11. `Operating revenue` growth — dinamica del fatturato.
12. `NACE main section` + `paese` — ciclicità settoriale e contesto.

> Validazione esterna del target come "early warning": collegare `Legal_info.Status` /
> `StatusHistory` (dissolved/bankruptcy) per verificare che valori within bassi
> precedano l'uscita (event-study a 1–3 anni dalla cessazione).

### Doppio scenario decisionale (DSS)
- **Scenario A — Between (banca creditrice / agenzia di rating / supervisore).**
  *"Quali imprese sono strutturalmente più fragili e vanno classificate in fasce di
  rischio?"* L'albero **between** spiega il livello di copertura/solvibilità in funzione
  di leva, dimensione, settore → **rating, pricing del credito e capital allocation
  prudenziale** per segmento.
- **Scenario B — Within (CRO / CFO / credit monitoring).**
  *"La mia copertura si sta deteriorando rispetto al mio normale: quali leve la stanno
  spingendo giù quest'anno?"* L'albero **within** segnala gli scostamenti dalla media
  d'impresa (caduta di EBITDA, allungamento dei crediti, calo di cassa) → **sistema di
  early-warning e intervento tempestivo** prima dell'insolvenza.

---

## Sintesi comparativa

| # | Target | ICC atteso | Taglio decisionale | Tabelle ORBIS chiave |
|---|--------|-----------|--------------------|----------------------|
| 1 | ROA / redditività | 0,80–0,90 | Corporate finance, equity screening | Industry G F&R, Detailed |
| 2 | Gearing / leva | 0,82–0,92 | Credito, capital structure | Industry G F&R, Detailed |
| 3 | Valore aggiunto / addetto | 0,85–0,93 | Produttività, policy industriale | Industry G F&R, Cash flow |
| 4 | Interest cover / solvibilità | 0,78–0,88 | Early-warning distress, ESG-G | Industry G F&R, Cash flow, Legal/Status |

### Esito empirico (implementato)

L'applicazione realizzata è l'Idea 4, su universo europeo (EU-27 + EFTA/EEA + UK, non solo
Italia), panel bilanciato 2013-2022. Verifica sui dati ORBIS reali:

- Interest cover scartato come target: ICC ≈ 0,35 (within-dominato), non soddisfa la
  condizione di necessità del metodo.
- Target adottato: Solvency ratio (asset-based) = equity/total assets, ICC 0,826 (pool di
  227.844 imprese) / 0,839 (campione di 800), sopra la soglia 0,8. È il ratio di outcome
  più persistente disponibile e un indicatore di distress e capitalizzazione più pulito
  dell'interest cover. `gearing` (debito/equity) è tolto dai predittori perché complemento
  quasi-meccanico del target; `icover` resta come predittore (capacità di servizio del
  debito). Paese e settore restano etichette anagrafiche, non predittori. Vedi
  `data/DATA_README.md` e `output/icc_by_var.csv`.

### Raccomandazione

Per *Decision Support Systems* la combinazione più forte è l'Idea 4 (early-warning di
distress) come applicazione principale: è la più nativamente "decision support", si aggancia
a un outcome verificabile (cessazione via `StatusHistory`), e i due alberi mappano bene sui
due decisori (banca/supervisore e CRO/CFO). L'Idea 1 (ROA) resta come secondo caso di
robustezza più pulito sul piano metodologico (ICC alto ma non estremo, predittori standard,
narrativa investitore vs CFO immediata).

Le quattro idee condividono la stessa pipeline ORBIS (stessa estrazione, stesse regole di
pulizia, stesso `panel_e2tree()`): cambia solo il vettore target, quindi sono
intercambiabili a costo quasi nullo una volta costruito il panel.

---

## Note di mappatura sul codice esistente (riuso del repository)

La pipeline del paper precedente è riutilizzabile quasi 1:1 sostituendo la sorgente dati:

| Script (DSS) | Adattamento per DSS_ORBIS |
|---|---|
| `01_build_harmonized.R` | sostituire WDI/WHO con estrazione ORBIS (BvD ID × anno) e join multi-tabella |
| `02_predict_rf.R` | invariato (random forest sul nuovo panel) |
| `03_explain_e2tree.R` | invariato (surrogato e2tree pooled) |
| `06_panel_e2tree.R` | invariato (between/within e2tree) — cuore del metodo |
| `15/16_*necessity*` | ri-stimare ICC e necessity check sul target ORBIS scelto |
| `17_shap_panel.R` | SHAP between/within sul nuovo target |

I file dati grezzi ORBIS **non vanno copiati** in `data/` (decine di GB sul volume esterno
`Crucial X6`); `01_build_harmonized.R` li legge dal path esterno e salva in
`data/model_data.rds` il solo panel modellabile, come nel design originale.

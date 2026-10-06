# Report dati — Panel e2tree su ORBIS (manuscript DSS)

> Documento di provenienza dei dati per il replication package del lavoro *"An
> explainable decision-support framework for longitudinal machine-learning models"*
> (target: *Decision Support Systems*): da dove vengono, come sono stati estratti e
> costruiti, cosa sono.
>
> Fonti: `README.md`, `data/DATA_README.md`, `code/00_extract_orbis.sh`,
> `code/01_build_orbis.R`, `code/02_compare_ideas.R`, `code/03_panel_e2tree_orbis.R`,
> `reports/ORBIS_APPLICATION_IDEAS.md`, `reports/IDEA_COMPARISON.md`,
> `output/{icc_by_var,descriptives,coverage_by_var,firms_by_country}.csv`,
> `manuscript/e2tree_panel_DSS.tex` (§ *Application to European corporate solvency*).

---

## 1. In una frase

Il lavoro usa **un panel bilanciato di 800 imprese europee osservate per 10 anni
(2013–2022) = 8.000 osservazioni impresa-anno**, estratto dal database **ORBIS di
Bureau van Dijk** (snapshot Dicembre 2025). Il target è la **solvibilità**
(equity/total assets); ci sono 11 predittori finanziari (dimensione ed età inclusi come
controlli strutturali quasi time-invariant). Paese e settore NACE accompagnano il panel
come etichette anagrafiche ma **non entrano nei modelli**. A questo si affianca un
**secondo dataset sintetico simulato** con driver noti, usato solo per validare il metodo.

---

## 2. La fonte primaria: ORBIS / Bureau van Dijk

- Cos'è: ORBIS è il database globale di bilanci d'impresa di Bureau van Dijk (gruppo
  Moody's). Lo snapshot usato è quello di Dicembre 2025, con copertura globale di circa
  53,8 milioni di entità e panel bilanci-anno dal 2000 al 2024 circa.
- Perché si presta a questo metodo: i bilanci d'impresa hanno una forte inerzia
  strutturale (dimensione, settore, paese, capitalizzazione quasi-fisse nel tempo →
  varianza *between* dominante) ma anche fluttuazioni annuali genuine (margini, liquidità,
  cash flow che oscillano col ciclo → varianza *within*). Ne risultano ICC alti, il regime
  in cui una spiegazione pooled fallisce e panel e2tree serve.
- Valuta: si usa la versione EUR delle tabelle finanziarie, così i paesi sono
  confrontabili tra loro.
- Licenza: ORBIS è dato proprietario. I file grezzi e i derivati `.rds` non sono
  ridistribuibili; del replication package è pubblicabile solo il codice, e riprodurre da
  zero richiede una subscription ORBIS valida.

**Dove risiedono i grezzi**: un volume esterno non versionato —
`/Volumes/Crucial X6/Orbis Univ Dec_Work In Progress` — file `.txt` da decine di GB
ciascuno.

Tre tabelle sorgente usate:

| Tabella ORBIS | File | Cosa fornisce |
|---|---|---|
| `Industry_Global_financials_and_ratios` (EUR) | `Industry-Global_financials_and_ratios-EUR.txt` | ~90 voci di bilancio + ~40 ratio già calcolati (il cuore) |
| `Legal_info` | `Legal_info.txt` | anagrafica: status, forma legale, data di costituzione (→ età), quotazione |
| `Industry_classifications` | `Industry_classifications.txt` | settore NACE Rev.2, sezione, codice core |

---

## 3. Come sono stati estratti — `code/00_extract_orbis.sh`

Uno script bash che **streama una sola volta** i `.txt` giganti con `awk` e scrive
estratti compatti tab-separated (`.tsv`) in `data/raw/output/`. Filtri applicati già in
estrazione:

- **Universo geografico**: EU-27 + EFTA/EEA (CH, NO, IS, LI) + UK = **32 prefissi
  ISO-2**. Il paese è codificato nei primi 2 caratteri del BvD ID
  (`substr($1,1,2)`), quindi si seleziona senza join. Russia/Ucraina/Balcani non-UE
  esclusi per comparabilità istituzionale.
- **Solo conti UNCONSOLIDATED** (codice consolidamento `U*`): il regime standard per
  imprese non-gruppo, evita doppi conteggi.
- **Anni fiscali di 12 mesi** (colonna "number of months" = 12).
- **Total assets non vuoto** e **closing year in [2013, 2022]**.

Gli indici di colonna nello script si riferiscono allo schema Dicembre 2025 (verificato
sugli header). Output: `orbis_eu_financials.tsv`, `orbis_eu_legal.tsv`,
`orbis_eu_nace.tsv`.

> Nota: questi `.tsv` sono stati **cancellati** dal repo per spazio e per licenza; si
> rigenerano con lo script solo avendo il volume esterno.

---

## 4. Come sono stati costruiti/trasformati — `code/01_build_orbis.R`

Lo script R riduce gli estratti a un panel modellabile. Ha **due entry point che
producono output bit-identico** (`set.seed(1)`):

- **FULL build**: se i `.tsv` grezzi sono presenti, li rilegge e ricostruisce da zero il
  pool, riscrivendo `orbis_eu_pool.rds`.
- **FAST build**: se i grezzi sono assenti (caso normale nel repo), riparte da
  `orbis_eu_pool.rds` già salvato → stesso `model_data.rds`.

Passaggi di costruzione:

1. **Un record per impresa-anno**: dedup tenendo la closing date più recente; filtri di
   validità (total assets > 0, turnover > 0).
2. **Predittori derivati** con significato economico: `ln_assets = log(total assets)`,
   `cf_margin = 100·cashflow/turnover`, `wc_gap = collection_days − credit_days`
   (squilibrio del working capital), `rev_per_emp`, `icover` (EBIT/interessi).
3. **Join** con anagrafica legale (→ `firm_age = year − incorp_year`) e NACE (settore).
4. **Winsorizzazione all'1%/99%** su tutti i ratio a coda pesante.
5. **Complete-case + panel bilanciato**: si tengono **solo le imprese con tutti e 10 gli
   anni** e con varianza within non nulla nel target. → **pool completo: 227.844 imprese
   / 2.278.440 impresa-anno**.
6. **Campionamento riproducibile** di **800 imprese** (`sample`, seed 1): serve a tenere
   trattabile la matrice di prossimità O(N²) di e2tree.
7. **Decomposizione di Mundlak** between/within: per ogni variabile crea colonne `_bw`
   (media di lungo periodo dell'impresa x̄ᵢ) e `_wn` (deviazione annuale xᵢₜ − x̄ᵢ).

**Perché il panel bilanciato a 10 anni?** Alza l'ICC del target sopra la soglia 0,8 (con
≥7 anni scenderebbe a ~0,795) e dà una decomposizione within pulita.

---

## 5. I file dati prodotti

| File | Righe | Peso | Cos'è |
|---|---|---|---|
| `data/model_data.rds` | 8.000 | 321 KB | **Input principale di modellazione.** Panel bilanciato 800×10, complete-case. Colonne: `bvdid`, `country`, `sector` (sezione NACE), `year`, target `solvency_asset` + 11 predittori |
| `data/harmonized_panel.rds` | 8.000 | 987 KB | Stesso panel + colonne Mundlak `_bw`/`_wn` per ogni variabile |
| `data/harmonized_cross.rds` | 800 | 50 KB | Cross-section per impresa: medie di lungo periodo + etichette strutturali (paese, NACE, forma legale, età, n. anni) |
| `data/orbis_eu_pool.rds` | 2.278.440 | 145 MB | **Pool europeo completo** (227.844 imprese × 10 anni) da cui si campiona; conservato per riproducibilità senza i grezzi e per robustezza/re-sampling |
| `data/raw/output/orbis_eu_*.tsv` | — | — | Estratti grezzi pre-filtrati; **cancellati** (spazio + licenza), rigenerabili con `00_extract_orbis.sh` |

---

## 6. Cosa sono le variabili (target + 11 predittori)

**Target** — `solvency_asset` = **Solvency ratio asset-based (%)** = shareholders'
funds / total assets. Il buffer di capitalizzazione dell'impresa: struttura chi è
ben-capitalizzato vs sotto-capitalizzato, molto persistente ma con drift annuale reale.
È un indicatore di distress/early-warning.

**Predittori** (con statistiche descrittive dal campione, `output/descriptives.csv`):

| Variabile | Significato | media | sd | mediana |
|---|---|---:|---:|---:|
| `roa_pbt` | ROA su utile ante imposte (%) — redditività | 4,08 | 8,60 | 2,68 |
| `ebitda_margin` | Margine EBITDA (%) | 8,01 | 10,14 | 6,23 |
| `cf_margin` | Cash flow / ricavi (%) | 5,83 | 10,17 | 4,36 |
| `icover` | Interest coverage (EBIT/interessi) — servizio del debito | 21,13 | 64,83 | 4,14 |
| `current_ratio` | Liquidità corrente | 2,12 | 2,65 | 1,52 |
| `liquidity_ratio` | Quick ratio | 1,38 | 1,93 | 1,00 |
| `wc_gap` | Giorni incasso − giorni pagamento (squilibrio WC) | 16,25 | 73,82 | 7,49 |
| `stock_turnover` | Rotazione magazzino | 32,87 | 74,52 | 7,86 |
| `ln_assets` | log(total assets) — dimensione | 14,45 | 1,68 | 14,30 |
| `firm_age` | Età impresa (anni) | 24,12 | 15,90 | 21 |
| `rev_per_emp` | Ricavi operativi per addetto (migliaia) | 306,35 | 566,79 | 155,80 |

Le etichette anagrafiche `country` (prefisso ISO-2 del BvD ID) e `sector` (sezione NACE)
accompagnano il panel ma **non entrano come predittori** nei modelli: il set di feature
effettivo passato a *ogni* modello (`code/03_panel_e2tree_orbis.R`) è esattamente le 11
variabili qui sopra. `country` è usato solo per contare i paesi del campione (15) e per la
composizione geografica (§8); `sector` per profilare il campione. **Copertura: 100% su
tutte le variabili** (complete-case per costruzione, `output/coverage_by_var.csv`). Nota:
`gearing` (debito/equity) fu **escluso** dai predittori perché complemento quasi-meccanico
del target.

---

## 7. Perché proprio la solvibilità — il necessity check (ICC)

Il metodo serve **solo quando l'ICC del target è alta** (varianza dominata dal between).
Valori sul panel finale (stimatore *between-share*, coerente con la definizione del paper
ICC = σ²B/(σ²B+σ²W); sono i valori che entrano nel manoscritto, in `numbers.tex` e in
`tab_icc.tex`):

| Variabile | ICC pool | ICC campione |
|---|---:|---:|
| **solvency_asset (target)** | **0,826** | **0,839** |
| ln_assets, firm_age (strutturali) | ~0,97 | ~0,97 |
| predittori (roa, margini, liquidità, efficienza…) | 0,41–0,72 | 0,39–0,75 |

> **Convenzione unica dell'ICC.** In tutta la pipeline l'ICC è la **quota between della
> varianza osservata**, σ²B/(σ²B+σ²W) con σ²B+σ²W = varianza totale — la stessa usata dal
> manoscritto (`icc_of` in `03_panel_e2tree_orbis.R` e `supplement_tables.R`) e dagli
> script di build (`01_build_orbis.R`, `02_compare_ideas.R`), e l'unica coerente con la
> decomposizione additiva su cui poggiano la Proposizione 1 e il bound (1−ICC)·ρ_W.
> `output/icc_by_var.csv` è rigenerato con questa convenzione e coincide esattamente con
> `tab_icc.tex` (solvency 0,826 / 0,839). *[Storico: fino a luglio 2026 gli script di build
> usavano un ICC a componenti di varianza random-effects, ~2 punti più basso — 0,806/0,821
> per la solvibilità — poi sostituito per uniformare la convenzione.]*

Il target è fortemente *between-dominated* mentre i predittori portano vera varianza
*within* — esattamente il regime in cui un SHAP/surrogato pooled confonde "chi è
l'impresa" con "cosa le succede quest'anno".

**Confronto equo a 4 idee** (`code/02_compare_ideas.R`, `reports/IDEA_COMPARISON.md`),
ognuna valutata sul proprio universo naturale:

| Idea | Target | Imprese | ICC target | Esito |
|---|---|---:|---:|---|
| 1 Redditività | ROA | 241.893 | 0,396 | scartata (troppo within) |
| 2 Leva | Gearing | 183.526 | 0,617 | caso di contrasto |
| 3 Produttività | VA/addetto | 243.653 | 0,574 | ben spiegata già in pooled |
| 4 Solvibilità | equity/assets | 214.144 | 0,811 | scelta |

Solo l'Idea 4 supera la soglia. **Interest cover**, l'altro candidato distress, fu
scartato: ICC ~0,35, within-dominato.

> Nota: questi ICC sono **storici**, allo stadio di selezione (`02_compare_ideas.R`),
> ciascuno sul proprio universo naturale (numerosità diverse per idea) e prodotti col
> vecchio stimatore random-effects su estratti grezzi **non più su disco** (quindi non
> rigenerabili). Il codice ora usa la convenzione unica between-share: una re-run coi
> grezzi alzerebbe leggermente questi valori senza cambiare la conclusione (solo la
> solvibilità supera 0,8). I valori definitivi del panel scelto sono quelli sopra
> (0,826 / 0,839).

---

## 8. Caratteristiche del campione — composizione geografica

Il filtro balanced complete-case fa dominare il **Sud Europa**
(`output/firms_by_country.csv`, 15 paesi nel campione):

- **IT 366, ES 224, PT 96** (Sud Europa ≈ 86% delle imprese) — poi SE 37, DE 23, GB 13,
  BG 12, GR 10, BE 7, HU 4, SI 4, FR/IE/PL/RO 1 ciascuno.

Attenzione: questo sbilanciamento non è una scelta di campionamento ma riflette la
copertura ORBIS dei conti non consolidati di PMI con reporting completo dei ratio, molto
più denso in IT/ES/PT. È un limite di generalizzabilità geografica da tenere presente.

---

## 9. Il secondo dataset: la simulazione con driver noti

La sezione *"Simulation with known drivers"* del paper **non usa dati reali** ma un
**panel sintetico** (`code/03_panel_e2tree_orbis.R`, `set.seed(123)`), per verificare che
la decomposizione recuperi la *verità* e non solo migliori un numero:

- **N = 60 unità × T = 12 periodi**.
- Variabile **A** guida solo il between (costante within), variabile **B** guida solo il
  within, **Z1, Z2** sono rumore.
- Modello: `y = 50 + 2.5·A + 3·Bw + rumore`; ICC risultante moderato (sopra soglia della
  dicotomia).
- Variante con **interazione between×within** (parametro γ) per quantificare il costo
  della violazione dell'additività.

Risultato atteso: il between tree usa solo A, il within tree solo B (nessuna
contaminazione), mentre il pooled tree confonde i ruoli pur avendo alta fidelity —
dimostrando che *la fidelity da sola non diagnostica il problema panel*.

---

## 10. Catena di riproducibilità (end-to-end)

```
ORBIS .txt grezzi (esterni, licenziati)
  │  code/00_extract_orbis.sh   (awk streaming, filtri UE + U* + 12m + 2013-22)
  ▼
data/raw/output/orbis_eu_*.tsv  (NON spediti)
  │  code/01_build_orbis.R  (FULL)
  ▼
data/orbis_eu_pool.rds  (pool bilanciato complete-case: 227.844 imprese)
  │  code/01_build_orbis.R  (FAST, ricampiona identico)  + seed 1
  ▼
data/model_data.rds + harmonized_panel.rds + harmonized_cross.rds
  │  code/03_panel_e2tree_orbis.R, 17_shap_panel.R, paper_reproduce.R
  ▼
output/repro/*  →  numbers_to_tex.R  →  manuscript/numbers.tex (macro LaTeX auto-generate)
```

Ogni numero nel paper è prodotto da un'unica pipeline scriptata: `numbers.tex` è
**generato, mai editato a mano**, così la prosa non può divergere dai calcoli. Il motore
e2tree è pinnato a un commit esatto (`agostinognasso/e2tree@644d06f`) — il CRAN non
riproduce i numeri.

---

## Sintesi dei numeri-chiave

| | Valore |
|---|---|
| Fonte | ORBIS / Bureau van Dijk, snapshot Dic-2025, tabelle EUR |
| Universo | EU-27 + EFTA/EEA + UK, conti unconsolidated, 12 mesi, 2013–2022 |
| Pool completo | 227.844 imprese / 2.278.440 impresa-anno |
| **Campione di modellazione** | **800 imprese × 10 anni = 8.000 impresa-anno, 15 paesi** |
| Target | Solvency ratio = equity/total assets, ICC 0,826 (pool) / 0,839 (campione), between-share |
| Predittori | 11 finanziari (paese e settore restano etichette, non predittori) |
| Simulazione | 60 unità × 12 periodi, sintetica, driver A/B noti |
| Licenza | proprietaria, non ridistribuibile — solo codice pubblico |

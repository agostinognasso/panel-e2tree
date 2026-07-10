#!/usr/bin/env bash
# 00 — Pre-filter the heavy ORBIS .txt exports down to the European universe.
#
# The raw ORBIS exports (Bureau van Dijk, Dec-2025 snapshot) are tens-to-hundreds
# of GB each and live on an external volume; they are NOT shipped with this
# package. This script streams them once with awk and writes compact, shippable
# tab-separated extracts to data/raw/output/, which 01_build_orbis.R then turns
# into the modelling panel.
#
# Universe: EUROPEAN firms — EU-27 + EFTA/EEA (CH, NO, IS, LI) + UK (GB), 32
# countries. The BvD ID number encodes the country in its first two characters,
# so substr($1,1,2) selects the country without needing the Contact_info table.
# (Russia, Ukraine and the non-EU Balkans are excluded for institutional
# comparability and to keep the volume tractable; add their ISO-2 codes to EU
# below to include them.)
#
# Financials are further restricted to: 12-month fiscal years (col 5),
# UNCONSOLIDATED accounts (consolidation code col 2 starting with "U", the
# standard regime for non-group firms, avoids double counting), non-empty Total
# assets (col 21), and closing year in [2013, 2022]. We read the EUR-denominated
# financials file so all countries are directly comparable.
#
# Usage: point ORBIS_ROOT at the directory holding your licensed ORBIS .txt
# exports (the folder that contains Industry_Global_financials_and_ratios/,
# Legal_info/ and Industry_classifications/), then run:
#   ORBIS_ROOT="/path/to/orbis_export" bash code/00_extract_orbis.sh
#
# Column indices below refer to the Dec-2025 schema (verified against the file
# headers). If the schema changes, re-check with:  head -1 FILE | tr '\t' '\n' | nl
set -euo pipefail

: "${ORBIS_ROOT:?set ORBIS_ROOT to the directory of your licensed ORBIS .txt exports (see usage note above)}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT="$(cd "$HERE/.." && pwd)/data/raw/output"
mkdir -p "$OUT"

# European country set (ISO-2 BvD ID prefixes)
EU="AT BE BG HR CY CZ DK EE FI FR DE GR HU IE IT LV LT LU MT NL PL PT RO SK SI ES SE IS LI NO CH GB"
Y0=2013; Y1=2022

FIN="$ORBIS_ROOT/Industry_Global_financials_and_ratios/Industry-Global_financials_and_ratios-EUR.txt"
LEG="$ORBIS_ROOT/Legal_info/Legal_info.txt"
NAC="$ORBIS_ROOT/Industry_classifications/Industry_classifications.txt"

echo "[1/3] Financials (Industry G F&R, EUR, $Y0-$Y1) -> orbis_eu_financials.tsv"
LC_ALL=C awk -F'\t' -v eu="$EU" -v y0="$Y0" -v y1="$Y1" 'BEGIN{OFS="\t";
  n=split(eu,a," "); for(i=1;i<=n;i++) E[a[i]]=1;
  print "bvdid","consol","closdate","origcur","fixed_assets","intang_fixed","tang_fixed","current_assets","stock","debtors","cash","total_assets","shfunds","noncurr_liab","lt_debt","provisions","curr_liab","employees","turnover","gross_profit","ebit","pbt","net_income","cost_employees","deprec","interest_paid","cashflow","added_value","ebitda","roa_pbt","profit_margin","gross_margin","ebitda_margin","ebit_margin","interest_cover","stock_turnover","collection_days","credit_days","current_ratio","liquidity_ratio","solvency_asset","gearing","oprev_per_emp"}
  (substr($1,1,2) in E) && $5=="12" && $2 ~ /^U/ && $21!="" && substr($4,1,4)+0>=y0 && substr($4,1,4)+0<=y1 {print $1,$2,$4,$10,$12,$13,$14,$16,$17,$18,$20,$21,$22,$25,$26,$28,$29,$37,$38,$41,$43,$47,$53,$56,$57,$58,$60,$61,$62,$65,$69,$70,$71,$72,$77,$78,$79,$80,$83,$84,$86,$88,$90}' \
  "$FIN" > "$OUT/orbis_eu_financials.tsv"

echo "[2/3] Legal info -> orbis_eu_legal.tsv"
LC_ALL=C awk -F'\t' -v eu="$EU" 'BEGIN{OFS="\t"; n=split(eu,a," "); for(i=1;i<=n;i++) E[a[i]]=1;
  print "bvdid","status","status_date","legal_form","incorp_date","category","listed"}
  (substr($1,1,2) in E) {print $1,$5,$6,$7,$10,$14,$15}' \
  "$LEG" > "$OUT/orbis_eu_legal.tsv"

echo "[3/3] Industry classifications (NACE) -> orbis_eu_nace.tsv"
LC_ALL=C awk -F'\t' -v eu="$EU" 'BEGIN{OFS="\t"; n=split(eu,a," "); for(i=1;i<=n;i++) E[a[i]]=1;
  print "bvdid","nace_section","nace_core","nace_core_text","bvd_sector"}
  (substr($1,1,2) in E) {print $1,$7,$8,$9,$26}' \
  "$NAC" > "$OUT/orbis_eu_nace.tsv"

echo "Done. Extracts written to $OUT"
wc -l "$OUT"/orbis_eu_*.tsv

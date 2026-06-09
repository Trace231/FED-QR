#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RAW_DIR="$ROOT/data/raw/household_power"
PROCESSED_DIR="$ROOT/data/processed"
mkdir -p "$RAW_DIR" "$PROCESSED_DIR"

URL="https://cdn.uci-ics-mlr-prod.aws.uci.edu/235/individual%2Bhousehold%2Belectric%2Bpower%2Bconsumption.zip"
FALLBACK_URL="https://archive.ics.uci.edu/ml/machine-learning-databases/00235/household_power_consumption.zip"
ZIP_FILE="$RAW_DIR/individual_household_power_consumption.zip"
TXT_FILE="$RAW_DIR/household_power_consumption.txt"
OUT_FILE="$PROCESSED_DIR/household_power_consumption.csv"

if [ ! -f "$ZIP_FILE" ]; then
  echo "Downloading UCI Household Power Consumption data..."
  if ! curl -fL "$URL" -o "$ZIP_FILE"; then
    echo "Primary URL failed; trying UCI legacy mirror..."
    curl -fL "$FALLBACK_URL" -o "$ZIP_FILE"
  fi
fi

if [ ! -f "$TXT_FILE" ]; then
  echo "Unzipping household power data..."
  unzip -o "$ZIP_FILE" -d "$RAW_DIR"
fi

echo "Converting semicolon text to CSV..."
awk 'BEGIN{FS=";"; OFS=","} {print $1,$2,$3,$4,$5,$6,$7,$8,$9}' "$TXT_FILE" > "$OUT_FILE"

echo "Wrote:"
echo "$OUT_FILE"
wc -l "$OUT_FILE"

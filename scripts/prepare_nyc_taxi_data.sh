#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RAW_DIR="$ROOT/data/raw/nyc_taxi"
PROCESSED_DIR="$ROOT/data/processed"
mkdir -p "$RAW_DIR" "$PROCESSED_DIR"

MONTHS="${MONTHS:-2024-01 2024-02 2024-03}"
MAX_ROWS="${MAX_ROWS:-1000000}"
MIN_ZONE_ROWS="${MIN_ZONE_ROWS:-5000}"
SEED="${SEED:-20260526}"
OUT="$PROCESSED_DIR/nyc_taxi_yellow_qr_2024q1_sample.csv"
DESIGN_OUT="$PROCESSED_DIR/nyc_taxi_yellow_qr_2024q1_design_summary.csv"

if ! command -v duckdb >/dev/null 2>&1; then
  echo "duckdb CLI is required. Install it with: brew install duckdb" >&2
  exit 1
fi

ZONE_URL="https://d37ci6vzurychx.cloudfront.net/misc/taxi_zone_lookup.csv"
ZONE_FILE="$RAW_DIR/taxi_zone_lookup.csv"
if [ ! -f "$ZONE_FILE" ]; then
  echo "Downloading taxi zone lookup..."
  curl -L "$ZONE_URL" -o "$ZONE_FILE"
fi

PARQUETS=()
for month in $MONTHS; do
  file="$RAW_DIR/yellow_tripdata_${month}.parquet"
  url="https://d37ci6vzurychx.cloudfront.net/trip-data/yellow_tripdata_${month}.parquet"
  if [ ! -f "$file" ]; then
    echo "Downloading $url..."
    curl -L "$url" -o "$file"
  fi
  PARQUETS+=("$file")
done

PARQUET_LIST="$(printf "'%s'," "${PARQUETS[@]}")"
PARQUET_LIST="[${PARQUET_LIST%,}]"

SQL_FILE="$RAW_DIR/prepare_nyc_taxi.sql"
cat > "$SQL_FILE" <<SQL
CREATE OR REPLACE TEMP TABLE clean AS
SELECT
  CAST(PULocationID AS INTEGER) AS client_id,
  z.Borough AS pickup_borough,
  z.Zone AS pickup_zone,
  EXTRACT(month FROM tpep_pickup_datetime) AS pickup_month,
  EXTRACT(hour FROM tpep_pickup_datetime) AS pickup_hour,
  EXTRACT(dow FROM tpep_pickup_datetime) AS pickup_dow,
  CAST(passenger_count AS DOUBLE) AS passenger_count,
  CAST(trip_distance AS DOUBLE) AS trip_distance,
  CAST(fare_amount AS DOUBLE) AS fare_amount,
  CAST(total_amount AS DOUBLE) AS total_amount,
  CAST(tip_amount AS DOUBLE) AS tip_amount,
  CAST(tolls_amount AS DOUBLE) AS tolls_amount,
  CAST(congestion_surcharge AS DOUBLE) AS congestion_surcharge,
  CAST(airport_fee AS DOUBLE) AS airport_fee,
  CAST(VendorID AS INTEGER) AS vendor_id,
  CAST(payment_type AS INTEGER) AS payment_type,
  CAST(RatecodeID AS INTEGER) AS rate_code,
  date_diff('second', tpep_pickup_datetime, tpep_dropoff_datetime) / 60.0 AS duration_min,
  ln(1 + CAST(total_amount AS DOUBLE)) AS log_total_amount,
  ln(1 + CAST(fare_amount AS DOUBLE)) AS log_fare_amount,
  ln(1 + CAST(trip_distance AS DOUBLE)) AS log_trip_distance
FROM read_parquet($PARQUET_LIST) AS t
LEFT JOIN read_csv_auto('$ZONE_FILE') AS z
  ON CAST(t.PULocationID AS INTEGER) = CAST(z.LocationID AS INTEGER)
WHERE
  passenger_count BETWEEN 1 AND 6
  AND trip_distance > 0.05 AND trip_distance < 80
  AND fare_amount > 0 AND fare_amount < 300
  AND total_amount > 0 AND total_amount < 400
  AND tpep_dropoff_datetime > tpep_pickup_datetime
  AND date_diff('second', tpep_pickup_datetime, tpep_dropoff_datetime) BETWEEN 60 AND 14400
  AND PULocationID IS NOT NULL
  AND z.Borough IS NOT NULL
  AND z.Borough != 'Unknown';

CREATE OR REPLACE TEMP TABLE eligible_clients AS
SELECT client_id
FROM clean
GROUP BY client_id
HAVING count(*) >= $MIN_ZONE_ROWS;

COPY (
  SELECT
    client_id,
    pickup_borough,
    pickup_zone,
    pickup_month,
    pickup_hour,
    pickup_dow,
    passenger_count,
    trip_distance,
    duration_min,
    fare_amount,
    total_amount,
    tip_amount,
    tolls_amount,
    congestion_surcharge,
    airport_fee,
    vendor_id,
    payment_type,
    rate_code,
    log_total_amount,
    log_fare_amount,
    log_trip_distance
  FROM clean
  WHERE client_id IN (SELECT client_id FROM eligible_clients)
  ORDER BY hash(client_id, pickup_hour, pickup_dow, trip_distance, total_amount, $SEED)
  LIMIT $MAX_ROWS
) TO '$OUT' (HEADER, DELIMITER ',');

COPY (
  SELECT
    count(*) AS raw_clean_rows,
    count(DISTINCT client_id) AS clean_clients,
    (SELECT count(*) FROM eligible_clients) AS eligible_clients,
    $MAX_ROWS AS requested_rows,
    $MIN_ZONE_ROWS AS min_zone_rows
  FROM clean
) TO '$DESIGN_OUT' (HEADER, DELIMITER ',');
SQL

echo "Preparing NYC Taxi modeling CSV with DuckDB..."
duckdb < "$SQL_FILE"

echo "Wrote:"
echo "$OUT"
echo "$DESIGN_OUT"

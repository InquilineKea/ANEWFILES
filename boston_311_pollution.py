"""
Boston 311 + Pollution Data Explorer
Sources:
  - 311 complaints:         https://311.boston.gov/open311/v2/requests.json
  - BERDO building emissions: resource 87521565-7f15-4b8d-a225-ac4df9e3f309
  - Community GHG emissions:  resource bd8dd4bb-867e-4ca2-b6c7-6c3bd9e6c290
  - Gov't operations GHG:     resource d2b83ed7-ce3c-448d-b20a-971ed911d9e7
  - Food establishment inspections: resource 4582bec6-2b4f-4f9e-bc55-cbaa73117f4c
  - Big Belly smart trash:    resource 15e7fa44-b9a8-42da-82e1-304e43460095

Outputs:
  311_environmental.csv       — filtered 311 env complaints
  berdo_emissions.csv         — BERDO building emissions
  boston_311_map.html         — interactive Folium map
  cross_correlation.png       — time-series + spatial correlation charts
"""

import requests
import pandas as pd
import numpy as np
from datetime import datetime, timedelta

CKAN_BASE  = "https://data.boston.gov/api/3/action/datastore_search"
OPEN311    = "https://311.boston.gov/open311/v2/requests.json"

ENVIRONMENTAL_311_TYPES = {
    "Rodent Activity", "Needle Pickup", "Needle Cleanup",
    "Litter", "Overflowing Trash Can", "Residential Trash out Illegally",
    "Illegal Dumping", "Graffiti Removal", "Encampment",
    "Dead Animal Collection", "Dead Animal Pickup",
}

def fetch_311(days_back=730, max_records=2000):
    all_rows = []
    chunk = 90
    end   = datetime.now()
    start = end - timedelta(days=days_back)
    cursor = start
    while cursor < end and len(all_rows) < max_records:
        chunk_end   = min(cursor + timedelta(days=chunk), end)
        params = {
            "start_date": cursor.strftime("%Y-%m-%dT00:00:00Z"),
            "end_date":   chunk_end.strftime("%Y-%m-%dT00:00:00Z"),
            "page_size":  100,
            "status":     "open,closed",
            "page":       1,
        }
        while len(all_rows) < max_records:
            try:
                r = requests.get(OPEN311, params=params, timeout=15)
                r.raise_for_status()
                text = r.text.strip()
                if not text:
                    break
                batch = r.json()
            except Exception:
                break
            if not batch:
                break
            all_rows.extend(batch)
            if len(batch) < 100:
                break
            params["page"] += 1
        cursor = chunk_end

    if not all_rows:
        csv_path = "311_environmental.csv"
        try:
            df = pd.read_csv(csv_path)
            print(f"  (loaded {len(df)} rows from {csv_path})")
        except FileNotFoundError:
            return pd.DataFrame()
    else:
        df = pd.DataFrame(all_rows)

    if not df.empty:
        for col in ["requested_datetime", "updated_datetime"]:
            if col in df.columns:
                df[col] = pd.to_datetime(df[col], errors="coerce", utc=True)
        for col in ["lat", "long"]:
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def fetch_ckan(resource_id, limit=5000):
    rows, offset = [], 0
    while True:
        r = requests.get(CKAN_BASE, params={
            "resource_id": resource_id,
            "limit": min(limit - len(rows), 1000),
            "offset": offset,
        }, timeout=20)
        r.raise_for_status()
        result = r.json()["result"]
        batch  = result["records"]
        rows.extend(batch)
        if not batch or len(rows) >= result["total"] or len(rows) >= limit:
            break
        offset += len(batch)
    return pd.DataFrame(rows)


if __name__ == "__main__":
    print("Run boston_311_pollution.py for full analysis.")

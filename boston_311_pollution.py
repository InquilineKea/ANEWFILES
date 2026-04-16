"""
Boston 311 + Pollution Data Explorer
Sources:
  - 311: https://311.boston.gov/open311/v2/requests.json (Open311 API)
  - BERDO building emissions: Boston data portal (resource 87521565-7f15-4b8d-a225-ac4df9e3f309)
  - TreeTect urban canopy: resource 90b33f18-fb6f-4552-8c1d-28f442f15797
  - Big Belly smart trash: resource 15e7fa44-b9a8-42da-82e1-304e43460095
"""

import requests
import pandas as pd
import json
from datetime import datetime, timedelta
from collections import Counter

CKAN_BASE = "https://data.boston.gov/api/3/action/datastore_search"
OPEN311_BASE = "https://311.boston.gov/open311/v2/requests.json"

# ── 311 ───────────────────────────────────────────────────────────────────────

ENVIRONMENTAL_311_TYPES = [
    "Rodent Activity",
    "Needle Pickup",
    "Needle Cleanup",
    "Litter",
    "Overflowing Trash Can",
    "Residential Trash out Illegally",
    "Illegal Dumping",
    "Graffiti Removal",
    "Encampment",
    "Dead Animal Collection",
]

def fetch_311(days_back=90, max_records=1000):
    """Pull recent 311 requests via Open311 API."""
    start_date = (datetime.now() - timedelta(days=days_back)).strftime("%Y-%m-%dT00:00:00Z")
    params = {
        "start_date": start_date,
        "page_size": 100,
        "status": "open,closed",
    }
    all_requests = []
    page = 1
    while len(all_requests) < max_records:
        params["page"] = page
        r = requests.get(OPEN311_BASE, params=params, timeout=15)
        r.raise_for_status()
        batch = r.json()
        if not batch:
            break
        all_requests.extend(batch)
        page += 1
        if len(batch) < 100:
            break
    return pd.DataFrame(all_requests)


def fetch_ckan(resource_id, limit=5000):
    """Pull a dataset from Boston's CKAN data portal."""
    all_records = []
    offset = 0
    while True:
        r = requests.get(CKAN_BASE, params={
            "resource_id": resource_id,
            "limit": min(limit, 1000),
            "offset": offset,
        }, timeout=15)
        r.raise_for_status()
        result = r.json()["result"]
        records = result["records"]
        all_records.extend(records)
        if len(all_records) >= result["total"] or len(records) == 0 or len(all_records) >= limit:
            break
        offset += len(records)
    return pd.DataFrame(all_records)


# ── ANALYSIS ──────────────────────────────────────────────────────────────────

def analyze_311(df):
    print("\n=== 311 SERVICE REQUESTS ===")
    print(f"Total records fetched: {len(df)}")
    if df.empty:
        return df

    # Filter to env-related
    env_df = df[df["service_name"].isin(ENVIRONMENTAL_311_TYPES)].copy()
    print(f"Environment-related requests: {len(env_df)}")

    print("\nTop 10 request types:")
    print(df["service_name"].value_counts().head(10).to_string())

    print("\nEnvironmental request breakdown:")
    print(env_df["service_name"].value_counts().to_string())

    # Status breakdown
    print("\nStatus (open vs closed):")
    print(df["status"].value_counts().to_string())

    # Resolution time for closed env requests
    closed = env_df[env_df["status"] == "closed"].copy()
    if not closed.empty and "requested_datetime" in closed and "updated_datetime" in closed:
        closed["requested_datetime"] = pd.to_datetime(closed["requested_datetime"], errors="coerce")
        closed["updated_datetime"] = pd.to_datetime(closed["updated_datetime"], errors="coerce")
        closed["resolution_days"] = (
            closed["updated_datetime"] - closed["requested_datetime"]
        ).dt.total_seconds() / 86400
        print("\nAvg resolution time by type (days):")
        print(
            closed.groupby("service_name")["resolution_days"]
            .mean()
            .sort_values()
            .round(1)
            .to_string()
        )

    return env_df


def analyze_berdo(df):
    print("\n=== BERDO BUILDING EMISSIONS ===")
    print(f"Total buildings: {len(df)}")
    if df.empty:
        return

    print(f"Columns: {list(df.columns[:15])} ...")

    # Find the GHG/emissions column dynamically
    ghg_col = next((c for c in df.columns if "ghg" in c.lower() or "emission" in c.lower()), None)
    energy_col = next((c for c in df.columns if "energy" in c.lower() and "use" in c.lower()), None)
    prop_type_col = next((c for c in df.columns if "property" in c.lower() and "type" in c.lower()), None)

    if ghg_col:
        df[ghg_col] = pd.to_numeric(df[ghg_col], errors="coerce")
        print(f"\nTop 10 highest-emissions buildings ({ghg_col}):")
        top = df.nlargest(10, ghg_col)[
            [c for c in ["Property Name", "Address", prop_type_col, ghg_col] if c and c in df.columns]
        ]
        print(top.to_string(index=False))

        if prop_type_col and prop_type_col in df.columns:
            print(f"\nAvg emissions by property type:")
            print(
                df.groupby(prop_type_col)[ghg_col]
                .mean()
                .sort_values(ascending=False)
                .head(10)
                .round(0)
                .to_string()
            )


def analyze_big_belly(df):
    print("\n=== BIG BELLY SMART TRASH LOCATIONS ===")
    print(f"Total units: {len(df)}")
    if df.empty:
        return
    print(f"Columns: {list(df.columns)}")
    if "neighborhood" in [c.lower() for c in df.columns]:
        nb_col = next(c for c in df.columns if c.lower() == "neighborhood")
        print("\nUnits per neighborhood:")
        print(df[nb_col].value_counts().head(15).to_string())


def cross_analysis(df_311_env, df_berdo):
    """Simple cross-analysis: high-emissions zip codes vs 311 env complaints."""
    print("\n=== CROSS ANALYSIS: Emissions vs 311 Environmental Complaints ===")

    if df_311_env.empty or df_berdo is None or df_berdo.empty:
        print("Insufficient data for cross-analysis.")
        return

    zip_col = next((c for c in df_berdo.columns if "zip" in c.lower()), None)
    ghg_col = next((c for c in df_berdo.columns if "ghg" in c.lower() or "emission" in c.lower()), None)
    addr_col = next((c for c in df_311_env.columns if "address" in c.lower()), None)

    if zip_col and ghg_col:
        df_berdo[ghg_col] = pd.to_numeric(df_berdo[ghg_col], errors="coerce")
        high_emission_zips = (
            df_berdo.groupby(zip_col)[ghg_col]
            .sum()
            .sort_values(ascending=False)
            .head(5)
        )
        print("Top 5 zip codes by total building emissions (metric tons CO2e):")
        print(high_emission_zips.round(0).to_string())

    if addr_col:
        print("\n311 environmental complaints by street (top 10):")
        street_counts = (
            df_311_env[addr_col]
            .dropna()
            .str.extract(r"(?:\d+\s+)?(.+)")[0]
            .value_counts()
            .head(10)
        )
        print(street_counts.to_string())


# ── MAIN ──────────────────────────────────────────────────────────────────────

def main():
    print("Fetching Boston 311 data (last 90 days)...")
    try:
        df_311 = fetch_311(days_back=90, max_records=500)
        df_311_env = analyze_311(df_311)
    except Exception as e:
        print(f"311 fetch failed: {e}")
        df_311_env = pd.DataFrame()

    print("\nFetching BERDO building emissions (2024)...")
    try:
        # 2024 BERDO energy & emissions metrics
        df_berdo = fetch_ckan("87521565-7f15-4b8d-a225-ac4df9e3f309", limit=3000)
        analyze_berdo(df_berdo)
    except Exception as e:
        print(f"BERDO fetch failed: {e}")
        df_berdo = pd.DataFrame()

    print("\nFetching Big Belly smart trash locations...")
    try:
        df_belly = fetch_ckan("15e7fa44-b9a8-42da-82e1-304e43460095")
        analyze_big_belly(df_belly)
    except Exception as e:
        print(f"Big Belly fetch failed: {e}")

    cross_analysis(df_311_env, df_berdo)

    # Save outputs
    if not df_311_env.empty:
        df_311_env.to_csv("/home/user/ANEWFILES/311_environmental.csv", index=False)
        print("\nSaved: 311_environmental.csv")
    if not df_berdo.empty:
        df_berdo.to_csv("/home/user/ANEWFILES/berdo_emissions.csv", index=False)
        print("Saved: berdo_emissions.csv")


if __name__ == "__main__":
    main()

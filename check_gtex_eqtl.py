#!/usr/bin/env python3
"""
Query GTEx v8 API for eQTL associations of rs2802292 and rs4945816
across genes: FOXO3, SESN1, HACE1, ATG5

Tests whether the FOXO3 longevity haplotype affects the broader
stress-response neighborhood (trans/cis eQTL spread).

GTEx API v2 docs: https://gtexportal.org/api/v2/docs
singleTissueEqtl returns only significant eQTL associations (FDR ≤ 0.05)
"""

import json
import time
import urllib.request
import urllib.parse

GTEX_BASE = "https://gtexportal.org/api/v2"

# Pre-resolved from GTEx API (GRCh38/hg38, GTEx v8, Gencode v26)
VARIANT_IDS = {
    "rs2802292": "chr6_108587315_G_T_b38",   # chr6:108,587,315 (within FOXO3)
    "rs4945816": "chr6_108680839_C_T_b38",   # chr6:108,680,839 (3'end of FOXO3)
}

GENE_IDS = {
    "FOXO3": "ENSG00000118689.14",  # chr6:108,559,835-108,684,774
    "SESN1": "ENSG00000080546.13",  # chr6:108,986,437-109,094,819  (~300 kb downstream)
    "HACE1": "ENSG00000085382.11",  # chr6:104,728,093-104,859,919  (~3.8 Mb upstream)
    "ATG5":  "ENSG00000057663.14",  # chr6:106,045,423-106,325,791  (~2.5 Mb upstream)
}

GENES   = ["FOXO3", "SESN1", "HACE1", "ATG5"]
VARIANTS = ["rs2802292", "rs4945816"]


def gtex_get(path, params=None):
    url = f"{GTEX_BASE}{path}"
    if params:
        url += "?" + urllib.parse.urlencode(params)
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    try:
        with urllib.request.urlopen(req, timeout=30) as r:
            return json.loads(r.read())
    except Exception as e:
        return {"error": str(e)}


def get_eqtl_hits(variant_id, gencode_id, dataset="gtex_v8"):
    """
    Query singleTissueEqtl — returns only *significant* associations
    (those that passed the tissue-specific FDR ≤ 0.05 threshold in GTEx v8).
    Paginates automatically.
    """
    results = []
    page = 0
    while True:
        data = gtex_get("/association/singleTissueEqtl", {
            "variantId":  variant_id,
            "gencodeId":  gencode_id,
            "datasetId":  dataset,
            "page":       page,
            "itemsPerPage": 250,
        })
        if "error" in data:
            return None, data["error"]
        batch = data.get("data", [])
        results.extend(batch)
        total = data.get("paging_info", {}).get("totalNumberOfItems", 0)
        if len(results) >= total:
            break
        page += 1
        time.sleep(0.15)
    return results, None


def main():
    print("=" * 72)
    print("GTEx v8 eQTL lookup")
    print("Variants : rs2802292, rs4945816")
    print("Genes    : FOXO3, SESN1, HACE1, ATG5")
    print("=" * 72)
    print()
    print("Genomic context (GRCh38/hg38, all on chr6):")
    print("  HACE1    : 104.7 Mb  |  ATG5     : 106.0 Mb")
    print("  FOXO3    : 108.56–108.68 Mb")
    print("  rs2802292: 108.587 Mb  (within FOXO3)")
    print("  rs4945816: 108.681 Mb  (3' end / near FOXO3)")
    print("  SESN1    : 108.99–109.09 Mb  (~300 kb downstream of FOXO3)")
    print()

    results = {}
    for rsid in VARIANTS:
        vid = VARIANT_IDS[rsid]
        results[rsid] = {}
        for gene in GENES:
            gencode_id = GENE_IDS[gene]
            print(f"Querying {rsid} × {gene} ...", end=" ", flush=True)
            hits, err = get_eqtl_hits(vid, gencode_id)
            if err:
                print(f"ERROR: {err}")
                results[rsid][gene] = []
            else:
                results[rsid][gene] = hits
                print(f"{len(hits)} sig tissue(s)")
            time.sleep(0.2)

    # ── Pretty-print results ──────────────────────────────────────────────────
    print()
    print("=" * 72)
    print("SIGNIFICANT eQTL ASSOCIATIONS  (GTEx v8, FDR ≤ 0.05)")
    print("=" * 72)

    for rsid in VARIANTS:
        vid = VARIANT_IDS[rsid]
        print(f"\nVariant: {rsid}  [{vid}]")
        print(f"  {'Gene':<8}  {'# sig tissues':<16}  Top tissues (NES, p-value)")
        print(f"  {'-'*8}  {'-'*16}  {'-'*44}")
        for gene in GENES:
            hits = results[rsid][gene]
            if not hits:
                print(f"  {gene:<8}  {'0':<16}  —")
            else:
                hits_sorted = sorted(hits, key=lambda x: x.get("pValue", 1))
                top3 = hits_sorted[:3]
                top_str = ";  ".join(
                    f"{h['tissueSiteDetailId']} (β={h['nes']:+.3f}, p={h['pValue']:.2e})"
                    for h in top3
                )
                print(f"  {gene:<8}  {len(hits):<16}  {top_str}")

    # ── Summary matrix ────────────────────────────────────────────────────────
    print()
    print("=" * 72)
    print("SUMMARY MATRIX  (number of significant tissues)")
    print("=" * 72)
    print(f"  {'':12}  {'FOXO3':>8}  {'SESN1':>8}  {'HACE1':>8}  {'ATG5':>8}")
    for rsid in VARIANTS:
        row = "  {:<12}".format(rsid)
        for gene in GENES:
            n = len(results[rsid][gene])
            row += f"  {n:>8}"
        print(row)

    # ── Interpretation ────────────────────────────────────────────────────────
    print()
    print("=" * 72)
    print("INTERPRETATION")
    print("=" * 72)

    neighborhood_hits_rs2 = {
        g: results["rs2802292"][g]
        for g in ["SESN1", "HACE1", "ATG5"]
        if results["rs2802292"][g]
    }
    foxo3_rs2 = results["rs2802292"]["FOXO3"]

    print(f"\nrs2802292 → FOXO3 : {len(foxo3_rs2)} sig tissue(s)")

    if neighborhood_hits_rs2:
        print(f"rs2802292 → {', '.join(neighborhood_hits_rs2)} : SIGNIFICANT")
        print()
        print("★  POSITIVE RESULT  ★")
        print("  rs2802292 is an eQTL for FOXO3 AND neighboring genes.")
        print("  The FOXO3 longevity haplotype regulates the WHOLE")
        print("  stress-response neighborhood (SESN1 / HACE1 / ATG5).")
        print("  This is a BIGGER story than FOXO3 alone.")
        print("  Shared chromatin architecture has functional consequences.")
    else:
        print("rs2802292 → SESN1 / HACE1 / ATG5 : NO significant associations")
        print()
        print("◆  NEGATIVE RESULT  ◆")
        print("  rs2802292 shows NO significant eQTL effect on SESN1, HACE1, or ATG5.")
        print("  The FOXO3 haplotype effect appears relatively CONTAINED.")
        print("  Despite shared chromatin architecture, functional impact is")
        print("  likely gene-specific (FOXO3-centric).")

    # ── Save JSON ─────────────────────────────────────────────────────────────
    output = {
        "query": {
            "variants": VARIANT_IDS,
            "genes":    GENE_IDS,
            "dataset":  "gtex_v8",
            "endpoint": "singleTissueEqtl",
            "note":     "Returns only significant eQTLs (FDR ≤ 0.05 per tissue)"
        },
        "results": results,
        "summary": {
            rsid: {gene: len(results[rsid][gene]) for gene in GENES}
            for rsid in VARIANTS
        }
    }
    out_path = "/home/user/ANEWFILES/gtex_eqtl_results.json"
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2)
    print(f"\nFull results saved → {out_path}")


if __name__ == "__main__":
    main()

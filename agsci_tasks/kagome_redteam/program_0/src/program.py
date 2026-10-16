"""Starting point: one parameter set where the page's Weyl detector is known to miss nodes.

solution() returns a list (at most 12) of dicts with keys from t, breath, flux, mA, tz, dz and
optionally gap ('12' or '23'). Each distinct verified disagreement with the referee scores 1.
kagome_core.py in this folder has the detector (weyl_detect_like_js) and the scan (weyl_gap_scan).
"""


def solution():
    return [
        # Sensitive: rounding these to 2 decimals moves the case off the failure.
        {"t": 1.0, "breath": 0.08479690148362495, "flux": 306.94782185303643, "mA": 0.37176407241713605,
         "tz": -0.28788306271533215, "dz": 1.0196445630942264, "gap": "12"},
    ]

"""
Parse CELLxGENE / HsapDv `development_stage` labels into a single monotonic
numeric age axis, in years, with birth (40 post-conception weeks) at 0.0.

Why this module exists
----------------------
CELLxGENE's `development_stage` is an ontology *label*, not a number. Across
datasets the same underlying quantity shows up in ~8 mutually incompatible
surface forms ("45-year-old stage", "45-49 year-old stage", "third decade
human stage", "20 post-conception weeks", "Carnegie stage 13", "human adult
stage", ...). Regressing anything on age requires collapsing these to one
float, and doing it wrong is silent -- you get a number, it's just the wrong
number.

Two decisions worth knowing about:

1. Prenatal ages are NEGATIVE, measured from term birth. 20 pcw -> -0.385 yr.
   This keeps the axis monotonic across birth, so a model fit spanning
   gestation to adulthood does not hit a discontinuity at 0.

2. Nothing is silently dropped or silently guessed. Every parse returns a
   precision tag alongside the number. Coarse labels ("human adult stage")
   parse to a midpoint but are tagged *_VERY_IMPRECISE, because treating a
   label that spans ~40 years as a point estimate is the main way this kind
   of analysis goes wrong. Unrecognised labels return NaN, never 0.0.

Run `python3 hsapdv_age.py` to execute the self-test (no network required).
"""

from __future__ import annotations

import math
import re
from typing import NamedTuple

# Term birth in post-conception weeks. Prenatal ages are measured back from
# this point, so that age 0.0 == birth.
TERM_PCW = 40.0
WEEKS_PER_YEAR = 52.0


class AgeParse(NamedTuple):
    """Result of parsing one development_stage label."""

    years: float  # NaN when unparseable. Negative before term birth.
    precision: str  # tag describing how much to trust `years`
    span_years: float  # width of the interval `years` summarises; NaN if unknown

    @property
    def is_imprecise(self) -> bool:
        """True for labels too coarse to use as a point estimate."""
        return self.precision.endswith("VERY_IMPRECISE")

    @property
    def ok(self) -> bool:
        return not math.isnan(self.years)


NOT_PARSED = AgeParse(float("nan"), "unparsed", float("nan"))


# Carnegie stage -> (min, max) days post-fertilisation.
# Source: standard Carnegie staging tables. These are APPROXIMATE -- staging is
# morphological, not chronological, and sources disagree by a few days at the
# later stages. Only relevant if you keep embryonic samples; tagged accordingly.
_CARNEGIE_DAYS: dict[int, tuple[float, float]] = {
    1: (1, 1), 2: (2, 3), 3: (4, 5), 4: (6, 6), 5: (7, 12),
    6: (13, 15), 7: (15, 17), 8: (17, 19), 9: (19, 21), 10: (22, 23),
    11: (23, 26), 12: (26, 30), 13: (28, 32), 14: (31, 35), 15: (35, 38),
    16: (37, 42), 17: (42, 44), 18: (44, 48), 19: (48, 51), 20: (51, 53),
    21: (53, 54), 22: (54, 56), 23: (56, 60),
}

# Coarse life-stage labels with no numeric content. Ranges are the conventional
# HsapDv definitions where one exists. Anything spanning more than a decade is
# flagged VERY_IMPRECISE.
_COARSE_STAGES: dict[str, tuple[float, float]] = {
    "newborn": (0.0, 28.0 / 365.0),
    "neonate": (0.0, 28.0 / 365.0),
    "infant": (28.0 / 365.0, 1.0),
    "toddler": (1.0, 3.0),
    "child": (3.0, 13.0),
    "juvenile": (3.0, 13.0),
    "adolescent": (13.0, 19.0),
    "young adult": (19.0, 35.0),
    "prime adult": (19.0, 45.0),
    "middle aged": (45.0, 65.0),
    "late adult": (65.0, 90.0),
    "aged": (65.0, 90.0),
    "adult": (19.0, 90.0),
    "mature": (19.0, 90.0),
}

# Coarse labels wider than this (years) get the VERY_IMPRECISE tag.
_IMPRECISE_SPAN_THRESHOLD = 12.0

_ORDINAL_WORDS: dict[str, int] = {
    "first": 1, "second": 2, "third": 3, "fourth": 4, "fifth": 5,
    "sixth": 6, "seventh": 7, "eighth": 8, "ninth": 9, "tenth": 10,
    "eleventh": 11, "twelfth": 12,
}


def _pcw_to_years(pcw: float) -> float:
    """Post-conception weeks -> years relative to term birth."""
    return (pcw - TERM_PCW) / WEEKS_PER_YEAR


def _midpoint(lo: float, hi: float) -> float:
    return (lo + hi) / 2.0


def parse_development_stage(label: object) -> AgeParse:
    """Parse one `development_stage` label into years-from-birth.

    Returns an AgeParse whose `years` is NaN when the label carries no usable
    age information. Never raises on unexpected input.
    """
    if label is None:
        return NOT_PARSED
    text = str(label).strip().lower()
    if not text or text in {"na", "nan", "unknown", "none", "unspecified"}:
        return NOT_PARSED

    # Normalise unicode dashes so "45–49" parses like "45-49".
    text = text.replace("–", "-").replace("—", "-")

    # ---- ORDER MATTERS ----------------------------------------------------
    # Ranges must be tried before single years. "45-49 year-old" contains the
    # substring "49 year-old", so a single-year pattern tested first matches
    # the tail of the range and silently returns the upper bound instead of
    # the midpoint. That is a real bug this ordering prevents; the self-test
    # pins it.

    # -- year ranges: "45-49 year-old", "45 to 49 year old"
    m = re.search(r"(\d+)\s*(?:-|to)\s*(\d+)\s*[-\s]*year", text)
    if m:
        lo, hi = float(m.group(1)), float(m.group(2))
        if lo > hi:
            lo, hi = hi, lo
        return AgeParse(_midpoint(lo, hi), "range_midpoint", hi - lo)

    # -- month ages: "6-month-old", "3 to 5 month old". Before single-year so
    #    the digits are not mistaken for years.
    m = re.search(r"(\d+)\s*(?:-|to)\s*(\d+)\s*[-\s]*month", text)
    if m:
        lo, hi = float(m.group(1)) / 12.0, float(m.group(2)) / 12.0
        return AgeParse(_midpoint(lo, hi), "range_midpoint", abs(hi - lo))
    m = re.search(r"(\d+)\s*[-\s]*month[-\s]*old", text)
    if m:
        return AgeParse(float(m.group(1)) / 12.0, "exact_month", 1.0 / 12.0)

    # -- post-conception / post-fertilisation weeks
    m = re.search(
        r"(\d+)\s*(?:-|to)\s*(\d+)\s*(?:post[-\s]?conception|post[-\s]?fertilization|"
        r"post[-\s]?fertilisation|pcw|pf)\s*week",
        text,
    )
    if m:
        lo, hi = float(m.group(1)), float(m.group(2))
        return AgeParse(
            _pcw_to_years(_midpoint(lo, hi)), "pcw_range", (hi - lo) / WEEKS_PER_YEAR
        )
    m = re.search(
        r"(\d+)\s*(?:post[-\s]?conception|post[-\s]?fertilization|"
        r"post[-\s]?fertilisation|pcw|pf)\s*week",
        text,
    )
    if m:
        return AgeParse(_pcw_to_years(float(m.group(1))), "pcw", 1.0 / WEEKS_PER_YEAR)
    # "20th week post-fertilization"
    m = re.search(r"(\d+)\s*(?:st|nd|rd|th)?\s*week\s*post[-\s]?(?:conception|fertilization|fertilisation)", text)
    if m:
        return AgeParse(_pcw_to_years(float(m.group(1))), "pcw", 1.0 / WEEKS_PER_YEAR)

    # -- Carnegie stages
    m = re.search(r"carnegie\s*stage\s*(\d+)", text)
    if m:
        cs = int(m.group(1))
        if cs in _CARNEGIE_DAYS:
            lo_d, hi_d = _CARNEGIE_DAYS[cs]
            pcw = _midpoint(lo_d, hi_d) / 7.0
            return AgeParse(
                _pcw_to_years(pcw), "carnegie_APPROX", (hi_d - lo_d) / 7.0 / WEEKS_PER_YEAR
            )
        return NOT_PARSED

    # -- single year: "45-year-old", "45 year old human stage"
    m = re.search(r"(\d+)\s*[-\s]*year[-\s]*old", text)
    if m:
        return AgeParse(float(m.group(1)), "exact_year", 1.0)

    # -- decades: "third decade human stage" -> ages 20-29 -> 25
    m = re.search(r"(\w+)\s+decade", text)
    if m:
        n = _ORDINAL_WORDS.get(m.group(1))
        if n is None and m.group(1).isdigit():
            n = int(m.group(1))
        if n:
            lo = (n - 1) * 10.0
            return AgeParse(lo + 5.0, "decade_midpoint", 10.0)

    # -- bare "N weeks" in an obviously embryonic/fetal context
    m = re.search(r"(\d+)\s*week", text)
    if m and any(k in text for k in ("embryo", "fetal", "fetus", "gestation", "prenatal")):
        return AgeParse(_pcw_to_years(float(m.group(1))), "pcw", 1.0 / WEEKS_PER_YEAR)

    # -- coarse life stages. Longest key first so "young adult" beats "adult"
    #    and "late adult" is not swallowed by "adult".
    for key in sorted(_COARSE_STAGES, key=len, reverse=True):
        if key in text:
            lo, hi = _COARSE_STAGES[key]
            span = hi - lo
            tag = "coarse_%s" % key.replace(" ", "_")
            if span > _IMPRECISE_SPAN_THRESHOLD:
                tag += "_VERY_IMPRECISE"
            return AgeParse(_midpoint(lo, hi), tag, span)

    # -- generic embryonic / fetal with no number at all
    if "embryo" in text:
        return AgeParse(_pcw_to_years(6.0), "coarse_embryonic_VERY_IMPRECISE", 8.0 / WEEKS_PER_YEAR)
    if "fetal" in text or "fetus" in text:
        return AgeParse(_pcw_to_years(20.0), "coarse_fetal_VERY_IMPRECISE", 24.0 / WEEKS_PER_YEAR)

    return NOT_PARSED


# --------------------------------------------------------------------------
# Self-test. Runs offline; no network, no census, no data files.
# --------------------------------------------------------------------------

_CASES: list[tuple[str, float, str]] = [
    # (label, expected years, expected precision tag)
    ("45-year-old human stage", 45.0, "exact_year"),
    ("45-year-old stage", 45.0, "exact_year"),
    ("2-year-old human stage", 2.0, "exact_year"),
    ("89-year-old human stage", 89.0, "exact_year"),
    # The regression case: must be the MIDPOINT 47, not the upper bound 49.
    ("45-49 year-old stage", 47.0, "range_midpoint"),
    ("45-49 year-old human stage", 47.0, "range_midpoint"),
    ("65-69 year-old stage", 67.0, "range_midpoint"),
    ("30 to 39 year-old stage", 34.5, "range_midpoint"),
    ("45–49 year-old stage", 47.0, "range_midpoint"),
    # Decades.
    ("third decade human stage", 25.0, "decade_midpoint"),
    ("fourth decade human stage", 35.0, "decade_midpoint"),
    ("eighth decade human stage", 75.0, "decade_midpoint"),
    # Prenatal: negative, measured back from term.
    ("20 post-conception week human stage", (20.0 - 40.0) / 52.0, "pcw"),
    ("8 post-conception weeks human stage", (8.0 - 40.0) / 52.0, "pcw"),
    ("Carnegie stage 13", ((28 + 32) / 2 / 7 - 40.0) / 52.0, "carnegie_APPROX"),
    # Months.
    ("6-month-old human stage", 0.5, "exact_month"),
    # Coarse labels: parsed, but flagged.
    ("human adult stage", 54.5, "coarse_adult_VERY_IMPRECISE"),
    ("newborn human stage", 28.0 / 365.0 / 2, "coarse_newborn"),
    ("adolescent stage", 16.0, "coarse_adolescent"),
    # Unparseable must be NaN, never 0.0.
    ("", float("nan"), "unparsed"),
    ("unknown", float("nan"), "unparsed"),
    (None, float("nan"), "unparsed"),
]


def run_self_test(verbose: bool = True) -> tuple[int, int]:
    """Run the offline parser test suite. Returns (passed, total)."""
    passed = 0
    failures: list[str] = []
    for label, want_years, want_tag in _CASES:
        got = parse_development_stage(label)
        years_ok = (
            math.isnan(want_years) and math.isnan(got.years)
        ) or (
            not math.isnan(want_years)
            and not math.isnan(got.years)
            and abs(got.years - want_years) < 1e-6
        )
        tag_ok = got.precision == want_tag
        if years_ok and tag_ok:
            passed += 1
        else:
            failures.append(
                "  FAIL %-40r  want (%.6f, %s)  got (%.6f, %s)"
                % (label, want_years, want_tag, got.years, got.precision)
            )

    total = len(_CASES)
    if verbose:
        print("age parser self-test: %d/%d passed" % (passed, total))
        for line in failures:
            print(line)
        # Monotonicity across birth is the property the sign convention exists
        # to guarantee; assert it explicitly rather than trusting the cases.
        ordered = [
            "Carnegie stage 13",
            "8 post-conception weeks human stage",
            "20 post-conception week human stage",
            "newborn human stage",
            "6-month-old human stage",
            "2-year-old human stage",
            "adolescent stage",
            "third decade human stage",
            "45-49 year-old stage",
            "89-year-old human stage",
        ]
        vals = [parse_development_stage(x).years for x in ordered]
        mono = all(a < b for a, b in zip(vals, vals[1:]))
        print("monotonic across birth: %s" % ("PASS" if mono else "FAIL"))
        if not mono:
            for lab, v in zip(ordered, vals):
                print("    %-45s %+.4f" % (lab, v))
    return passed, total


if __name__ == "__main__":
    import sys

    ok, total = run_self_test()
    sys.exit(0 if ok == total else 1)

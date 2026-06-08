"""
Regression + robustness harness for IRSFunctions._verify_binary_lens_roots.

Goals
-----
1. In-scope mass ratios (q > 1e-6): outputs of _get_n_images / _image_calculator
   must be unchanged after switching to the tolerance-based fallback.
2. Out-of-scope mass ratios (q <= 1e-6): cases that previously raised
   "Wrong number of solutions" must now return a physically meaningful result
   (3 or 5 images, finite magnification that agrees with an independent estimate).

Usage
-----
    python tests/annulus_root_verification_test.py baseline   # run on OLD code, saves JSON
    python tests/annulus_root_verification_test.py current    # run on NEW code, saves JSON + diff

The two JSON files are compared by compare_baseline().
"""
import os
import sys
import json
import argparse

os.environ.setdefault("MPLCONFIGDIR", "/tmp/mplconfig_nblens")
os.makedirs(os.environ["MPLCONFIGDIR"], exist_ok=True)

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)

import numpy as np
from IRSMicroLensing.IRSFunctions import IRSFunctions as F
import VBMicrolensing

OUT_DIR = os.path.join(REPO, "tests", "_verification_artifacts")
os.makedirs(OUT_DIR, exist_ok=True)

# Independent reference: VBMicrolensing point-source binary magnification. This uses a
# completely different algorithm from the WM95 polynomial + root-verification path, so it
# validates the fixed out-of-scope points where MulensModel itself raises. Verified to agree
# with MulensModel WM95 to ~1e-10 with the identity mapping BinaryMag0(s, q, x_cm, y_cm).
_VBM = VBMicrolensing.VBMicrolensing()


def vbm_reference(x, y, q, s):
    return _VBM.BinaryMag0(s, q, x, y)


# (label, q, s, in_scope)
SCENARIOS = [
    ("inscope_q1e-3_s0.8",  1.0e-3, 0.80, True),
    ("inscope_q1e-3_s1.2",  1.0e-3, 1.20, True),
    ("inscope_q1e-2_s1.0",  1.0e-2, 1.00, True),
    ("inscope_q1e-4_s1.5",  1.0e-4, 1.50, True),
    ("inscope_q5e-3_s0.5",  5.0e-3, 0.50, True),
    ("inscope_q1e-3_s2.0",  1.0e-3, 2.00, True),
    # out-of-scope (extreme mass ratio)
    ("oos_mercury",         2.357392042097786e-07, 0.07683837479958969, False),
    ("oos_q1e-7_s0.5",      1.0e-7, 0.50, False),
    ("oos_q1e-8_s1.0",      1.0e-8, 1.00, False),
    ("oos_q1e-6_s0.3",      1.0e-6, 0.30, False),
    ("oos_q5e-7_s2.0",      5.0e-7, 2.00, False),
]

U_VALUES = [0.1, 0.2426, 0.5, 1.0, 1.1137, 1.5, 2.0]
N_THETA = 144


def get_n_images_point(x, y, q, s):
    """Return (count, magnification) or ('ERROR', msg)."""
    try:
        mag, imgs = F._get_n_images(x, y, q, s)
        return (len(imgs), float(mag))
    except Exception as e:  # noqa
        return ("ERROR", str(e).splitlines()[0])


def run():
    results = {}
    for (label, q, s, in_scope) in SCENARIOS:
        cm_x = q * s / (1.0 + q)  # center-of-mass x offset used by _image_calculator
        star_x = -s / (1.0 + q)   # primary (star) x-position in the CM frame
        results[label] = {"q": q, "s": s, "in_scope": in_scope, "by_u": {}}
        for u in U_VALUES:
            thetas = np.linspace(0, 2 * np.pi, N_THETA, endpoint=False)
            points = []          # per-theta: [count_or_"ERROR", mag_or_None]
            vbm_relerr = []
            for th in thetas:
                # source on ring of radius u about origin (CM frame), as _image_calculator does
                x = u * np.cos(th) - cm_x
                y = u * np.sin(th)
                res = get_n_images_point(x, y, q, s)
                if res[0] == "ERROR":
                    points.append(["ERROR", None])
                else:
                    cnt, mag = res
                    points.append([cnt, mag])
                    # independent cross-check against VBMicrolensing (different algorithm)
                    ref = vbm_reference(x, y, q, s)
                    if np.isfinite(ref) and ref > 0:
                        vbm_relerr.append(abs(mag - ref) / ref)

            # _image_calculator aggregate
            try:
                yp, ym, mm = F._image_calculator(u, q, s)
                imgcalc = [float(yp), float(ym), float(mm)]
            except Exception as e:  # noqa
                imgcalc = ["ERROR", str(e).splitlines()[0]]

            n_error = sum(1 for p in points if p[0] == "ERROR")
            counts = {}
            for p in points:
                if p[0] != "ERROR":
                    counts[p[0]] = counts.get(p[0], 0) + 1

            results[label]["by_u"][str(u)] = {
                "points": points,
                "counts": {str(k): v for k, v in sorted(counts.items())},
                "n_error": n_error,
                "n_points": N_THETA,
                "vbm_relerr_max": (float(np.max(vbm_relerr)) if vbm_relerr else None),
                "image_calculator": imgcalc,
            }
    return results


def summarize(results):
    print("\n================ SUMMARY ================")
    for label, d in results.items():
        tag = "IN " if d["in_scope"] else "OOS"
        total_err = sum(v["n_error"] for v in d["by_u"].values())
        all_counts = set()
        max_vbm = 0.0
        for v in d["by_u"].values():
            all_counts |= set(int(k) for k in v["counts"].keys())
            if v["vbm_relerr_max"] is not None:
                max_vbm = max(max_vbm, v["vbm_relerr_max"])
        print(f"[{tag}] {label:24s} q={d['q']:.2e} s={d['s']:.3f}  "
              f"counts={sorted(all_counts)} total_errors={total_err}  vbm_relerr<= {max_vbm:.2e}")


def compare_baseline(baseline, current):
    """
    Per-point comparison.

    Rule: any point that produced a VALID result in the baseline must produce a
    byte-identical (count, magnification) result now. Points that ERRORED in the
    baseline are allowed to become valid (that is the fix) but must never silently
    change into a *different* valid number.
    """
    print("\n================ BASELINE vs CURRENT (per-point) ================")
    regressions = 0
    fixed_points = 0
    newly_broken = 0
    for label in baseline:
        d_base = baseline[label]
        d_cur = current[label]
        in_scope = d_base["in_scope"]
        label_fixed = 0
        for u in d_base["by_u"]:
            bpts = d_base["by_u"][u]["points"]
            cpts = d_cur["by_u"][u]["points"]
            for idx, (bp, cp) in enumerate(zip(bpts, cpts)):
                b_ok = bp[0] != "ERROR"
                c_ok = cp[0] != "ERROR"
                if b_ok and c_ok:
                    # must be identical
                    same = (bp[0] == cp[0]) and (abs(bp[1] - cp[1]) <= 1e-12 * max(1.0, abs(bp[1])))
                    if not same:
                        regressions += 1
                        print(f"  REGRESSION [{ 'IN ' if in_scope else 'OOS'}] {label} u={u} i={idx}: "
                              f"{bp} -> {cp}")
                elif b_ok and not c_ok:
                    newly_broken += 1
                    print(f"  NEWLY-BROKEN [{ 'IN ' if in_scope else 'OOS'}] {label} u={u} i={idx}: {bp} -> ERROR")
                elif (not b_ok) and c_ok:
                    label_fixed += 1
                    fixed_points += 1
        if label_fixed:
            vbm = max((d_cur["by_u"][u]["vbm_relerr_max"] or 0.0) for u in d_cur["by_u"])
            tag = "IN " if in_scope else "OOS"
            print(f"  FIXED [{tag}] {label}: {label_fixed} previously-erroring point(s) now valid "
                  f"(max VBM relerr across this scenario = {vbm:.2e})")

    print("\n  ---- verdict ----")
    print(f"  in/out-of-scope previously-valid points changed (regressions): {regressions}")
    print(f"  previously-valid points newly broken:                          {newly_broken}")
    print(f"  previously-erroring points now fixed:                          {fixed_points}")
    if regressions == 0 and newly_broken == 0:
        print("  PASS: no previously-valid point changed; only error->valid improvements.")
    else:
        print("  FAIL: at least one previously-valid point changed.")
    return regressions + newly_broken


def validate(results, vbm_tol=1e-6):
    """
    Self-contained pass/fail (no stored baseline needed).

    Asserts:
      * no point raises "Wrong number of solutions" (zero errors anywhere),
      * every image count is physical (3 or 5),
      * every magnification agrees with the independent VBMicrolensing reference
        to within vbm_tol.
    """
    print("\n================ SELF-CONTAINED VALIDATION ================")
    failures = 0
    for label, d in results.items():
        n_err = sum(v["n_error"] for v in d["by_u"].values())
        counts = set()
        max_vbm = 0.0
        for v in d["by_u"].values():
            counts |= set(int(k) for k in v["counts"].keys())
            if v["vbm_relerr_max"] is not None:
                max_vbm = max(max_vbm, v["vbm_relerr_max"])
        bad_count = counts - {3, 5}
        ok = (n_err == 0) and (not bad_count) and (max_vbm <= vbm_tol)
        if not ok:
            failures += 1
        print(f"  [{'PASS' if ok else 'FAIL'}] {label:24s} errors={n_err} "
              f"counts={sorted(counts)} vbm_relerr<= {max_vbm:.2e}")
    print(f"\n  {'ALL PASS' if failures == 0 else str(failures) + ' SCENARIO(S) FAILED'}")
    return failures


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", nargs="?", default="validate",
                    choices=["validate", "baseline", "current"],
                    help="validate: self-contained VBM check (default). "
                         "baseline/current: capture JSON for a byte-identical per-point diff "
                         "(run 'baseline' on the pre-fix code, then 'current' on the new code).")
    args = ap.parse_args()

    results = run()
    summarize(results)

    if args.mode == "validate":
        sys.exit(1 if validate(results) else 0)

    path = os.path.join(OUT_DIR, f"{args.mode}.json")
    with open(path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nsaved -> {path}")

    if args.mode == "current":
        base_path = os.path.join(OUT_DIR, "baseline.json")
        if os.path.exists(base_path):
            with open(base_path) as f:
                baseline = json.load(f)
            issues = compare_baseline(baseline, results)
            sys.exit(1 if issues else 0)
        else:
            print("No baseline.json found; run with 'baseline' first.")


if __name__ == "__main__":
    main()

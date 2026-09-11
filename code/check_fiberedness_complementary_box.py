#!/usr/bin/env python3
"""HFK checks for a finite box in the complementary range p < r < p+q."""

from __future__ import annotations

import argparse
import csv
from concurrent.futures import ThreadPoolExecutor, as_completed
from math import gcd
import json
from pathlib import Path
import subprocess
import sys

from generalized_ttk_volume_experiment import Parameters, alexander_coefficients


FIELDS = (
    "p", "q", "r", "s", "hfk_status", "fibered", "top_alexander_rank",
    "seifert_genus", "alexander_breadth", "genus_matches_breadth", "tau",
    "l_space_knot", "total_rank", "simplified_crossings", "elapsed_seconds",
    "error",
)


def parameter_box(p_max: int, twists: tuple[int, ...]):
    for p in range(3, p_max + 1):
        for q in range(2, p):
            if gcd(p, q) != 1:
                continue
            for r in range(p + 1, p + q):
                for s in twists:
                    yield Parameters(p, q, r, s)


def run_one(parameters: Parameters, timeout: float, simplify: str):
    worker = Path(__file__).with_name("generalized_ttk_hfk_experiment.py")
    command = [
        sys.executable, str(worker), "--worker-json",
        str(parameters.p), str(parameters.q), str(parameters.r), str(parameters.s),
        simplify,
    ]
    try:
        completed = subprocess.run(
            command, capture_output=True, text=True, timeout=timeout, check=False
        )
        if completed.returncode:
            result = {
                "hfk_status": "error",
                "error": f"worker exit {completed.returncode}: {completed.stderr[-1000:]}",
            }
        else:
            lines = [line for line in completed.stdout.splitlines() if line.strip()]
            result = json.loads(lines[-1])
    except subprocess.TimeoutExpired:
        result = {
            "hfk_status": "timeout",
            "elapsed_seconds": str(timeout),
            "error": f"exceeded {timeout} seconds",
        }

    coefficients = alexander_coefficients(parameters)
    breadth = len(coefficients) - 1
    genus = result.get("seifert_genus", "")
    row = {
        "p": parameters.p,
        "q": parameters.q,
        "r": parameters.r,
        "s": parameters.s,
        "alexander_breadth": breadth,
        **result,
    }
    row["genus_matches_breadth"] = (
        str(2 * int(genus) == breadth).lower() if genus != "" else ""
    )
    return parameters, row


def write_csv(path: Path, rows: dict[Parameters, dict]):
    temporary = path.with_suffix(".csv.tmp")
    with temporary.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows[key] for key in sorted(rows))
    temporary.replace(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--p-max", type=int, default=6)
    parser.add_argument("--twists", type=int, nargs="+", default=(-2, -1, 1, 2))
    parser.add_argument("--timeout", type=float, default=180.0)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--simplify", choices=("none", "basic", "pickup"), default="pickup")
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).with_name("output")
        / "complementary_hfk_p3-6_s-2--1_1-2.csv",
    )
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)

    parameters = list(parameter_box(args.p_max, tuple(args.twists)))
    rows: dict[Parameters, dict] = {}
    if args.output.exists():
        with args.output.open(newline="", encoding="utf-8") as stream:
            for row in csv.DictReader(stream):
                key = Parameters(*(int(row[name]) for name in ("p", "q", "r", "s")))
                rows[key] = row
    pending = [key for key in parameters if rows.get(key, {}).get("hfk_status") != "computed"]
    print(f"cases={len(parameters)} pending={len(pending)}", flush=True)

    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {
            pool.submit(run_one, key, args.timeout, args.simplify): key for key in pending
        }
        for index, future in enumerate(as_completed(futures), 1):
            key, row = future.result()
            rows[key] = row
            write_csv(args.output, rows)
            print(
                f"{index}/{len(pending)} {key} {row['hfk_status']} "
                f"fibered={row.get('fibered', '')} time={row.get('elapsed_seconds', '')}",
                flush=True,
            )

    selected = [rows[key] for key in parameters]
    summary = {
        "total": len(selected),
        "computed": sum(row["hfk_status"] == "computed" for row in selected),
        "fibered": sum(row.get("fibered") == "true" for row in selected),
        "nonfibered": sum(row.get("fibered") == "false" for row in selected),
        "top_rank_not_one": sum(
            row["hfk_status"] == "computed" and row.get("top_alexander_rank") != "1"
            for row in selected
        ),
        "genus_mismatches": sum(row.get("genus_matches_breadth") == "false" for row in selected),
        "timeouts": sum(row["hfk_status"] == "timeout" for row in selected),
        "errors": sum(row["hfk_status"] == "error" for row in selected),
    }
    summary_path = args.output.with_suffix(".json")
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()

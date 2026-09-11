#!/usr/bin/env python3
"""Compare knot Floer homology inside equal-Alexander TTK groups.

This script consumes ``knots.csv`` from generalized_ttk_volume_experiment.py.
Each HFK computation runs in its own subprocess with a hard timeout, so a
difficult diagram cannot stall a long experiment.  The computed invariant is
the bigraded hat knot Floer homology over the coefficient field reported by
Spherogram (normally F_2).
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
from hashlib import sha256
import json
import os
from pathlib import Path
import subprocess
import sys
import time
from typing import Sequence

from generalized_ttk_volume_experiment import (
    Parameters,
    alexander_coefficients,
    generalized_braid_word,
)


HFK_FIELDS = (
    "p",
    "q",
    "r",
    "s",
    "alexander_sha256",
    "raw_crossings",
    "simplified_crossings",
    "source_alexander_breadth",
    "source_known_type",
    "source_volume_status",
    "source_volume",
    "source_isometry_signature",
    "hfk_status",
    "hfk_sha256",
    "coefficient_modulus",
    "seifert_genus",
    "tau",
    "fibered",
    "l_space_knot",
    "epsilon",
    "nu",
    "top_alexander_rank",
    "total_rank",
    "ranks_json",
    "elapsed_seconds",
    "error",
)


PAIR_FIELDS = (
    "alexander_sha256",
    "knot_1",
    "knot_2",
    "hfk_sha256_1",
    "hfk_sha256_2",
    "genus_1",
    "genus_2",
    "tau_1",
    "tau_2",
    "fibered_1",
    "fibered_2",
    "total_rank_1",
    "total_rank_2",
    "source_volume_1",
    "source_volume_2",
    "source_isometry_signature_1",
    "source_isometry_signature_2",
)


def _normalized_euler_characteristic(hfk: dict) -> tuple[int, ...]:
    coefficients: Counter[int] = Counter()
    for (alexander, maslov), rank in hfk["ranks"].items():
        sign = -1 if int(maslov) % 2 else 1
        coefficients[int(alexander)] += sign * int(rank)
    low, high = min(coefficients), max(coefficients)
    result = tuple(coefficients[index] for index in range(low, high + 1))
    if result[-1] < 0:
        result = tuple(-value for value in result)
    return result


def compute_hfk(parameters: Parameters, simplify: str) -> dict[str, str]:
    started = time.monotonic()
    result = {
        "hfk_status": "error",
        "simplified_crossings": "",
        "hfk_sha256": "",
        "coefficient_modulus": "",
        "seifert_genus": "",
        "tau": "",
        "fibered": "",
        "l_space_knot": "",
        "epsilon": "",
        "nu": "",
        "top_alexander_rank": "",
        "total_rank": "",
        "ranks_json": "",
        "elapsed_seconds": "",
        "error": "",
    }
    try:
        from spherogram import Link

        word, _ = generalized_braid_word(parameters)
        link = Link(braid_closure=word)
        if simplify in ("basic", "pickup"):
            link.simplify("basic")
        if simplify == "pickup":
            link.simplify("pickup")
        result["simplified_crossings"] = str(len(link.crossings))
        hfk = link.knot_floer_homology()
        expected = alexander_coefficients(parameters)
        computed = _normalized_euler_characteristic(hfk)
        if computed != expected:
            raise ArithmeticError(
                f"Euler characteristic does not match Alexander formula: {computed} != {expected}"
            )

        ranks = sorted(
            [int(alexander), int(maslov), int(rank)]
            for (alexander, maslov), rank in hfk["ranks"].items()
        )
        payload = {
            "modulus": int(hfk.get("modulus", 2)),
            "ranks": ranks,
        }
        encoded = json.dumps(payload, separators=(",", ":"), sort_keys=True)
        genus = int(hfk["seifert_genus"])
        result.update(
            {
                "hfk_status": "computed",
                "hfk_sha256": sha256(encoded.encode("ascii")).hexdigest(),
                "coefficient_modulus": str(payload["modulus"]),
                "seifert_genus": str(genus),
                "tau": str(hfk.get("tau", "")),
                "fibered": str(bool(hfk.get("fibered", False))).lower(),
                "l_space_knot": str(bool(hfk.get("L_space_knot", False))).lower(),
                "epsilon": str(hfk.get("epsilon", "")),
                "nu": str(hfk.get("nu", "")),
                "top_alexander_rank": str(
                    sum(rank for alexander, _, rank in ranks if alexander == genus)
                ),
                "total_rank": str(int(hfk["total_rank"])),
                "ranks_json": encoded,
            }
        )
    except Exception as error:
        result["error"] = f"{type(error).__name__}: {error}"
    result["elapsed_seconds"] = format(time.monotonic() - started, ".6f")
    return result


def _worker_main(argv: Sequence[str]) -> int:
    if len(argv) != 6:
        raise SystemExit("worker usage: --worker-json P Q R S SIMPLIFY")
    parameters = Parameters(*(int(value) for value in argv[1:5]))
    print(json.dumps(compute_hfk(parameters, argv[5]), separators=(",", ":")))
    return 0


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def _write_csv(path: Path, rows: list[dict[str, str]], fields: Sequence[str]) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def _write_json(path: Path, data: object) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(path)


def _key(row: dict[str, str]) -> Parameters:
    return Parameters(*(int(row[name]) for name in ("p", "q", "r", "s")))


def _run_subprocess(task: tuple[dict[str, str], argparse.Namespace]) -> tuple[Parameters, dict[str, str]]:
    source, args = task
    parameters = _key(source)
    command = [
        sys.executable,
        str(Path(__file__).resolve()),
        "--worker-json",
        str(parameters.p),
        str(parameters.q),
        str(parameters.r),
        str(parameters.s),
        args.diagram_simplify,
    ]
    environment = os.environ.copy()
    environment.setdefault("PYTHONPYCACHEPREFIX", "/tmp/ttk-hfk-pycache")
    try:
        completed = subprocess.run(
            command,
            capture_output=True,
            text=True,
            timeout=args.timeout,
            env=environment,
            check=False,
        )
        if completed.returncode != 0:
            return parameters, {
                "hfk_status": "error",
                "error": f"worker exit {completed.returncode}: {completed.stderr[-1000:]}",
            }
        lines = [line for line in completed.stdout.splitlines() if line.strip()]
        return parameters, json.loads(lines[-1])
    except subprocess.TimeoutExpired:
        return parameters, {
            "hfk_status": "timeout",
            "elapsed_seconds": str(args.timeout),
            "error": f"exceeded {args.timeout} seconds",
        }
    except Exception as error:
        return parameters, {
            "hfk_status": "error",
            "error": f"{type(error).__name__}: {error}",
        }


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path, help="knots.csv or its containing directory")
    parser.add_argument(
        "--selection",
        choices=("unresolved", "all-collisions", "all-knots"),
        default="unresolved",
        help="unresolved collision groups, all collision groups, or every source knot",
    )
    parser.add_argument("--max-crossings", type=int, default=60)
    parser.add_argument("--timeout", type=float, default=60.0)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--checkpoint-every", type=int, default=5)
    parser.add_argument(
        "--diagram-simplify", choices=("none", "basic", "pickup"), default="pickup"
    )
    parser.add_argument("--retry-failures", action="store_true")
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args(argv)
    if args.max_crossings < 1 or args.timeout <= 0 or args.workers < 1:
        parser.error("crossing bound, timeout, and workers must be positive")
    source = args.source.resolve()
    args.source_csv = source / "knots.csv" if source.is_dir() else source
    if not args.source_csv.exists():
        parser.error(f"source CSV not found: {args.source_csv}")
    if args.output_dir is None:
        args.output_dir = args.source_csv.parent / (
            f"hfk_{args.selection}_crossings{args.max_crossings}"
        )
    return args


def run(args: argparse.Namespace) -> dict:
    output_dir: Path = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    source_rows = _read_csv(args.source_csv)
    groups: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in source_rows:
        groups[row["alexander_sha256"]].append(row)

    selected_groups = []
    for members in groups.values():
        if args.selection == "all-knots":
            selected_groups.append(members)
        elif len(members) >= 2:
            if args.selection == "unresolved" and not any(
                row["volume_status"] == "unresolved" for row in members
            ):
                continue
            selected_groups.append(members)

    candidates = [
        row
        for members in selected_groups
        for row in members
        if int(row["braid_crossings_before_simplification"]) <= args.max_crossings
    ]
    if args.selection == "all-knots":
        eligible_digests = {row["alexander_sha256"] for row in candidates}
    else:
        eligible_digests = {
            digest
            for digest, members in groups.items()
            if sum(member in candidates for member in members) >= 2
        }
        candidates = [
            row for row in candidates if row["alexander_sha256"] in eligible_digests
        ]
    candidates.sort(key=lambda row: _key(row))

    result_path = output_dir / "hfk_knots.csv"
    previous = {_key(row): row for row in _read_csv(result_path)} if result_path.exists() else {}
    records: dict[Parameters, dict[str, str]] = {}
    pending: list[dict[str, str]] = []
    for source in candidates:
        parameters = _key(source)
        base = {
            "p": source["p"],
            "q": source["q"],
            "r": source["r"],
            "s": source["s"],
            "alexander_sha256": source["alexander_sha256"],
            "raw_crossings": source["braid_crossings_before_simplification"],
            "simplified_crossings": "",
            "source_alexander_breadth": source["alexander_breadth"],
            "source_known_type": source["known_type"],
            "source_volume_status": source["volume_status"],
            "source_volume": source["volume"],
            "source_isometry_signature": source["isometry_signature"],
            "hfk_status": "pending",
            "hfk_sha256": "",
            "coefficient_modulus": "",
            "seifert_genus": "",
            "tau": "",
            "fibered": "",
            "l_space_knot": "",
            "epsilon": "",
            "nu": "",
            "top_alexander_rank": "",
            "total_rank": "",
            "ranks_json": "",
            "elapsed_seconds": "",
            "error": "",
        }
        if parameters in previous:
            base.update(previous[parameters])
        records[parameters] = base
        failed = base["hfk_status"] in ("timeout", "error")
        if base["hfk_status"] in ("pending", "") or (args.retry_failures and failed):
            base["hfk_status"] = "pending"
            pending.append(source)

    ordered = lambda: [records[_key(row)] for row in candidates]
    _write_csv(result_path, ordered(), HFK_FIELDS)
    print(
        f"Selected groups: {len(eligible_digests)}; HFK candidates: {len(candidates)}; "
        f"pending: {len(pending)}",
        flush=True,
    )

    completed_count = 0
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(_run_subprocess, (row, args)): row for row in pending}
        for future in as_completed(futures):
            parameters, result = future.result()
            records[parameters].update(result)
            completed_count += 1
            if completed_count % args.checkpoint_every == 0 or completed_count == len(pending):
                _write_csv(result_path, ordered(), HFK_FIELDS)
                print(f"  HFK phase: {completed_count}/{len(pending)}", flush=True)

    computed_groups: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in ordered():
        if row["hfk_status"] == "computed":
            computed_groups[row["alexander_sha256"]].append(row)
    different_pairs: list[dict[str, str]] = []
    compared_groups = 0
    for digest, members in computed_groups.items():
        if len(members) < 2:
            continue
        compared_groups += 1
        for index, left in enumerate(members):
            for right in members[index + 1 :]:
                if left["hfk_sha256"] == right["hfk_sha256"]:
                    continue
                different_pairs.append(
                    {
                        "alexander_sha256": digest,
                        "knot_1": f"T({left['p']},{left['q']};{left['r']},{left['s']})",
                        "knot_2": f"T({right['p']},{right['q']};{right['r']},{right['s']})",
                        "hfk_sha256_1": left["hfk_sha256"],
                        "hfk_sha256_2": right["hfk_sha256"],
                        "genus_1": left["seifert_genus"],
                        "genus_2": right["seifert_genus"],
                        "tau_1": left["tau"],
                        "tau_2": right["tau"],
                        "fibered_1": left["fibered"],
                        "fibered_2": right["fibered"],
                        "total_rank_1": left["total_rank"],
                        "total_rank_2": right["total_rank"],
                        "source_volume_1": left["source_volume"],
                        "source_volume_2": right["source_volume"],
                        "source_isometry_signature_1": left["source_isometry_signature"],
                        "source_isometry_signature_2": right["source_isometry_signature"],
                    }
                )
    _write_csv(output_dir / "same_alexander_different_hfk.csv", different_pairs, PAIR_FIELDS)

    completed_rows = [row for row in ordered() if row["hfk_status"] == "computed"]
    non_l_space = [row for row in completed_rows if row["l_space_knot"] == "false"]
    nonfibered = [row for row in completed_rows if row["fibered"] == "false"]
    genus_gap = [
        row
        for row in completed_rows
        if 2 * int(row["seifert_genus"]) != int(row["source_alexander_breadth"])
    ]
    _write_csv(output_dir / "non_l_space_knots.csv", non_l_space, HFK_FIELDS)
    _write_csv(output_dir / "nonfibered_knots.csv", nonfibered, HFK_FIELDS)
    _write_csv(output_dir / "alexander_genus_gap.csv", genus_gap, HFK_FIELDS)

    status_counts = Counter(row["hfk_status"] for row in ordered())
    computed_new_region = [row for row in completed_rows if int(row["r"]) > int(row["p"])]
    eligible_collision_groups = sum(
        sum(row["alexander_sha256"] == digest for row in candidates) >= 2
        for digest in eligible_digests
    )
    summary = {
        "source_csv": str(args.source_csv),
        "selection": args.selection,
        "max_raw_crossings": args.max_crossings,
        "timeout_seconds": args.timeout,
        "eligible_alexander_groups": len(eligible_digests),
        "eligible_alexander_collision_groups": eligible_collision_groups,
        "hfk_candidates": len(candidates),
        "hfk_compared_groups": compared_groups,
        "same_alexander_different_hfk_pairs": len(different_pairs),
        "computed_new_region_knots": len(computed_new_region),
        "non_l_space_knots": len(non_l_space),
        "non_l_space_new_region_knots": sum(
            int(row["r"]) > int(row["p"]) for row in non_l_space
        ),
        "nonfibered_knots": len(nonfibered),
        "alexander_breadth_less_than_twice_hfk_genus": len(genus_gap),
        **{f"status_{key}": value for key, value in sorted(status_counts.items())},
    }
    _write_json(output_dir / "hfk_summary.json", summary)
    print(json.dumps(summary, indent=2, sort_keys=True), flush=True)
    print(f"Results: {output_dir}", flush=True)
    return summary


def main(argv: Sequence[str] | None = None) -> int:
    actual = list(sys.argv[1:] if argv is None else argv)
    if actual and actual[0] == "--worker-json":
        return _worker_main(actual)
    return 0 if run(parse_args(actual)) else 1


if __name__ == "__main__":
    raise SystemExit(main())

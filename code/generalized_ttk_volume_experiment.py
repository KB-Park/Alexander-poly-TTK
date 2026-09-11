#!/usr/bin/env python3
"""Compare Alexander polynomials and volumes of generalized twisted torus knots.

The parameter convention is T(p,q;r,s), where

    p > q > 0, gcd(p,q) = 1, 0 < r < p+q, and s > 0.

For r <= p, the knot is represented by the usual p-strand braid

    (sigma_1 ... sigma_{p-1})^q (sigma_1 ... sigma_{r-1})^(rs).

For p < r < p+q, we use the positive (p+q)-strand Lorenz shuffle braid
for T(p,q), followed by the full twist on its first r strands.  This avoids
the expensive breadth-first Markov-destabilization search used by the
exploratory scripts in this repository.

The program is designed for experiments, not for proving hyperbolicity.
Without ``--certify``, "hyperbolic_numerical" means that SnapPy found a
positively oriented geometric triangulation.  Run the same script under
Sage with ``--certify`` to use interval arithmetic::

    DOT_SAGE=/tmp/sage-ttk sage generalized_ttk_volume_experiment.py ... --certify

Outputs are checkpointed in CSV files, so an interrupted volume run resumes
when invoked again with the same output directory and configuration.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
import csv
from dataclasses import asdict, dataclass
from hashlib import sha256
import json
from math import gcd
from pathlib import Path
import sys
import time
from typing import Iterable, Iterator, Sequence


CSV_FIELDS = (
    "p",
    "q",
    "r",
    "s",
    "alexander_sha256",
    "alexander_breadth",
    "alexander_nonzero_terms",
    "braid_strands",
    "braid_crossings_before_simplification",
    "braid_crossings_after_simplification",
    "known_type",
    "volume_status",
    "solution_type",
    "volume",
    "volume_accuracy_digits",
    "volume_verified",
    "isometry_signature",
    "attempts",
    "error",
)


@dataclass(frozen=True, order=True)
class Parameters:
    p: int
    q: int
    r: int
    s: int


def validate_parameters(parameters: Parameters) -> None:
    p, q, r, s = asdict(parameters).values()
    if not (p > q > 0):
        raise ValueError("require p > q > 0")
    if gcd(p, q) != 1:
        raise ValueError("require gcd(p,q) = 1")
    if not (0 < r < p + q):
        raise ValueError("require 0 < r < p+q")


def _convolve(left: Counter[int], right: Counter[int]) -> Counter[int]:
    result: Counter[int] = Counter()
    for left_exp, left_coeff in left.items():
        for right_exp, right_coeff in right.items():
            result[left_exp + right_exp] += left_coeff * right_coeff
    # Unary ``+`` on a Counter discards negative coefficients, which are
    # essential here.
    return Counter({exponent: value for exponent, value in result.items() if value})


def _divide_by_geometric_sum(
    numerator: Counter[int], length: int
) -> tuple[int, ...]:
    """Exactly divide a Laurent polynomial by 1+t+...+t^(length-1).

    The ascending recurrence is linear in the exponent span and avoids
    SymPy's much heavier symbolic polynomial division.
    """
    if not numerator:
        raise ArithmeticError("Alexander numerator unexpectedly vanished")
    low = min(numerator)
    high = max(numerator)
    quotient: list[int] = []
    previous_window = 0
    for exponent in range(low, high + 1):
        coefficient = numerator.get(exponent, 0) - previous_window
        quotient.append(coefficient)
        previous_window += coefficient
        if len(quotient) >= length:
            previous_window -= quotient[-length]

    if previous_window != 0:
        raise ArithmeticError("nonzero remainder in Alexander division")
    while len(quotient) > 1 and quotient[-1] == 0:
        quotient.pop()
    while len(quotient) > 1 and quotient[0] == 0:
        quotient.pop(0)
    if quotient[-1] < 0:
        quotient = [-coefficient for coefficient in quotient]
    return tuple(quotient)


def torus_alexander_coefficients(p: int, q: int) -> tuple[int, ...]:
    """Normalized dense coefficients of the Alexander polynomial of T(p,q)."""
    if p <= 0 or q <= 0 or gcd(p, q) != 1:
        raise ValueError("torus parameters must be positive and coprime")
    # (1-t)(1-t^(pq)) / ((1-t^p)(1-t^q)).
    numerator: Counter[int] = Counter()
    for exponent, coefficient in ((0, 1), (1, -1), (p * q, -1), (p * q + 1, 1)):
        numerator[exponent] += coefficient
        if numerator[exponent] == 0:
            del numerator[exponent]
    denominator = _convolve(
        Counter({0: 1, p: -1}), Counter({0: 1, q: -1})
    )
    degree = (p - 1) * (q - 1)
    remainder = [numerator.get(i, 0) for i in range(p * q + 2)]
    quotient: list[int] = []
    denominator_degree = max(denominator)
    for exponent in range(degree + 1):
        coefficient = remainder[exponent]
        quotient.append(coefficient)
        for offset, value in denominator.items():
            remainder[exponent + offset] -= coefficient * value
    if any(remainder):
        raise ArithmeticError("torus Alexander polynomial division was not exact")
    if len(quotient) != degree + 1 or denominator_degree != p + q:
        raise AssertionError("unexpected torus Alexander polynomial degree")
    if quotient[-1] < 0:
        quotient = [-coefficient for coefficient in quotient]
    return tuple(quotient)


def alexander_coefficients(parameters: Parameters) -> tuple[int, ...]:
    """Return normalized dense coefficients from the unified formula.

    The lowest exponent is zero and the leading coefficient is positive.
    The r=1 case is handled separately because it is just T(p,q), including
    the q=1 unknot boundary that is outside the formula's stated hypotheses.
    """
    validate_parameters(parameters)
    p, q, r, s = (
        parameters.p,
        parameters.q,
        parameters.r,
        parameters.s,
    )
    if r == 1:
        return torus_alexander_coefficients(p, q)

    j_star = (-r * pow(p, -1, q)) % q
    s_values = [(-j * p) % q for j in range(q)]
    lengths_l = [(r - 1 - value) // q + 1 for value in s_values]

    # Use the integer formula, not merely its congruence modulo p.  The
    # distinction matters at r=q, where i_star=p rather than 0.
    numerator_i = p * j_star + r
    if numerator_i % q:
        raise ArithmeticError("the unified index i_star is not integral")
    i_star = p + 1 - numerator_i // q
    t_values = [(-i * q) % p for i in range(p)]
    lengths_m = [(r - 1 - value) // p + 1 for value in t_values]

    cumulative_l: list[int] = []
    total = 0
    for value in lengths_l:
        total += value
        cumulative_l.append(total)

    cumulative_m: list[int] = []
    total = 0
    for value in lengths_m:
        total += value
        cumulative_m.append(total)

    lambda_star = (lengths_l[j_star] - 1) * s + 1
    e_a = (
        r * s * (cumulative_l[j_star - 1] if j_star else 0)
        + p * j_star
        + r * lambda_star
    )
    e_c = r * r * s + p * q

    x_a = Counter(
        r * s * cumulative_l[m] + p * m for m in range(j_star)
    )
    x_c = Counter(r * s * cumulative_l[m] + p * m for m in range(q))

    suffix_m = [0] * (p + 1)
    for index in range(p - 1, -1, -1):
        suffix_m[index] = suffix_m[index + 1] + lengths_m[index]
    z_a = Counter(
        e_a - q * h - r * s * suffix_m[p - h]
        for h in range(p - i_star + 1)
    )
    z_c = Counter(
        e_c - r * s * cumulative_m[m] - q * m for m in range(p)
    )

    numerator = _convolve(z_a, x_c)
    for exponent, coefficient in _convolve(x_a, z_c).items():
        numerator[exponent] -= coefficient
        if numerator[exponent] == 0:
            del numerator[exponent]

    coefficients = _divide_by_geometric_sum(numerator, r)
    if coefficients != coefficients[::-1]:
        raise ArithmeticError(
            f"Alexander polynomial is not symmetric for {parameters}"
        )
    return coefficients


def alexander_digest(coefficients: Sequence[int]) -> str:
    encoded = ",".join(map(str, coefficients)).encode("ascii")
    return sha256(encoded).hexdigest()


def polynomial_string(coefficients: Sequence[int]) -> str:
    """Compact ascending-exponent polynomial string for collision reports."""
    terms: list[str] = []
    for exponent, coefficient in enumerate(coefficients):
        if coefficient == 0:
            continue
        magnitude = abs(coefficient)
        if exponent == 0:
            body = str(magnitude)
        elif exponent == 1:
            body = "t" if magnitude == 1 else f"{magnitude}*t"
        else:
            body = f"t^{exponent}" if magnitude == 1 else f"{magnitude}*t^{exponent}"
        if not terms:
            terms.append(body if coefficient > 0 else f"-{body}")
        else:
            terms.append((" + " if coefficient > 0 else " - ") + body)
    return "".join(terms) or "0"


def full_twist_word(r: int, s: int) -> list[int]:
    if r <= 1:
        return []
    if s > 0:
        return list(range(1, r)) * (r * s)
    if s < 0:
        return list(range(-1, -r, -1)) * (-r * s)
    return []


def lorenz_shuffle_word(p: int, q: int) -> list[int]:
    """The standard positive (p+q)-strand Lorenz braid for T(p,q)."""
    word: list[int] = []
    for j in range(q):
        word.extend(range(p + j, j, -1))
    return word


def generalized_braid_word(parameters: Parameters) -> tuple[list[int], int]:
    """Return a deterministic braid word and its explicitly intended strand count."""
    validate_parameters(parameters)
    p, q, r, s = parameters.p, parameters.q, parameters.r, parameters.s
    twist = full_twist_word(r, s)
    if r <= p:
        return list(range(1, p)) * q + twist, p
    return lorenz_shuffle_word(p, q) + twist, p + q


def known_torus_type(parameters: Parameters) -> str:
    """Return a torus-knot description for elementary non-hyperbolic cases."""
    p, q, r, s = parameters.p, parameters.q, parameters.r, parameters.s
    if r == 1:
        return f"T({p},{q})"
    if q == 1:
        return f"T({r},{r * s + 1})"
    if r == q:
        return f"T({q},{p + q * s})"
    if r == p:
        return f"T({p},{q + p * s})"
    return ""


def parameter_range(args: argparse.Namespace) -> Iterator[Parameters]:
    for p in range(args.p_min, args.p_max + 1):
        q_high = min(p - 1, args.q_max if args.q_max is not None else p - 1)
        for q in range(args.q_min, q_high + 1):
            if gcd(p, q) != 1:
                continue
            if args.region == "old":
                low, high = 2, p - 1
            elif args.region == "new":
                low, high = p + 1, p + q - 1
            else:
                low, high = 1, p + q - 1
            low = max(low, args.r_min if args.r_min is not None else low)
            high = min(high, args.r_max if args.r_max is not None else high)
            for r in range(low, high + 1):
                for s in range(args.s_min, args.s_max + 1):
                    yield Parameters(p, q, r, s)


def _empty_record(parameters: Parameters, coefficients: Sequence[int]) -> dict[str, str]:
    word, strands = generalized_braid_word(parameters)
    torus_type = known_torus_type(parameters)
    status = "known_nonhyperbolic_torus" if torus_type else "not_requested"
    return {
        "p": str(parameters.p),
        "q": str(parameters.q),
        "r": str(parameters.r),
        "s": str(parameters.s),
        "alexander_sha256": alexander_digest(coefficients),
        "alexander_breadth": str(len(coefficients) - 1),
        "alexander_nonzero_terms": str(sum(value != 0 for value in coefficients)),
        "braid_strands": str(strands),
        "braid_crossings_before_simplification": str(len(word)),
        "braid_crossings_after_simplification": "",
        "known_type": torus_type,
        "volume_status": status,
        "solution_type": "",
        "volume": "",
        "volume_accuracy_digits": "",
        "volume_verified": "false",
        "isometry_signature": "",
        "attempts": "0",
        "error": "",
    }


def _has_sage() -> bool:
    try:
        from sage.all import RealIntervalField  # noqa: F401
    except Exception:
        return False
    return True


def _compute_geometry(task: tuple[Parameters, dict]) -> tuple[Parameters, dict[str, str]]:
    parameters, options = task
    result = {
        "braid_crossings_after_simplification": "",
        "volume_status": "error",
        "solution_type": "",
        "volume": "",
        "volume_accuracy_digits": "",
        "volume_verified": "false",
        "isometry_signature": "",
        "attempts": "0",
        "error": "",
    }
    try:
        import snappy

        word, _ = generalized_braid_word(parameters)
        link = snappy.Link(braid_closure=word)
        if len(link.link_components) != 1:
            raise ArithmeticError(
                f"braid closure has {len(link.link_components)} components"
            )
        simplify = options["diagram_simplify"]
        if simplify in ("basic", "pickup"):
            link.simplify("basic")
        if simplify == "pickup":
            link.simplify("pickup")
        result["braid_crossings_after_simplification"] = str(len(link.crossings))

        manifold = link.exterior()
        numerical_candidate: tuple[str, float, int] | None = None
        verified = False
        attempts = options["retries"] + 1
        for attempt in range(attempts):
            result["attempts"] = str(attempt + 1)
            solution_type = manifold.solution_type()
            result["solution_type"] = solution_type
            volume, accuracy = manifold.volume(accuracy=True)
            volume_float = float(volume)
            if (
                solution_type == "all tetrahedra positively oriented"
                and volume_float > options["volume_epsilon"]
            ):
                numerical_candidate = (solution_type, volume_float, int(accuracy))
                if options["certify"]:
                    try:
                        verified, _ = manifold.verify_hyperbolicity(
                            bits_prec=options["bits_prec"]
                        )
                    except Exception:
                        verified = False
                    if verified:
                        interval = manifold.volume(
                            verified=True, bits_prec=options["bits_prec"]
                        )
                        result["volume"] = format(float(interval), ".17g")
                        result["volume_accuracy_digits"] = str(options["bits_prec"])
                        result["volume_verified"] = "true"
                        result["volume_status"] = "hyperbolic_verified"
                        break
                else:
                    high_precision = manifold.high_precision()
                    hp_volume, hp_accuracy = high_precision.volume(accuracy=True)
                    result["volume"] = format(float(hp_volume), ".17g")
                    result["volume_accuracy_digits"] = str(int(hp_accuracy))
                    result["volume_status"] = "hyperbolic_numerical"
                    break
            if attempt + 1 < attempts:
                manifold.randomize()
                manifold.simplify()

        if options["certify"] and not verified:
            result["volume_status"] = "unresolved"
            if numerical_candidate is not None:
                _, volume_float, accuracy = numerical_candidate
                result["volume"] = format(volume_float, ".17g")
                result["volume_accuracy_digits"] = str(accuracy)
        elif result["volume_status"] == "error":
            result["volume_status"] = "unresolved"

        if (
            options["isometry_signatures"]
            and result["volume_status"].startswith("hyperbolic_")
        ):
            try:
                result["isometry_signature"] = manifold.isometry_signature(
                    of_link=True,
                    verified=options["certify"],
                )
            except Exception as error:
                result["error"] = f"isometry signature: {type(error).__name__}: {error}"
    except Exception as error:
        result["error"] = f"{type(error).__name__}: {error}"
    return parameters, result


def _read_csv(path: Path) -> dict[Parameters, dict[str, str]]:
    if not path.exists():
        return {}
    rows: dict[Parameters, dict[str, str]] = {}
    with path.open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            key = Parameters(*(int(row[name]) for name in ("p", "q", "r", "s")))
            rows[key] = row
    return rows


def _write_csv(path: Path, rows: Iterable[dict[str, str]], fields=CSV_FIELDS) -> None:
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


def _hyperbolic(row: dict[str, str]) -> bool:
    return row["volume_status"] in ("hyperbolic_numerical", "hyperbolic_verified")


def _pair_row(left: dict[str, str], right: dict[str, str]) -> dict[str, str]:
    return {
        "knot_1": f"T({left['p']},{left['q']};{left['r']},{left['s']})",
        "knot_2": f"T({right['p']},{right['q']};{right['r']},{right['s']})",
        "alexander_sha256_1": left["alexander_sha256"],
        "alexander_sha256_2": right["alexander_sha256"],
        "volume_1": left["volume"],
        "volume_2": right["volume"],
        "volume_difference": format(abs(float(left["volume"]) - float(right["volume"])), ".17g"),
        "isometry_signature_1": left["isometry_signature"],
        "isometry_signature_2": right["isometry_signature"],
    }


PAIR_FIELDS = (
    "knot_1",
    "knot_2",
    "alexander_sha256_1",
    "alexander_sha256_2",
    "volume_1",
    "volume_2",
    "volume_difference",
    "isometry_signature_1",
    "isometry_signature_2",
)


def write_comparison_reports(
    output_dir: Path,
    records: Sequence[dict[str, str]],
    coefficients_by_digest: dict[str, tuple[int, ...]],
    tolerance: float,
    volume_scope: str,
) -> dict[str, int | float | str]:
    groups: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in records:
        groups[row["alexander_sha256"]].append(row)

    collision_rows: list[dict[str, str]] = []
    same_alex_different_volume: list[dict[str, str]] = []
    same_alex_different_signature: list[dict[str, str]] = []
    mixed_geometry_groups = 0
    for digest, members in sorted(groups.items()):
        if len(members) < 2:
            continue
        hyperbolic = [row for row in members if _hyperbolic(row)]
        known_torus = [
            row for row in members if row["volume_status"] == "known_nonhyperbolic_torus"
        ]
        volumes = [float(row["volume"]) for row in hyperbolic]
        spread = max(volumes) - min(volumes) if volumes else 0.0
        signatures = {
            row["isometry_signature"]
            for row in hyperbolic
            if row["isometry_signature"]
        }
        if known_torus and hyperbolic:
            verdict = "known_torus_and_hyperbolic"
            mixed_geometry_groups += 1
        elif len(signatures) > 1:
            verdict = "different_link_signatures"
        elif len(hyperbolic) < 2:
            verdict = "insufficient_hyperbolic_members"
        elif spread > tolerance:
            verdict = "different_volumes"
        else:
            verdict = "same_volume_within_tolerance"
        collision_rows.append(
            {
                "alexander_sha256": digest,
                "alexander_polynomial": polynomial_string(coefficients_by_digest[digest]),
                "group_size": str(len(members)),
                "hyperbolic_members": str(len(hyperbolic)),
                "known_torus_members": str(len(known_torus)),
                "distinct_link_signatures": str(len(signatures)),
                "volume_spread": format(spread, ".17g"),
                "verdict": verdict,
                "members": json.dumps(
                    [
                        {
                            "parameters": [row[name] for name in ("p", "q", "r", "s")],
                            "status": row["volume_status"],
                            "volume": row["volume"],
                            "isometry_signature": row["isometry_signature"],
                        }
                        for row in members
                    ],
                    separators=(",", ":"),
                ),
            }
        )
        for index, left in enumerate(hyperbolic):
            for right in hyperbolic[index + 1 :]:
                if abs(float(left["volume"]) - float(right["volume"])) > tolerance:
                    same_alex_different_volume.append(_pair_row(left, right))
                if (
                    left["isometry_signature"]
                    and right["isometry_signature"]
                    and left["isometry_signature"] != right["isometry_signature"]
                ):
                    same_alex_different_signature.append(_pair_row(left, right))

    collision_fields = (
        "alexander_sha256",
        "alexander_polynomial",
        "group_size",
        "hyperbolic_members",
        "known_torus_members",
        "distinct_link_signatures",
        "volume_spread",
        "verdict",
        "members",
    )
    _write_csv(output_dir / "alexander_collision_groups.csv", collision_rows, collision_fields)
    _write_csv(
        output_dir / "same_alexander_different_volume.csv",
        same_alex_different_volume,
        PAIR_FIELDS,
    )
    _write_csv(
        output_dir / "same_alexander_different_signature.csv",
        same_alex_different_signature,
        PAIR_FIELDS,
    )

    same_volume_different_alexander: list[dict[str, str]] = []
    if volume_scope == "all":
        hyperbolic_records = sorted(
            (row for row in records if _hyperbolic(row)), key=lambda row: float(row["volume"])
        )
        for index, left in enumerate(hyperbolic_records):
            left_volume = float(left["volume"])
            for right in hyperbolic_records[index + 1 :]:
                difference = float(right["volume"]) - left_volume
                if difference > tolerance:
                    break
                if left["alexander_sha256"] != right["alexander_sha256"]:
                    same_volume_different_alexander.append(_pair_row(left, right))
    _write_csv(
        output_dir / "same_volume_different_alexander.csv",
        same_volume_different_alexander,
        PAIR_FIELDS,
    )

    counts = Counter(row["volume_status"] for row in records)
    verdict_counts = Counter(row["verdict"] for row in collision_rows)
    summary: dict[str, int | float | str] = {
        "total_knots": len(records),
        "unique_alexander_polynomials": len(groups),
        "alexander_collision_groups": len(collision_rows),
        "knots_in_alexander_collision_groups": sum(
            len(members) for members in groups.values() if len(members) > 1
        ),
        "same_alexander_different_volume_pairs": len(same_alex_different_volume),
        "same_alexander_different_signature_pairs": len(same_alex_different_signature),
        "alexander_groups_mixing_known_torus_and_hyperbolic": mixed_geometry_groups,
        "same_volume_different_alexander_pairs": len(same_volume_different_alexander),
        "volume_tolerance": tolerance,
        "volume_scope": volume_scope,
    }
    summary.update({f"status_{key}": value for key, value in sorted(counts.items())})
    summary.update(
        {f"collision_verdict_{key}": value for key, value in sorted(verdict_counts.items())}
    )
    _write_json(output_dir / "summary.json", summary)
    return summary


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--p-min", type=int, default=3)
    parser.add_argument("--p-max", type=int, default=10)
    parser.add_argument("--q-min", type=int, default=1)
    parser.add_argument("--q-max", type=int)
    parser.add_argument("--r-min", type=int)
    parser.add_argument("--r-max", type=int)
    parser.add_argument("--s-min", type=int, default=1)
    parser.add_argument("--s-max", type=int, default=3)
    parser.add_argument(
        "--region",
        choices=("all", "old", "new"),
        default="all",
        help="all: 0<r<p+q; old: 1<r<p; new: p<r<p+q",
    )
    parser.add_argument(
        "--volume-scope",
        choices=("collisions", "all"),
        default="collisions",
        help="compute only Alexander-collision groups, or every candidate",
    )
    parser.add_argument("--volume-tolerance", type=float, default=1e-8)
    parser.add_argument("--volume-epsilon", type=float, default=1e-7)
    parser.add_argument("--retries", type=int, default=5)
    parser.add_argument("--bits-prec", type=int, default=100)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--checkpoint-every", type=int, default=10)
    parser.add_argument(
        "--diagram-simplify", choices=("none", "basic", "pickup"), default="pickup"
    )
    parser.add_argument("--certify", action="store_true")
    parser.add_argument("--isometry-signatures", action="store_true")
    parser.add_argument("--dry-run", action="store_true", help="skip every SnapPy computation")
    parser.add_argument("--retry-errors", action="store_true")
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args(argv)

    if args.p_min < 2 or args.p_max < args.p_min:
        parser.error("require 2 <= p-min <= p-max")
    if args.q_min < 1:
        parser.error("require q-min >= 1")
    if args.s_min < 1 or args.s_max < args.s_min:
        parser.error("require 1 <= s-min <= s-max")
    if args.workers < 1 or args.checkpoint_every < 1:
        parser.error("workers and checkpoint-every must be positive")
    if args.volume_tolerance <= 0 or args.volume_epsilon <= 0:
        parser.error("volume tolerances must be positive")
    if args.certify and not _has_sage():
        parser.error("--certify requires running this script with Sage")
    if args.output_dir is None:
        suffix = "certified" if args.certify else "numerical"
        args.output_dir = Path(__file__).resolve().parent / "output" / (
            f"generalized_p{args.p_min}-{args.p_max}_s{args.s_min}-{args.s_max}_"
            f"{args.region}_{args.volume_scope}_{suffix}"
        )
    return args


def _configuration(args: argparse.Namespace) -> dict:
    names = (
        "p_min",
        "p_max",
        "q_min",
        "q_max",
        "r_min",
        "r_max",
        "s_min",
        "s_max",
        "region",
        "volume_scope",
        "volume_tolerance",
        "volume_epsilon",
        "retries",
        "bits_prec",
        "diagram_simplify",
        "certify",
        "isometry_signatures",
        "dry_run",
    )
    return {name: getattr(args, name) for name in names}


def run(args: argparse.Namespace) -> dict[str, int | float | str]:
    output_dir: Path = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    config_path = output_dir / "config.json"
    configuration = _configuration(args)
    if config_path.exists():
        previous = json.loads(config_path.read_text(encoding="utf-8"))
        if previous != configuration:
            raise RuntimeError(
                f"configuration differs from {config_path}; choose another --output-dir"
            )
    else:
        _write_json(config_path, configuration)

    started = time.monotonic()
    parameters = list(parameter_range(args))
    if not parameters:
        raise ValueError("the requested parameter range is empty")
    print(f"Computing exact Alexander polynomials for {len(parameters)} knots...", flush=True)

    records_by_parameter: dict[Parameters, dict[str, str]] = {}
    coefficient_cache: dict[str, tuple[int, ...]] = {}
    members_by_digest: dict[str, list[Parameters]] = defaultdict(list)
    for index, item in enumerate(parameters, 1):
        coefficients = alexander_coefficients(item)
        digest = alexander_digest(coefficients)
        cached = coefficient_cache.setdefault(digest, coefficients)
        if cached != coefficients:
            raise ArithmeticError("SHA-256 collision between distinct polynomials")
        members_by_digest[digest].append(item)
        records_by_parameter[item] = _empty_record(item, coefficients)
        if index % 1000 == 0:
            print(f"  Alexander phase: {index}/{len(parameters)}", flush=True)

    results_path = output_dir / "knots.csv"
    previous_rows = _read_csv(results_path)
    for item, old in previous_rows.items():
        if item not in records_by_parameter:
            continue
        current = records_by_parameter[item]
        if old.get("alexander_sha256") != current["alexander_sha256"]:
            raise RuntimeError(f"stored Alexander polynomial changed for {item}")
        for field in CSV_FIELDS:
            if field in old and field not in (
                "alexander_sha256",
                "alexander_breadth",
                "alexander_nonzero_terms",
                "braid_strands",
                "braid_crossings_before_simplification",
                "known_type",
            ):
                current[field] = old[field]

    collision_digests = {
        digest for digest, members in members_by_digest.items() if len(members) > 1
    }
    targets: list[Parameters] = []
    for item in parameters:
        row = records_by_parameter[item]
        if row["known_type"]:
            continue
        requested = (
            args.volume_scope == "all"
            or row["alexander_sha256"] in collision_digests
        )
        if not requested or args.dry_run:
            if row["volume_status"] in ("pending", "not_requested"):
                row["volume_status"] = "not_requested"
            continue
        finished = row["volume_status"] not in ("pending", "not_requested", "")
        if finished and not (args.retry_errors and row["volume_status"] in ("error", "unresolved")):
            continue
        row["volume_status"] = "pending"
        targets.append(item)

    ordered_records = lambda: [records_by_parameter[item] for item in parameters]
    _write_csv(results_path, ordered_records())
    print(
        f"Unique polynomials: {len(members_by_digest)}; "
        f"collision groups: {len(collision_digests)}; pending volumes: {len(targets)}",
        flush=True,
    )

    base_options = {
        "retries": args.retries,
        "volume_epsilon": args.volume_epsilon,
        "bits_prec": args.bits_prec,
        "diagram_simplify": args.diagram_simplify,
        "certify": args.certify,
    }

    def geometry_task(item: Parameters) -> tuple[Parameters, dict]:
        options = dict(base_options)
        # Canonical signatures can be substantially more expensive than
        # volumes for large diagrams.  For a full-volume run they are needed
        # at this stage only on equal-Alexander groups; any accidental
        # equal-volume/different-Alexander candidates can be rerun separately.
        options["isometry_signatures"] = (
            args.isometry_signatures
            and records_by_parameter[item]["alexander_sha256"] in collision_digests
        )
        return item, options

    completed = 0
    if targets and args.workers == 1:
        iterator = (_compute_geometry(geometry_task(item)) for item in targets)
        for item, geometry in iterator:
            records_by_parameter[item].update(geometry)
            completed += 1
            if completed % args.checkpoint_every == 0 or completed == len(targets):
                _write_csv(results_path, ordered_records())
                print(f"  Volume phase: {completed}/{len(targets)}", flush=True)
    elif targets:
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            futures = {
                pool.submit(_compute_geometry, geometry_task(item)): item for item in targets
            }
            for future in as_completed(futures):
                item, geometry = future.result()
                records_by_parameter[item].update(geometry)
                completed += 1
                if completed % args.checkpoint_every == 0 or completed == len(targets):
                    _write_csv(results_path, ordered_records())
                    print(f"  Volume phase: {completed}/{len(targets)}", flush=True)

    summary = write_comparison_reports(
        output_dir,
        ordered_records(),
        coefficient_cache,
        args.volume_tolerance,
        args.volume_scope,
    )
    summary["elapsed_seconds"] = round(time.monotonic() - started, 3)
    _write_json(output_dir / "summary.json", summary)
    print(json.dumps(summary, indent=2, sort_keys=True), flush=True)
    print(f"Results: {output_dir}", flush=True)
    return summary


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    try:
        run(args)
    except KeyboardInterrupt:
        print("Interrupted; the latest CSV checkpoint can be resumed.", file=sys.stderr)
        raise SystemExit(130)


if __name__ == "__main__":
    main()

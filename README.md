# Computations for twisted torus knots in the complementary range

This repository contains the code and data accompanying

> Adnan and Kyungbae Park, *The Alexander polynomial of twisted torus knots
> in the complementary range*.

The parameter convention is \(T(p,q;r,s)\). The main positive-twist
experiment uses

\[
3\le p\le20,\qquad 0<q<p,\qquad \gcd(p,q)=1,\qquad
0<r<p+q,\qquad 1\le s\le5.
\]

The scripts compute the Alexander polynomial from the unified formula in the
paper, construct braid diagrams, ask SnapPy for hyperbolic volumes and
decorated link isometry signatures, and compute hat knot Floer homology with
Spherogram and Szabó's HFK Calculator.

## Contents

- `code/alexander_polynomial.py`: standalone implementation of the
  unified formula for all nonzero relatively prime integers \(p,q\),
  \(1<r<|p|+|q|\), and every integer \(s\).
- `code/generalized_ttk_volume_experiment.py`: Alexander-polynomial and
  hyperbolic-volume experiment, with checkpointing and optional interval
  certification.
- `code/generalized_ttk_hfk_experiment.py`: HFK comparison inside
  equal-Alexander groups.
- `code/check_fiberedness_complementary_box.py`: the 52
  complementary-range HFK computations reported in the paper.
- `code/check_alexander_degree.py`: finite verification of the breadth
  formula.
- `code/check_non_l_space_family.py`: coefficient check for the infinite
  non-\(L\)-space family.
- `code/check_fiberedness_end_families.py`: independent implementation
  of the Alexander formula and endpoint braid models used by the checks.
- `code/test_generalized_ttk_volume_experiment.py`: regression tests.
- `results/`: the machine-readable data used in the computational
  discussion.

The files in `results/volume-all-p3-20-s1-5` come from the all-volume
run. The files in
`results/volume-collisions-with-signatures-p3-20-s1-5` come from the
collision-only run with decorated link isometry signatures enabled. The
directory `results/hfk-complementary-p3-6` contains all 52 HFK records
used as evidence for the fiberedness conjecture.

## Software

The reported computations used:

- Python 3.9.6;
- SnapPy 3.2;
- Spherogram 2.3;
- `knot_floer_homology` 1.2.2;
- SageMath 10.7 for optional interval-certified volume computations.

For an ordinary Python environment:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

SageMath is not required for the numerical experiments or the regression
tests. It is required only when `--certify` is requested.

## Quick verification

Run the regression tests:

```bash
python -m unittest discover -s code -p 'test_*.py'
```

Verify the Alexander-breadth formula on the default finite range:

```bash
python code/check_alexander_degree.py
```

Verify the coefficient obstruction for \(2\le k\le100\):

```bash
python code/check_non_l_space_family.py --k-max 100
```

These three commands do not perform a large volume or HFK search.

## Direct Alexander-polynomial calculation

The main reusable functions are

```python
# Run from the repository's code/ directory, or add code/ to sys.path.
from alexander_polynomial import (
    alexander_coefficients,
    alexander_polynomial,
)

coefficients = alexander_coefficients(p, q, r, s)
polynomial = alexander_polynomial(p, q, r, s)
```

Here \(p,q\) may be any nonzero relatively prime integers,
\(1<r<|p|+|q|\), and \(s\) may be any integer. The coefficient tuple is in
ascending exponent order, normalized to have lowest exponent zero and
positive leading coefficient. For example:

```bash
python code/alexander_polynomial.py 5 -3 6 2 --coefficients
```

The program uses

\[
T(p,q;r,s)\sim T(q,p;r,s),\qquad
T(-p,-q;r,s)\sim T(p,q;r,s),
\]

and

\[
\overline{T(p,q;r,s)}\sim T(p,-q;r,-s)
\]

to reduce to \(P>Q>0\). Since the Alexander polynomial is invariant under
mirroring, opposite signs of \(p,q\) replace \(s\) by \(-s\). The block
formula is then evaluated exactly over the integers. The implementation also
covers \(Q=1\), the boundary values \(r=P,Q\), and \(s=0\).

## Volume computations

A small smoke test of the generalized range is:

```bash
python code/generalized_ttk_volume_experiment.py \
  --p-max 8 --s-max 2 --region all --volume-scope collisions \
  --isometry-signatures --workers 2 \
  --output-dir output/smoke-volume
```

The full all-volume run underlying the count of 7,504 numerical hyperbolic
parameter choices is:

```bash
python code/generalized_ttk_volume_experiment.py \
  --p-max 20 --s-max 5 --region all --volume-scope all \
  --workers 4 --output-dir output/volume-all-p3-20-s1-5
```

The collision-only run with decorated link isometry signatures is:

```bash
python code/generalized_ttk_volume_experiment.py \
  --p-max 20 --s-max 5 --region all --volume-scope collisions \
  --isometry-signatures --workers 4 \
  --output-dir output/volume-collisions-with-signatures-p3-20-s1-5
```

The program scans 12,435 parameter quadruples because it retains the
elementary boundary cases \(r=1\) and \(r=p\) for consistency checks. The
range addressed by the two non-boundary formulas in the paper consists of
7,450 choices with \(1<r<p\) and 3,725 choices with \(p<r<p+q\).

For interval-certified hyperbolicity, run the program through Sage:

```bash
DOT_SAGE=/tmp/sage-ttk sage \
  code/generalized_ttk_volume_experiment.py \
  --p-max 12 --s-max 3 --region all --volume-scope collisions \
  --certify --bits-prec 100 --output-dir output/certified-check
```

Without `--certify`, a row marked `hyperbolic_numerical`
means that SnapPy found a positively oriented numerical geometric
triangulation. It is computational evidence rather than a proof of
hyperbolicity. Failed triangulations are recorded as `unresolved`,
not as nonhyperbolic.

## Knot Floer homology

Recompute the 52 cases reported in the paper:

```bash
python code/check_fiberedness_complementary_box.py \
  --p-max 6 --twists -2 -1 1 2 --workers 4 \
  --output output/complementary-hfk-p3-6.csv
```

To study HFK inside Alexander collision groups after a volume run:

```bash
python code/generalized_ttk_hfk_experiment.py \
  output/volume-collisions-with-signatures-p3-20-s1-5 \
  --selection unresolved --max-crossings 120 --timeout 60 --workers 4
```

Every accepted HFK computation is checked by comparing its graded Euler
characteristic with the unified Alexander-polynomial formula. Each HFK
calculation runs in a separate subprocess with a hard timeout.

## Reproducing and interpreting the stored data

Each volume-result directory contains:

- `config.json`: the exact command-line configuration;
- `knots.csv`: one row per parameter choice;
- `alexander_collision_groups.csv`: nontrivial equal-polynomial groups;
- `same_alexander_different_volume.csv`: possible counterexamples;
- `same_volume_different_alexander.csv`: converse volume collisions,
  when every volume was requested;
- `same_alexander_different_signature.csv`: signature discrepancies,
  when signatures were requested;
- `summary.json`: final counts.

Alexander polynomials are normalized to have lowest exponent zero and positive
leading coefficient. Equality is ultimately checked using the complete
integer coefficient tuple; SHA-256 is used only as an indexing key. Volume
agreement is tested with tolerance \(10^{-8}\).

Numerical equality of volumes does not prove isotopy. Likewise, a finite
search cannot prove the conjectures stated in the paper.

## License

The code is released under the MIT License; see `LICENSE`. The
computed CSV and JSON data may be reused with attribution to the accompanying
paper and this repository.

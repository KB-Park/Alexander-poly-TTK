# Computational review of the Alexander-polynomial conjecture

Date: 2026-09-10

## 1. Scope and normalization

The main positive-twist search used

\[
3\le p\le20,\qquad 0<q<p,\qquad \gcd(p,q)=1,
\qquad 0<r<p+q,\qquad 1\le s\le5.
\]

Alexander polynomials were normalized to have lowest exponent zero and
positive leading coefficient. Equality was checked using the complete dense
integer coefficient tuple; SHA-256 was only used as an indexing key, with an
exact tuple comparison guarding every hash collision group.

For `r<=p`, SnapPy received the standard `p`-strand braid. For `p<r<p+q`,
it received the `(p+q)`-strand Lorenz shuffle braid for `T(p,q)` followed by
the full twist on the first `r` strands. Small examples were independently
checked against the unified Alexander formula and against the previously
constructed destabilized braids.

Numerical hyperbolicity was accepted only when SnapPy found a triangulation
with all tetrahedra positively oriented. Failed triangulations were recorded
as `unresolved`, not as non-hyperbolic. Known torus-knot boundary cases were
removed before volume comparison. Decorated link isometry signatures were
computed in addition to volumes.

## 2. Positive-twist conjecture search

The complete parameter box contains 12,435 parameter quadruples and 10,156
distinct Alexander polynomials. There are 650 nontrivial equal-Alexander
groups containing 2,929 parameter quadruples.

Among the collision members:

- 2,281 are elementary known torus-knot presentations;
- 586 received numerical hyperbolic structures;
- 62 remained unresolved by the selected SnapPy triangulation attempts.

There are 258 collision groups with at least two numerical hyperbolic members.
In every one of these groups:

- all volumes agree within `10^-8` (in fact, the stored floating-point values
  agree at the displayed precision);
- all computed decorated link isometry signatures agree;
- no same-Alexander/different-volume or
  same-Alexander/different-signature pair occurs.

Thus this box supplies substantially stronger evidence than volume agreement
alone: the computed peripheral hyperbolic invariants agree as well. It is still
a finite numerical experiment, not a proof of isotopy or of the conjecture.

An additional all-volume run computed every non-elementary candidate in the
same parameter box. Its final status counts are:

- 7,504 numerical hyperbolic knots;
- 2,655 elementary known torus-knot presentations;
- 2,276 unresolved cases.

Among the 7,504 numerical hyperbolic knots, there is no pair whose volumes
agree within `10^-8` but whose Alexander polynomials differ. Thus the converse
observation from the earlier paper also persists throughout the successfully
solved part of this generalized box.

As a certified spot check, both `T(4,3;5,1)` and `T(6,5;3,1)` are hyperbolic
with 120-bit interval-certified volume

\[
4.0261201267429207\ldots
\]

and both have decorated link isometry signature
`fLAMcaccdeejsjqqn_bBRLrM`.

## 3. Parameter redundancy in the new range

Exactly 26 members with `p<r<p+q` occur in an Alexander collision group in
the selected box. Every one has `q=p-1` and `r=p+1`, and every collision is
of the form

\[
T(p,p-1;p+1,s)
\quad\longleftrightarrow\quad
T((p+1)s+1,p+1;p-1,1).
\]

All 20 new-range members for which SnapPy found a hyperbolic structure have
the same volume and decorated link signature as the displayed partner. The
six unresolved new-range members also fit the same parameter formula. The
four small unresolved pairs

\[
T(3,2;4,s),\quad T(4s+1,4;2,1),\qquad 1\le s\le4,
\]

have identical computed bigraded hat HFK. This uniform pattern strongly
suggests that the parameter transformation should be proved by an explicit
braid isotopy. A proof would remove a systematic source of parameter
nonuniqueness before formulating injectivity of the Alexander polynomial.

## 4. Same Alexander polynomial and different HFK

No such pair was found in the positive-twist search performed so far. For the
collision groups containing unresolved volume cases, HFK was attempted on all
members with at most 300 raw crossings:

- 227 candidates in 47 Alexander collision groups;
- 211 successful HFK computations;
- 41 groups with at least two successful computations;
- 12 failures caused by the HFK calculator's girth bound;
- 4 timeouts;
- zero same-Alexander/different-HFK pairs.

There is, however, a clean negative-twist example. The already known
same-Alexander/different-volume pair

\[
K_1=T(4,3;2,-5),\qquad K_2=T(6,5;4,-2)
\]

has common normalized Alexander coefficient tuple

\[
(1,-2,2,-1,1,-1,2,-2,1),
\]

but different bigraded hat knot Floer homology. More simply,

| knot | genus | tau | total HFK rank | certified volume |
|---|---:|---:|---:|---:|
| `T(4,3;2,-5)` | 4 | -2 | 13 | 4.854663338650145... |
| `T(6,5;4,-2)` | 4 | 0 | 15 | 4.553766624919092... |

Both knots are fibered and non-L-space knots. Their hyperbolicity and volumes
were verified with 120-bit interval arithmetic. Hence HFK, and already its
tau invariant, strictly refines the Alexander polynomial on negatively
twisted torus knots.

## 5. Other application computations

### Fiberedness, genus, and L-space status

For every one of the 52 parameter choices satisfying

\[
3\le p\le6,\qquad 2\le q<p,\qquad \gcd(p,q)=1,
\qquad p<r<p+q,
\qquad s\in\{-2,-1,1,2\},
\]

the HFK computation completed. Every knot was fibered, its top Alexander
grading had rank one, its HFK genus was half the Alexander breadth, and its
graded Euler characteristic agreed with the unified Alexander-polynomial
formula. The complete records are in
`complementary_hfk_p3-6_s-2--1_1-2.csv`.

For `3<=p<=8`, `1<=s<=2`, and raw braid crossing number at most 120, 271
HFK computations completed, including 35 knots in the new `p<r<p+q` range.

- every computed knot is fibered;
- every HFK genus equals half the Alexander breadth;
- 14 knots are non-L-space knots, all in the old `r<p` range;
- all 35 computed new-range knots are L-space knots in this small sample.

The 14 non-L-space examples are recorded in
`generalized_p3-8_s1-2_all_all_numerical/hfk_all-knots_crossings120/non_l_space_knots.csv`.
In each of these 14 cases, the Alexander coefficients already violate the
necessary alternating `0, +/-1` form for an L-space knot, so this particular
sample does not yet give a stronger HFK-only L-space obstruction.

### Volume growth versus Alexander degree

In the all-volume search with `3<=p<=10` and `1<=s<=5`, there are 103 fixed
triples `(p,q,r)` for which all five twists received numerical hyperbolic
structures. Volume is strictly increasing with `s=1,...,5` for all 103
triples, including all 37 such triples in the new range.

For the illustrative new-range family `T(4,3;5,s)`:

| s | Alexander genus | volume |
|---:|---:|---:|
| 1 | 13 | 4.026120126743 |
| 2 | 23 | 4.465921490293 |
| 3 | 33 | 4.592107095247 |
| 5 | 53 | 4.671283872468 |
| 10 | 103 | 4.710560847064 |
| 15 | 153 | 4.718596375812 |
| 20 | 203 | 4.721513105362 |

This illustrates a useful geometric application: genus and Alexander breadth
grow linearly with the number of twists while volume appears to converge to
the augmented-link volume. The data suggest studying eventual or full
monotonicity of volume under positive twisting, as well as effective
Dehn-filling error estimates.

## 6. Suggested mathematical next steps

1. Prove the parameter transformation in Section 3 by braid conjugation,
   exchange moves, or a flat-torus isotopy.
2. Certify all 586 numerical hyperbolic collision members with interval
   arithmetic; the present full-box result is numerical except for spot checks.
3. Resolve the 62 difficult geometries using alternative triangulations or
   recognized torus/satellite classifications.
4. Extend the negative-twist HFK collision search, using tau, epsilon,
   fiberedness, and total HFK rank as inexpensive first-stage filters.
5. Construct the augmented two-component link explicitly and compare its
   volume with the `s -> infinity` data for fixed `(p,q,r)`.

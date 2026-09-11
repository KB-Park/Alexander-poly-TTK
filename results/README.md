# Stored computational results

The stored files are unmodified outputs of the scripts in `../code`,
except that they have been grouped into descriptive directories for
publication.

## `volume-all-p3-20-s1-5`

This is the completed numerical run with `--volume-scope all` and
without decorated isometry signatures. Its summary records:

- 12,435 total parameter choices;
- 10,156 distinct normalized Alexander polynomials;
- 7,504 numerical hyperbolic structures;
- 2,655 elementary torus-knot presentations;
- 2,276 unresolved cases;
- no equal-volume/different-Alexander pair at tolerance \(10^{-8}\).

## `volume-collisions-with-signatures-p3-20-s1-5`

This run restricted volume and signature computation to the 650 nontrivial
Alexander collision groups. There are 258 groups with at least two numerical
hyperbolic members. No same-Alexander/different-volume or
same-Alexander/different-signature pair was found.

## `hfk-complementary-p3-6`

These are the 52 parameter choices

\[
3\le p\le6,\quad 2\le q<p,\quad \gcd(p,q)=1,\quad
p<r<p+q,\quad s\in\{-2,-1,1,2\}.
\]

All 52 HFK calculations completed. Each knot had rank one in its top
Alexander grading, and the HFK Euler characteristic agreed with the
Alexander-polynomial formula.

See `computational-report.md` for a broader account of the exploratory
computations and their limitations.

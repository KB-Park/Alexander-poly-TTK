# The two-parameter family T(2n+1, 2n; 3n, s)

The Alexander formula gives several uniform non-L-space families.  All
Alexander polynomials below are normalized to have lowest exponent zero and
positive constant term.

## Negative twists

Let `u >= 1`.  For every `n >= 3`,

\[
  [t^{4(3u-2)n}]\,
  \Delta_{T(2n+1,2n;3n,-u)}(t)=2.
\]

Consequently,

\[
  \{T(2n+1,2n;3n,-u): n\geq3,\ u\geq1\}
\]

is a two-parameter family of non-L-space knots.

To verify the coefficient, put

\[
 d=(3u-2)n,\quad
 g=2(3u-1)n-1,\quad
 \delta=3(2u-1)n.
\]

The block data in the paper's formula are independent of `u`:

\[
 j_*=i_*=n,
\]

\[
 (L_0,\ldots,L_{2n-1})
 =(2,\underbrace{1,\ldots,1}_{n},
      \underbrace{2,\ldots,2}_{n-1}),
\]

and

\[
 (M_0,\ldots,M_{2n})
 =(\underbrace{2,\ldots,2}_{n-1},
   \underbrace{1,\ldots,1}_{n+2}).
\]

Let

\[
 N(t)=Z_a(t)X_c(t)-X_a(t)Z_c(t)
     =(1+t+\cdots+t^{3n-1})\Delta(t)
\]

up to a unit.  Its lowest exponent is

\[
 E_0=(6-12u)n^2+(4-3u)n-1.
\]

After shifting by \(-E_0\), the terms of degree at most `4d` that can
affect the residue class `n` modulo `3n` are obtained from

\[
 hd\quad(0\leq h\leq4),\qquad g,\qquad g+d,
 \qquad \delta,\qquad \delta+d-1.
\]

The first three groups arise from \(Z_aX_c\), while the last two arise
from \(X_aZ_c\) and are subtracted.  In the coefficient recurrence

\[
 a_j=b_j-b_{j-1}+a_{j-3n},
\]

only the terms at `d`, `4d`, `g`, and `delta+d-1` contribute to
\(a_{4d}\).  Their contributions are respectively

\[
 1,\quad 1,\quad -1,\quad 1.
\]

Thus \(a_{4d}=2\).  When `u=1`, the middle two terms cancel at the
numerator level and the same conclusion follows from the terms at `d` and
`4d`.  The same recurrence also shows that all coefficients in degrees
strictly below `4d` belong to `0`, `1`, and `-1`.  Thus `4d` is the first
degree at which a coefficient has absolute value greater than one.

The breadth is

\[
 \operatorname{br}\Delta(t)
 =(9u-4)n^2-(3u+4)n+2.
\]

These two quantities also distinguish all members of the negative-twist
family.  Indeed, define

\[
 \nu(\Delta)=\min\{j:|[t^j]\Delta(t)|>1\}.
\]

For this family, \(\nu=4n(3u-2)\).  Put \(x=\nu/4\).  Rewriting the
breadth in terms of `x` gives

\[
 \operatorname{br}\Delta=(3n-1)x+2n^2-6n+2.
\]

For fixed `x`, the right side increases strictly with `n`, since its value
at `n+1` minus its value at `n` is

\[
 3x+4n-4>0.
\]

Thus the polynomial determines `n` uniquely, and then
\(u=(x/n+2)/3\).  Hence distinct pairs `(n,u)` in the negative-twist
family have distinct Alexander polynomials and therefore represent
distinct knots.

## Positive twists

Direct substitution in the same block formula gives two further uniform
obstructions:

\[
 [t^{(15u+9)n+3}]\,
 \Delta_{T(2n+1,2n;3n,u)}(t)=-2
 \quad(n\geq4,\ u\geq2),
\]

and

\[
 [t^{29n+4}]\,
 \Delta_{T(2n+1,2n;3n,1)}(t)=-2
 \quad(n\geq5).
\]

For reference, if

\[
 N(t)=\sum_j b_jt^{E_++j},\qquad
 E_+=(3u+2)n^2+(6u+3)n,
\]

then, for `u >= 2` and `J=(15u+9)n+3`, substitution of the displayed
`L`- and `M`-blocks shows that the only nonzero differences
\(b_j-b_{j-1}\) with \(j\equiv J\pmod{3n}\) and `j <= J` occur at

\[
 j=J-3n\quad\text{and}\quad j=J,
\]

and both differences equal `-1`.  Hence the quotient recurrence gives
\(a_J=-2\).  For `u=1`, setting `J=29n+4` gives directly
\(b_J-b_{J-1}=-2\), while all earlier differences in that residue class
vanish.  In both cases, the recurrence shows that this is the first
coefficient whose absolute value exceeds one.

The positive-twist breadth is

\[
 \operatorname{br}\Delta(t)
 =(9u+4)n^2-(3u+2)n.
\]

The positive-twist knots covered above are also pairwise distinct.  For
`u >= 2`, let \(\nu=(15u+9)n+3\) be the first obstructing degree and put
\(x=(\nu-3)/3=n(5u+3)\).  Then

\[
 \operatorname{br}\Delta
 =\frac{(9n-3)x-7n^2-n}{5}.
\]

For fixed `x`, this expression is strictly increasing through its
admissible values of `n`; the numerator of its forward difference is
\(9x-14n-8>0\).  Thus `n`, and then \(u=(x/n-3)/5\), are determined by
the polynomial.  In the `u=1` subfamily, \(\nu=29n+4\) determines `n`
directly.  Nor can a member with `u=1` coincide with one with `u >= 2`:
after equating both \(\nu\) and the breadth, one obtains

\[
 m^2(2925u^2-4059u-2311)+m(143-42u)+158=0,
\]

which is positive for \(m\geq4\) and \(u\geq2\).  Finally, a negative-
twist member cannot have the same Alexander polynomial as a positive-
twist member, because the first obstructing coefficient is respectively
`+2` and `-2`.

## Exceptional ranges

- `s=0` gives the torus knot \(T(2n+1,2n)\), hence an L-space knot.
- For positive twists, the above coefficient obstruction does not occur
  when `n=2`, when `n=3`, or for `(n,s)=(4,1)`.  Exact computations for
  `n=2,3` and `1 <= s <= 100`, as well as for `(n,s)=(4,1)`, give only
  coefficients `0`, `1`, and `-1`; this necessary Alexander condition
  alone does not decide their L-space status.
- For `s=-1`, `n=1` gives the unknot, while `n=2` gives the mirror of
  \(T(2,9)=9_1\), a negative L-space knot.  More generally, the displayed
  negative-twist obstruction begins at `n=3`.

Run the exact finite check with, for example,

```bash
python code/check_non_l_space_2n_family.py \
  --n-min 2 --n-max 30 --s-min -10 --s-max 10
```

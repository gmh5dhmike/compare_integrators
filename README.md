# Compare different quadrature rules for integration


Numerical Integration Error Assessment
(Landau §6.2.5 — Implementation and Error Analysis)
Overview

This project implements and compares three numerical integration methods:

Composite Trapezoid Rule

Composite Simpson’s Rule

Gauss–Legendre Quadrature using a custom high-precision implementation (HighPrecisionGaussInt) that computes Legendre polynomial roots and weights with decimal precision (~40 digits)

Two experiments were performed:

Error scaling on a smooth test function – results in Errors.png

Performance on a intentionally hard-to-integrate function – results in BadErrors.png

All code is in:

gqconstants.py — high-precision Gaussian quadrature

errors.py — smooth-function error study

bad_errors.py — pathological function error study

1. Smooth Test Function Error Study (Errors.png)

For the main convergence test, I used:

f(x)=sin⁡(x),x∈[0,π],
f(x)=sin(x),x∈[0,π],

which has the exact integral:

∫0πsin⁡(x) dx=2.
∫
0
π
	​

sin(x)dx=2.
Methods Compared

Trapezoid and Simpson’s rule were run for

N=10,20,40,80,160,320,640,1280.
N=10,20,40,80,160,320,640,1280.

Gaussian quadrature was evaluated with

n=2,4,6,8,10,12,16
n=2,4,6,8,10,12,16

high-precision Gauss–Legendre points.

The resulting absolute errors (plotted on log–log axes) appear in Errors.png.

Observed Behavior

Trapezoid Rule:
Converges as 
O(h2)
O(h
2
). The log–log slope is approximately –2.

Simpson’s Rule:
Converges as 
O(h4)
O(h
4
). The slope is approximately –4.

Gaussian Quadrature:
Errors decrease extremely rapidly with number of points.
For a smooth integrand like 
sin⁡x
sinx, even 6–10 Gauss points give very small errors, and round-off eventually dominates.

This matches the theoretical behavior described in Landau’s Fig. 6.3.

2. Hard-to-Integrate Function (BadErrors.png)

To stress-test the algorithms, I constructed a narrow spike on 
[0,1]
[0,1]:

f(x)={1	if ∣x−0.5∣<ε,
0	otherwise,ε=10−3.
f(x)={
1
0
	​

if ∣x−0.5∣<ε,
otherwise,
	​

ε=10
−3
.

The exact value of the integral is:

∫01f(x) dx=2ε=0.002.
∫
0
1
	​

f(x)dx=2ε=0.002.
Why This Function Is Hard

This integrand violates the smoothness assumptions used by all quadrature rules:

The spike is extremely localized (width 
2×10−3
2×10
−3
).

Most evaluation points miss it completely.

Trapezoid/Simple rules use evenly spaced points.

Gauss–Legendre rules use fixed nodes determined globally by the Legendre polynomial roots.

If no node falls inside 
[0.5−ε,0.5+ε]
[0.5−ε,0.5+ε], every method returns ~0, producing an absolute error ≈ 0.002 and no convergence.

This failure is shown clearly in BadErrors.png, where the error curves flatten rather than dropping with 
N
N.

3. Why the Algorithms Struggle

All three integration methods assume:

The integrand is reasonably smooth, and

It can be approximated well by low-order polynomials over the subintervals (trapezoid/Simpson) or by global orthogonal polynomials (Gaussian).

The narrow spike breaks these assumptions:

1. Localization

The important contribution to the integral exists only in a region too small for fixed sampling grids.

2. Sampling Failure

If none of the rule’s evaluation points land inside the spike, the method cannot detect it—it effectively integrates the zero function.

3. Lack of Adaptivity

None of these rules analyzes the behavior of the integrand and adjusts sampling accordingly.

Thus even Gaussian quadrature, which is excellent for smooth functions, can fail completely for sharply localized ones.

4. How to Improve the Tough Calculation

Even without external libraries, several strategies can dramatically improve accuracy:

A. Domain Decomposition (Most Effective)

Split the integration interval around the spike:

∫01f(x) dx=∫00.49f(x) dx+∫0.490.51f(x) dx+∫0.511f(x) dx.
∫
0
1
	​

f(x)dx=∫
0
0.49
	​

f(x)dx+∫
0.49
0.51
	​

f(x)dx+∫
0.51
1
	​

f(x)dx.

Then:

Use very fine resolution or a high-order Gauss rule only on the middle subinterval.

Keep coarse sampling elsewhere.

This guarantees sampling inside the spike.

B. Manual “Adaptive” Refinement

A simple refinement strategy using only our existing tools:

Compute trapezoid/Simpson with some 
N
N.

Double 
N
N.

Compare results; if they differ by more than a tolerance, subdivide that region further.

This directs resolution to difficult regions without uniformly increasing work.

C. Higher-Order Gaussian Quadrature on Smaller Intervals

Gaussian quadrature becomes extremely powerful when its nodes fall inside the feature.
Combining it with domain splitting yields excellent accuracy even for pathological integrands.

D. Change of Variables (for Other Hard Functions)

For endpoint singularities (not the spike example), a substitution like:

x=t2
x=t
2

can turn an integrable singularity such as 
1/x
1/
x
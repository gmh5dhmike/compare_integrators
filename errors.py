import math
from decimal import Decimal
import numpy as np
import matplotlib.pyplot as plt

from gqconstants import HighPrecisionGaussInt


# ---------------------------------------------------------
# Function to integrate 
# ---------------------------------------------------------

def f(x: float) -> float:
    return math.sin(x)


a, b = 0.0, math.pi
I_exact = 2.0  # analytical value of ∫0^π sin(x) dx


# ---------------------------------------------------------
# Double-precision integrators: trapezoid + Simpson
# ---------------------------------------------------------

def trap(f, a: float, b: float, N: int) -> float:
    """Composite trapezoid rule."""
    h = (b - a) / N
    s = 0.5 * (f(a) + f(b))
    for i in range(1, N):
        s += f(a + i * h)
    return h * s


def simpson(f, a: float, b: float, N: int) -> float:
    """Composite Simpson's rule (N must be even)."""
    if N % 2 != 0:
        raise ValueError("Simpson's rule requires even N.")
    h = (b - a) / N
    s = f(a) + f(b)
    for i in range(1, N):
        x = a + i * h
        if i % 2 == 0:
            s += 2 * f(x)
        else:
            s += 4 * f(x)
    return s * h / 3.0


# ---------------------------------------------------------
# Gaussian quadrature using high-precision weights
# ---------------------------------------------------------

def gauss_integral(npoints: int, f, a: float, b: float) -> float:
    """
    Use high-precision Gauss-Legendre quadrature with npoints nodes.
    Returns a standard Python float for plotting.
    """
    npoints = int(npoints)  # convert numpy.int64 → python int
    gauss_hp = HighPrecisionGaussInt(npoints=npoints, precision=40)
    I_dec = gauss_hp.integ(f, a, b)
    return float(I_dec)


# ---------------------------------------------------------
# Main error comparison
# ---------------------------------------------------------

def main():

    # N values for trapezoid and Simpson
    Ns = np.array([10, 20, 40, 80, 160, 320, 640, 1280])

    trap_err = []
    simp_err = []

    for N in Ns:
        # Trapezoid
        I_tr = trap(f, a, b, int(N))
        trap_err.append(abs(I_tr - I_exact))

        # Simpson (N must be even)
        N_simp = N if N % 2 == 0 else N + 1
        I_si = simpson(f, a, b, int(N_simp))
        simp_err.append(abs(I_si - I_exact))

    # Gaussian quadrature: vary node count
    gauss_points = np.array([2, 4, 6, 8, 10, 12, 16])
    gauss_err = []

    for npts in gauss_points:
        I_g = gauss_integral(int(npts), f, a, b)
        gauss_err.append(abs(I_g - I_exact))

    # -----------------------------------------------------
    # Create the log–log plot
    # -----------------------------------------------------

    plt.figure()

    plt.loglog(Ns, trap_err, "o-", label="Trapezoid")
    plt.loglog(Ns, simp_err, "s-", label="Simpson's Rule")
    plt.loglog(gauss_points, gauss_err, "^-", label="Gauss–Legendre")

    plt.xlabel("N (subintervals or quadrature points)")
    plt.ylabel("Absolute error")
    plt.title("Error Comparison of Numerical Integration Methods")
    plt.grid(True, which="both", linestyle="--", alpha=0.5)
    plt.legend()

    plt.savefig("Errors.png", dpi=200, bbox_inches="tight")
    print("Saved Errors.png")


if __name__ == "__main__":
    main()


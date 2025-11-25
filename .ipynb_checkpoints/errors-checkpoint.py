import math
from decimal import Decimal
import numpy as np
import matplotlib.pyplot as plt

from gqconstants import HighPrecisionGaussInt


# -----------------------------
# Test function and exact value
# -----------------------------

def f(x: float) -> float:
    """Function to integrate."""
    return math.sin(x)


a, b = 0.0, math.pi
I_exact = 2.0  # ∫_0^π sin(x) dx = 2


# -----------------------------
# Basic double-precision integrators
# -----------------------------

def trap(f, a: float, b: float, N: int) -> float:
    """Composite trapezoidal rule with N subintervals."""
    h = (b - a) / N
    s = 0.5 * (f(a) + f(b))
    for i in range(1, N):
        s += f(a + i * h)
    return h * s


def simpson(f, a: float, b: float, N: int) -> float:
    """
    Composite Simpson's rule with N subintervals.
    N must be even.
    """
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


# -----------------------------------------------
# Gaussian quadrature using HighPrecisionGaussInt
# -----------------------------------------------

def gauss_integral(npoints: int, f, a: float, b: float) -> float:
    """
    Use high-precision Gauss-Legendre with npoints nodes,
    but return a double-precision float.
    """
    gauss_hp = HighPrecisionGaussInt(npoints=npoints, precision=40)

    # We could call gauss_hp.integ directly, but that returns a Decimal.
    # Wrap once to get a float.
    I_dec = gauss_hp.integ(f, a, b)
    return float(I_dec)


# -----------------------------
# Error study and plot
# -----------------------------

def main():
    # N for trapezoid / Simpson
    Ns = np.array([10, 20, 40, 80, 160, 320, 640, 1280])

    trap_err = []
    simp_err = []

    for N in Ns:
        # Trapezoid
        I_trap = trap(f, a, b, N)
        trap_err.append(abs(I_trap - I_exact))

        # Simpson (ensure even N)
        N_simp = N if N % 2 == 0 else N + 1
        I_simp = simpson(f, a, b, N_simp)
        simp_err.append(abs(I_simp - I_exact))

    # Gaussian: vary number of points
    gauss_points = np.array([2, 4, 6, 8, 10, 12, 16])
    gauss_err = []

    for npts in gauss_points:
        I_g = gauss_integral(npts, f, a, b)
        gauss_err.append(abs(I_g - I_exact))

    # Convert to numpy arrays for plotting
    trap_err = np.array(trap_err)
    simp_err = np.array(simp_err)
    gauss_err = np.array(gauss_err)

    # ------------
    # Make plot
    # ------------

    plt.figure()

    # Trapezoid and Simpson vs N (number of subintervals)
    plt.loglog(Ns, trap_err, "o-", label="Trapezoid")
    plt.loglog(Ns, simp_err, "s-", label="Simpson")

    # Gaussian vs number of points
    plt.loglog(gauss_points, gauss_err, "^-", label="Gauss-Legendre")

    plt.xlabel("N / number of points")
    plt.ylabel("Absolute error")
    plt.legend()
    plt.grid(True, which="both", linestyle="--", alpha=0.5)

    plt.title("Error comparison of integration methods")
    plt.savefig("Errors.png", dpi=200, bbox_inches="tight")
    print("Saved Errors.png")


if __name__ == "__main__":
    main()

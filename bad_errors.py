import numpy as np
import matplotlib.pyplot as plt
from decimal import Decimal
from gqconstants import HighPrecisionGaussInt


# ---------------------------------------------------------
# Hard-to-integrate function: narrow top-hat spike
# ---------------------------------------------------------
# On [0, 1]:
#   f(x) = 1  if |x - 0.5| < eps
#          0  otherwise
# Exact integral = 2 * eps

eps = 1e-3
a, b = 0.0, 1.0
I_exact = 2.0 * eps


def bad_f(x: float) -> float:
    if abs(x - 0.5) < eps:
        return 1.0
    else:
        return 0.0


# ---------------------------------------------------------
# Double-precision trapezoid and Simpson
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
# Gaussian quadrature wrapper (high-precision weights)
# ---------------------------------------------------------

def gauss_integral(npoints: int, f, a: float, b: float) -> float:
    """
    Use high-precision Gauss-Legendre quadrature with npoints nodes.
    Returns a standard Python float.
    """
    npoints = int(npoints)
    gauss_hp = HighPrecisionGaussInt(npoints=npoints, precision=40)

    # Manual mapping to avoid numpy→Decimal issues
    a_dec = Decimal(str(a))
    b_dec = Decimal(str(b))
    c1 = (b_dec - a_dec) / Decimal(2)
    c2 = (b_dec + a_dec) / Decimal(2)

    total = Decimal(0)
    for w, xi in zip(gauss_hp.weight, gauss_hp.lroots):
        x_phys = c1 * xi + c2
        fx = Decimal(str(f(float(x_phys))))
        total += w * fx

    return float(c1 * total)


# ---------------------------------------------------------
# Main: compute errors and plot BadErrors.png
# ---------------------------------------------------------

def main():
    # N for trap/Simpson – go fairly large to show the failure
    Ns = np.array([10, 20, 40, 80, 160, 320, 640, 1280, 2560, 5120])

    trap_err = []
    simp_err = []

    for N in Ns:
        N = int(N)

        I_tr = trap(bad_f, a, b, N)
        trap_err.append(abs(I_tr - I_exact))

        N_simp = N if N % 2 == 0 else N + 1
        I_si = simpson(bad_f, a, b, N_simp)
        simp_err.append(abs(I_si - I_exact))

    # Gaussian: vary number of points
    gauss_points = np.array([2, 4, 6, 8, 10, 12, 16])
    gauss_err = []

    for npts in gauss_points:
        I_g = gauss_integral(int(npts), bad_f, a, b)
        gauss_err.append(abs(I_g - I_exact))

    trap_err = np.array(trap_err)
    simp_err = np.array(simp_err)
    gauss_err = np.array(gauss_err)

    # -----------------------------------------------------
    # Plot: BadErrors.png
    # -----------------------------------------------------

    plt.figure()

    plt.loglog(Ns, trap_err, "o-", label="Trapezoid (bad f)")
    plt.loglog(Ns, simp_err, "s-", label="Simpson (bad f)")
    plt.loglog(gauss_points, gauss_err, "^-", label="Gauss–Legendre (bad f)")

    plt.xlabel("N (subintervals / points)")
    plt.ylabel("Absolute error")
    plt.title("Errors for Hard-to-Integrate Function")
    plt.grid(True, which="both", linestyle="--", alpha=0.5)
    plt.legend()

    plt.savefig("BadErrors.png", dpi=200, bbox_inches="tight")
    print("Saved BadErrors.png")


if __name__ == "__main__":
    main()

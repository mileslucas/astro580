import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import fsolve

# Definitions
pi = np.pi
me = 9.1095e-28 # g
h = 6.6262e-27 # erg s
k = 1.3803e-16 # erg / K
chi1 = 3.939e-11 # erg
chi2 = 8.718e-11 # erg
NA = 6.02e23
rho = 1e-4 # g / cm^3
T = np.linspace(1e4, 2e5, 1000)


def eqns(p, T):
    z0, z1, z2 = p
    A = rho * NA * (z1 + 2 * z2) / 4
    C = (2 * pi * me * k * T) ** (3/2) / h**3

    eq1 = A * z1 / z0 - 4 * C * np.exp(-1 * chi1 / (k * T))
    eq2 = A * z2 / z1 - C * np.exp(-1 * chi2 / (k * T))
    eq3 = z0 + z1 + z2 - 1
    return (eq1, eq2, eq3)


if __name__ == '__main__':
    sols = np.asarray([fsolve(eqns, x0=(0.3, 0.3, 0.3), args=(t)) for t in T])
    z0, z1, z2 = sols.T
    ze = z1 + 2 * z2
    plt.style.use('seaborn')
    plt.title(rf'z1 and z2 for $\rho = {rho}$')

    plt.plot(T, z0, label='z0')
    plt.plot(T, z1, label='z1')
    plt.plot(T, z2, label='z2')
    plt.plot(T, ze, c='r', label='ze')

    plt.xlabel('T (K)')
    plt.ylim(0, 1)
    plt.legend()

    plt.show()

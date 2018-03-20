import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import least_squares

# Definitions
pi = np.pi
me = 9.1095e-28  # g
h = 6.6262e-27  # erg s
k = 1.3803e-16  # erg / K
chi1 = 3.939e-11  # erg
chi2 = 8.718e-11  # erg
NA = 6.02e23
rho = 1e-4  # g / cm^3
T = np.linspace(1e4, 2e5, 1000)


def eqns(p, T):
    z1, z2 = p

    A = rho * NA * (z1 + 2 * z2) / 4
    C = (2 * pi * me * k * T) ** (3 / 2) / h**3

    eq1 = A * z1 / (1 - z1 - z2) - 4 * C * np.exp(-1 * chi1 / (k * T))
    eq2 = A * z2 / z1 - C * np.exp(-1 * chi2 / (k * T))
    return (eq1, eq2)


if __name__ == '__main__':
    try:
        T, z0, z1, z2, ze = np.genfromtxt('model.csv', delimiter=',')
    except FileNotFoundError:
        z1 = []
        z2 = []
        for t in T:
            sol = least_squares(eqns, x0=(0.3, 0.3), args=(
                t, ), bounds=([0, 0], [1, 1]))
            z1.append(sol.x[0])
            z2.append(sol.x[1])
        z1, z2 = np.array([z1, z2])
        z0 = 1 - z1 - z2
        ze = z1 + 2 * z2
        np.savetxt('model.csv', [T, z0, z1, z2, ze], delimiter=',')

    # Find Crossing Points
    idx0 = np.argwhere(np.diff(np.sign(z0 - z1)) != 0)
    idx1 = np.argwhere(np.diff(np.sign(z1 - z2)) != 0)

    plt.style.use('seaborn')
    plt.title(rf'$z_i$ for $\rho = {rho:.2e}$')

    plt.semilogx(T, z0, label=r'$z_0$')
    plt.semilogx(T, z1, label=r'$z_1$')
    plt.semilogx(T, z2, label=r'$z_2$')

    plt.plot(T[idx0], z1[idx0], 'ko')
    plt.text(T[idx0] + 2000, z1[idx0], f'{int(T[idx0][0][0]):.2e}')
    plt.plot(T[idx1], z1[idx1], 'ko')
    plt.text(T[idx1] + 5000, z1[idx1], f'{int(T[idx1][0][0]):.2e}')
    # plt.semilogx(T, ze, c='r', label=r'z_e')

    plt.xlabel('T (K)')
    plt.ylabel('z')
    plt.legend()

    # plt.savefig('prob3.1.png')
    plt.show()

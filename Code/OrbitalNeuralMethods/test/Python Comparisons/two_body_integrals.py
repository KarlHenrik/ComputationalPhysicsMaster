def _trapz(f, x):
    n = len(x)
    delta_x = x[1] - x[0]
    val = 0

    for i in range(1, n):
        val += f[i - 1] + f[i]

    return 0.5 * val * delta_x

def _shielded_coulomb(x_1, x_2, alpha, a):
    return alpha / np.sqrt((x_1 - x_2) ** 2 + a ** 2)

def _compute_inner_integral(spf, l, num_grid_points, grid, alpha, a):
    inner_integral = np.zeros((l, l, num_grid_points), dtype=np.complex128)

    for i in range(num_grid_points):
        coulomb = _shielded_coulomb(grid[i], grid, alpha, a)
        for q in range(l):
            for s in range(l):
                inner_integral[q, s, i] = _trapz(
                    np.conjugate(spf[q]) * coulomb * spf[s],
                    grid,
                )

    return inner_integral


def _compute_orbital_integrals(spf, l, inner_integral, grid):
    u = np.zeros((l, l, l, l), dtype=np.complex128)
    for p in range(l):
        for q in range(l):
            for r in range(l):
                for s in range(l):
                    u[p, q, r, s] = _trapz(
                        np.conjugate(spf[p]) * inner_integral[q, s] * spf[r],
                        grid,
                    )

    return u

def add_spin_two_body(_u,):
    return np.kron(_u, np.einsum("pr, qs -> pqrs", np.eye(2), np.eye(2)))
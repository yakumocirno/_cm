def find_one_root_hill(c, tol=1e-10, max_iter=200000, restarts=60):
    R = cauchy_bound(c)

    best_z = None
    best_val = float("inf")

    for _ in range(restarts):
        # init
        r = R * (random.random() ** 0.5)
        theta = 2 * cmath.pi * random.random()
        z = r * cmath.exp(1j * theta)

        step = 0.2 * R  # initial neighborhood size

        for _it in range(max_iter):
            p, _ = poly_and_derivative(c, z)
            val = abs(p)
            if val < tol:
                return z

            if val < best_val:
                best_val = val
                best_z = z

            # try a few random neighbors, accept if improved
            improved = False
            for _ in range(20):
                dz = step * (random.random() - 0.5 + 1j * (random.random() - 0.5))
                z2 = z + dz
                p2, _ = poly_and_derivative(c, z2)
                if abs(p2) < val:
                    z = z2
                    improved = True
                    break

            # if stuck, shrink step; if too small, restart
            if not improved:
                step *= 0.7
                if step < 1e-10 * R:
                    break

    return best_z

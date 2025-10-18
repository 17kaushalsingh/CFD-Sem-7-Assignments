import argparse, numpy as np, time
from dataclasses import dataclass


# -----------------------------
# Inlined from solverLib/config
# -----------------------------
@dataclass
class Config:
    nx: int = 201
    ny: int = 51
    L: float = 10.0
    W: float = 1.0
    Re: float = 1000.0
    u0: float = 1.0


# ---------------------------
# Inlined from solverLib/utils
# ---------------------------
def create_grid(nx, ny, L, W):
    x = np.linspace(0.0, L, nx)
    y = np.linspace(0.0, W, ny)
    dx = x[1] - x[0]
    dy = y[1] - y[0]
    return x, y, dx, dy


# -----------------------------
# Inlined from solverLib/psisolve
# -----------------------------
def solve_psi_pointwise(psi, omega, dx, dy, tol=1e-8, maxiter=3000, omega_relax=1.0):
    nx, ny = psi.shape
    dx2, dy2 = dx * dx, dy * dy
    denom = 2.0 * (dx2 + dy2)
    for _ in range(maxiter):
        max_err = 0.0
        for i in range(1, nx - 1):
            for j in range(1, ny - 1):
                rhs = (
                    dx2 * (psi[i, j + 1] + psi[i, j - 1])
                    + dy2 * (psi[i + 1, j] + psi[i - 1, j])
                    + dx2 * dy2 * omega[i, j]
                )
                new_psi = rhs / denom
                diff = new_psi - psi[i, j]
                psi[i, j] += omega_relax * diff
                if abs(diff) > max_err:
                    max_err = abs(diff)
        if max_err < tol:
            break
    return psi


# ---------------------------------
# Inlined from solverLib/omegasolve
# ---------------------------------
def compute_velocity_from_psi(psi, dx, dy):
    nx, ny = psi.shape
    u, v = np.zeros_like(psi), np.zeros_like(psi)
    for i in range(1, nx - 1):
        for j in range(1, ny - 1):
            u[i, j] = (psi[i, j + 1] - psi[i, j - 1]) / (2 * dy)
            v[i, j] = -(psi[i + 1, j] - psi[i - 1, j]) / (2 * dx)
    return u, v


def advance_omega(omega, psi, dx, dy, dt, Re):
    nx, ny = omega.shape
    u, v = compute_velocity_from_psi(psi, dx, dy)
    new_omega = omega.copy()
    dx2, dy2 = dx * dx, dy * dy
    for i in range(1, nx - 1):
        for j in range(1, ny - 1):
            adv = (
                u[i, j] * (omega[i + 1, j] - omega[i - 1, j]) / (2 * dx)
                + v[i, j] * (omega[i, j + 1] - omega[i, j - 1]) / (2 * dy)
            )
            diff = (
                (omega[i + 1, j] - 2 * omega[i, j] + omega[i - 1, j]) / dx2
                + (omega[i, j + 1] - 2 * omega[i, j] + omega[i, j - 1]) / dy2
            )
            new_omega[i, j] = omega[i, j] + dt * (-adv + (1.0 / Re) * diff)
    return new_omega


# ---------------------------------
# Inlined from solverLib/boundary
# ---------------------------------
def set_initial_psi_from_inlet(psi, x, y, cfg):
    nx, ny = psi.shape
    W = cfg.W
    u0 = cfg.u0
    # parabolic velocity profile with mean u0
    u_inlet = 6.0 * u0 * (y / W) * (1.0 - y / W)
    psi_inlet = np.zeros_like(y)
    for j in range(1, ny):
        psi_inlet[j] = psi_inlet[j - 1] + 0.5 * (u_inlet[j] + u_inlet[j - 1]) * (y[j] - y[j - 1])
    for i in range(nx):
        psi[i, :] = psi_inlet
    psi[:, 0] = 0.0
    psi[:, -1] = psi_inlet[-1]


def apply_bcs_omega(omega, psi, cfg, dx, dy, dt=None, u_field=None):
    nx, ny = omega.shape

    # wall vorticity
    for i in range(1, nx - 1):
        omega[i, 0] = -2.0 * (psi[i, 1] - psi[i, 0]) / (dy * dy)
        omega[i, ny - 1] = -2.0 * (psi[i, ny - 2] - psi[i, ny - 1]) / (dy * dy)

    # inlet: set vorticity consistent with parabolic inlet
    W = cfg.W
    u0 = cfg.u0
    y = np.linspace(0, W, ny)
    u_in = 6.0 * u0 * (y / W) * (1.0 - y / W)
    du_dy = np.gradient(u_in, y)
    omega[0, :] = -du_dy

    # outlet: convective or Neumann
    if dt is not None and u_field is not None:
        for j in range(1, ny - 1):
            Ulocal = max(abs(u_field[-2, j]), 1e-8)
            omega[-1, j] = omega[-2, j] - (Ulocal * dt / dx) * (omega[-2, j] - omega[-3, j])
    else:
        omega[-1, :] = omega[-2, :]

    # corners
    omega[0, 0] = 0.5 * (omega[1, 0] + omega[0, 1])
    omega[0, ny - 1] = 0.5 * (omega[1, ny - 1] + omega[0, ny - 2])
    omega[nx - 1, 0] = 0.5 * (omega[nx - 2, 0] + omega[nx - 1, 1])
    omega[nx - 1, ny - 1] = 0.5 * (omega[nx - 2, ny - 1] + omega[nx - 1, ny - 2])


# ---------------------------------
# Inlined from solverLib/writer
# ---------------------------------
def write_tecplot(fname, x, y, psi, omega, u, v):
    nx, ny = x.size, y.size
    with open(fname, 'w') as f:
        f.write('TITLE = "Stream-Vorticity Solution"\n')
        f.write('VARIABLES = "X", "Y", "Psi", "Omega", "U", "V"\n')
        f.write(f'ZONE I={nx}, J={ny}, DATAPACKING=POINT\n')
        for j in range(ny):
            for i in range(nx):
                f.write(f"{x[i]:.6e} {y[j]:.6e} {psi[i,j]:.6e} {omega[i,j]:.6e} {u[i,j]:.6e} {v[i,j]:.6e}\n")


# -----------------
# CLI and main loop
# -----------------
def parse_args():
    p = argparse.ArgumentParser(
        description='2D streamfunction-vorticity CFD solver for a driven cavity / channel test case.'
    )
    p.add_argument('--nx', type=int, default=201, help='Number of grid points in x-direction (default: 201)')
    p.add_argument('--ny', type=int, default=51, help='Number of grid points in y-direction (default: 51)')
    p.add_argument('--L', type=float, default=10.0, help='Domain length in x (physical units, default: 10.0)')
    p.add_argument('--W', type=float, default=1.0, help='Domain height in y (physical units, default: 1.0)')
    p.add_argument('--Re', type=float, default=1000.0, help='Reynolds number (default: 1000)')
    p.add_argument('--u0', type=float, default=1.0, help='Inlet / lid velocity scale (default: 1.0)')
    p.add_argument('--itmax', type=int, default=2000, help='Maximum number of outer time-stepping iterations (default: 2000)')
    p.add_argument('--dt', type=float, default=None, help='Time step size; if omitted it is chosen from stability estimates')
    p.add_argument('--outfreq', type=int, default=10, help='Write output every N iterations (default: 10)')
    p.add_argument('--psi_tol', type=float, default=1e-8, help='Tolerance for iterative psi solver residual (default: 1e-8)')
    p.add_argument('--omega_tol', type=float, default=1e-6, help='Convergence tolerance for omega residual to stop simulation (default: 1e-6)')
    return p.parse_args()


def main():
    args = parse_args()
    cfg = Config(nx=args.nx, ny=args.ny, L=args.L, W=args.W, Re=args.Re, u0=args.u0)
    x, y, dx, dy = create_grid(cfg.nx, cfg.ny, cfg.L, cfg.W)
    psi, omega = np.zeros((cfg.nx, cfg.ny)), np.zeros((cfg.nx, cfg.ny))
    set_initial_psi_from_inlet(psi, x, y, cfg)

    if args.dt is None:
        conv_dt = min(dx, dy) / max(cfg.u0, 1e-8) * 0.5
        diff_dt = 0.25 * min(dx * dx, dy * dy) * cfg.Re
        dt = min(conv_dt, diff_dt)
    else:
        dt = args.dt

    apply_bcs_omega(omega, psi, cfg, dx, dy)
    outcount = 0

    logf = open('omegarsd.dat', 'w')
    logf.write('Iter Time[s] OmegaResidual\n')

    for it in range(1, args.itmax + 1):
        t0 = time.time()
        psi = solve_psi_pointwise(psi, omega, dx, dy, tol=args.psi_tol, maxiter=3000, omega_relax=1.0)
        new_omega = advance_omega(omega, psi, dx, dy, dt, cfg.Re)
        res = np.linalg.norm(new_omega - omega) / (np.linalg.norm(omega) + 1e-14)
        omega[:] = new_omega

        u, v = compute_velocity_from_psi(psi, dx, dy)
        apply_bcs_omega(omega, psi, cfg, dx, dy, dt=dt, u_field=u)
        t1 = time.time() - t0

        if it % args.outfreq == 0 or it == 1:
            outcount += 1
            fname = f'solution_{outcount:03d}.dat'
            write_tecplot(fname, x, y, psi, omega, u, v)

        print(f'V: Time {(dt * it):6.4f} | Iter {it:6d} | RunTime {t1:8.4f}s | OmegaRsd {res:10.4e}')
        logf.write(f'{it:6d} {(dt*it):6.4f} {res:12.6e}\n')

        if res < args.omega_tol:
            print(f'Omega solver converged at iteration {it}.')
            break

    logf.close()

    u, v = compute_velocity_from_psi(psi, dx, dy)
    write_tecplot('solution_final.dat', x, y, psi, omega, u, v)
    print('Simulation complete.')


if __name__ == '__main__':
    main()



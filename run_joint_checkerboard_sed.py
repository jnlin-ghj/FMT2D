#!/usr/bin/env python3
"""
FMT2D Joint Inversion Checkerboard Test — Pg + PsedP
=====================================================
Simultaneously invert velocity structure and sediment bottom interface
using Pg and PsedP phases.

Step 1: Generate initial velocity model (1D + flat sed bottom) + plot
Step 2: Generate true model (checkerboard velocity + undulating sed bottom) + plot
Step 3: Generate observation geometry (Pg + PsedP) + plot
Step 4: Forward modeling + plot ray paths and travel times
Step 5: Run 3 inversions and compare:
  - Pg only (velocity only)
  - Pg + PsedP (velocity only)
  - Pg + PsedP (velocity + sediment bottom interface)

Usage: python3 run_joint_checkerboard_sed.py [step]
  step: 1-5 or 'all' (default: all)
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import sys
import subprocess
from matplotlib.lines import Line2D

# ============================================================
# Configuration
# ============================================================
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
WORK_DIR = os.path.join(SCRIPT_DIR, 'example_inv_joint_sed')
FMT2D_BIN = os.path.join(SCRIPT_DIR, 'bin', 'fmt2d')
INVERT2D_BIN = os.path.join(SCRIPT_DIR, 'bin', 'invert2d')

# Model parameters
NX = 343
NZ = 35
VWATER = 1.5
VAIR = 1.5
X_MIN = 0.0
X_MAX = 114.0

# Seafloor
SEAFLOOR_BASE = 3.5
SEAFLOOR_AMP = 0.5
SEAFLOOR_WAVELENGTH = 30.0

# Velocity structure (1D reference model)
V_SEAFLOOR = 2.5
V_SED_BOTTOM = 4.5
V_UPPER_CRUST = 5.8
V_LOWER_CRUST = 6.8
V_MOHO = 7.0
V_BOTTOM = 8.1

# Depth levels (relative to seafloor)
Z_SED_BOTTOM = 2.0
Z_UPPER_CRUST = 4.5
Z_LOWER_CRUST = 7.0
Z_BOTTOM = 14.0

# Z positions (relative to seafloor)
ZPOS_REL = np.linspace(0, Z_BOTTOM, NZ)

# True Moho: flat (not inverted in this test)
MOHO_BASE = SEAFLOOR_BASE + Z_LOWER_CRUST  # = 10.5 km

# True sediment bottom: undulating
SED_BASE = SEAFLOOR_BASE + Z_SED_BOTTOM  # = 5.5 km
SED_AMP = 0.3
SED_WAVELENGTH = 25.0

# Checkerboard parameters
CB_DX = 18.0
CB_DZ = 5.0
CB_AMP = 0.05  # 5%

# Observation geometry
OBS_X = np.array([5, 12, 20, 28, 36, 44, 52, 60, 68, 76, 84, 92, 100, 110], dtype=float)
RCV_Z = 0.005
MIN_OFFSET_PG = 3.0
MAX_OFFSET_PG = 50.0
MIN_OFFSET_PSED = 3.0
MAX_OFFSET_PSED = 15.0
RCV_SPACING = 0.3

# Inversion parameters
INV_NITER = 20
INV_LAM_DAMP = 0.25
INV_LAM_SMOOTH = 0.4
INV_CRIT_CHI = 0.5
INV_LSQR_ITER = 1000
INV_TARGET_PERT = 0.015
INV_LAM_DAMP_IFACE = 0.1
INV_LAM_SMOOTH_IFACE = 1.5
INV_REFL_WEIGHT = 2.0
INV_PHASE_WT_PMP = 8.0
INV_PHASE_WT_PSED = 5.0
INV_CORR_X = 18.0
INV_CORR_Z = 5.0
INV_SMOOTH_TYPE = 1
INV_TARGET_DD = 0.2
INV_IFACE_STEP_RATIO = 1.0
INV_MAX_REL_PERT = 0.15
INV_ALT_MODE = 0

# Experiments
EXPERIMENTS = {
    'pg_vel': {
        'label': 'Pg only (vel)',
        'phases': [1],
        'inv_iface': False,
        'dir': 'pg_vel',
    },
    'pg_psed_vel': {
        'label': 'Pg+PsedP (vel)',
        'phases': [1, 2],
        'inv_iface': False,
        'dir': 'pg_psed_vel',
    },
    'pg_psed_vel_sed': {
        'label': 'Pg+PsedP (vel+sed)',
        'phases': [1, 2],
        'inv_iface': True,
        'dir': 'pg_psed_vel_sed',
    },
}

# Matplotlib global font settings for publication
plt.rcParams.update({
    'font.size': 18,
    'axes.labelsize': 18,
    'axes.titlesize': 18,
    'xtick.labelsize': 18,
    'ytick.labelsize': 18,
    'legend.fontsize': 16,
    'figure.titlesize': 22,
    'lines.linewidth': 2.0,
})


# ============================================================
# Utility functions
# ============================================================
def read_smesh(fname):
    with open(fname, 'r') as f:
        tokens = f.read().split()
    nx = int(tokens[0]); nz = int(tokens[1])
    vwater = float(tokens[2])
    x = np.array([float(tokens[4+i]) for i in range(nx)])
    topo = np.array([float(tokens[4+nx+i]) for i in range(nx)])
    z_rel = np.array([float(tokens[4+2*nx+i]) for i in range(nz)])
    vel = np.zeros((nz, nx))
    for ix in range(nx):
        for iz in range(nz):
            vel[iz, ix] = float(tokens[4+2*nx+nz+ix*nz+iz])
    return nx, nz, vwater, x, topo, z_rel, vel

def write_smesh(fname, nx, nz, vwater, vair, x, topo, z_rel, vel):
    with open(fname, 'w') as f:
        f.write(f' {nx:>8d} {nz:>8d} {vwater:>8.6f} {vair:>8.6f}\n')
        for ix in range(nx):
            f.write(f' {x[ix]:>13.6f}')
            if (ix + 1) % 8 == 0: f.write('\n')
        if nx % 8 != 0: f.write('\n')
        for ix in range(nx):
            f.write(f' {topo[ix]:>13.6f}')
            if (ix + 1) % 8 == 0: f.write('\n')
        if nx % 8 != 0: f.write('\n')
        for iz in range(nz):
            f.write(f' {z_rel[iz]:>13.6g}')
            if (iz + 1) % 8 == 0: f.write('\n')
        if nz % 8 != 0: f.write('\n')
        for ix in range(nx):
            for iz in range(nz):
                f.write(f' {vel[iz, ix]:>13.6f}')
                if (iz + 1) % 8 == 0: f.write('\n')
            f.write('\n')

def read_refl(fname):
    xs, zs = [], []
    with open(fname, 'r') as f:
        for line in f:
            parts = line.split()
            if len(parts) >= 2:
                xs.append(float(parts[0]))
                zs.append(float(parts[1]))
    return np.array(xs), np.array(zs)

def get_topo(x):
    return SEAFLOOR_BASE + SEAFLOOR_AMP * np.sin(2 * np.pi * x / SEAFLOOR_WAVELENGTH)

def build_vel_grid(x_smesh, topo, z_rel, vel, nz_grid=200):
    x_grid = x_smesh
    z_max = topo.max() + z_rel.max()
    z_grid = np.linspace(0, z_max, nz_grid)
    vel_grid = np.full((nz_grid, len(x_grid)), np.nan)
    for ix in range(len(x_grid)):
        for iz in range(nz_grid):
            z_abs = z_grid[iz]
            if z_abs < topo[ix]: continue
            z_r = z_abs - topo[ix]
            if z_r > z_rel[-1]: continue
            iz0 = np.searchsorted(z_rel, z_r) - 1
            iz0 = max(0, min(iz0, len(z_rel)-2))
            frac = (z_r - z_rel[iz0]) / (z_rel[iz0+1] - z_rel[iz0]) if z_rel[iz0+1] != z_rel[iz0] else 0
            vel_grid[iz, ix] = vel[iz0, ix] * (1-frac) + vel[iz0+1, ix] * frac
    return x_grid, z_grid, vel_grid

def plot_velocity_model(ax, x_smesh, topo, z_rel, vel, title, vmin=2.5, vmax=8.5,
                        plot_water=True, interfaces=None, sources=None):
    x_grid, z_grid, vel_grid = build_vel_grid(x_smesh, topo, z_rel, vel)
    im = ax.pcolormesh(x_grid, z_grid, vel_grid, cmap='jet_r', vmin=vmin, vmax=vmax, shading='auto')
    if plot_water:
        ax.fill_between(x_smesh, 0, topo, color='lightblue', alpha=0.3)
    ax.plot(x_smesh, topo, 'k-', linewidth=2)
    if interfaces:
        for iface_x, iface_z, style, color, label in interfaces:
            ax.plot(iface_x, iface_z, style, linewidth=2, label=label)   # removed color=color
    if sources is not None:
        for s in sources:
            ax.plot(s[0], s[1], 'r*', markersize=14, zorder=5)
    ax.set_ylim(topo.max() + z_rel.max(), 0)
    ax.set_xlabel('Distance (km)')
    ax.set_ylabel('Depth (km)')
    ax.set_title(title)
    return im

def run_cmd(cmd, cwd=None):
    print(f'  Running: {cmd}')
    result = subprocess.run(cmd, shell=True, cwd=cwd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f'  ERROR: {result.stderr[-500:]}')
    return result


# ============================================================
# Step 1: Generate initial velocity model (1D + flat sed bottom)
# ============================================================
def step1():
    print('\n' + '='*60)
    print('Step 1: Generate initial velocity model (1D + flat sed bottom)')
    print('='*60)

    xpos = np.linspace(X_MIN, X_MAX, NX)
    topo = get_topo(xpos)

    vel = np.zeros((NZ, NX))
    for ix in range(NX):
        for iz in range(NZ):
            z = ZPOS_REL[iz]
            if z < Z_SED_BOTTOM:
                vel[iz, ix] = V_SEAFLOOR + (V_SED_BOTTOM - V_SEAFLOOR) * z / Z_SED_BOTTOM
            elif z < Z_UPPER_CRUST:
                vel[iz, ix] = V_SED_BOTTOM + (V_UPPER_CRUST - V_SED_BOTTOM) * (z - Z_SED_BOTTOM) / (Z_UPPER_CRUST - Z_SED_BOTTOM)
            elif z < Z_LOWER_CRUST:
                vel[iz, ix] = V_UPPER_CRUST + (V_LOWER_CRUST - V_UPPER_CRUST) * (z - Z_UPPER_CRUST) / (Z_LOWER_CRUST - Z_UPPER_CRUST)
            else:
                vel[iz, ix] = V_MOHO + (V_BOTTOM - V_MOHO) * (z - Z_LOWER_CRUST) / (Z_BOTTOM - Z_LOWER_CRUST)

    smesh_file = os.path.join(WORK_DIR, 'init_vp.smesh')
    write_smesh(smesh_file, NX, NZ, VWATER, VAIR, xpos, topo, ZPOS_REL, vel)
    print(f'  Written: {smesh_file}')

    moho_flat = np.full(NX, MOHO_BASE)
    moho_file = os.path.join(WORK_DIR, 'init_moho.dat')
    with open(moho_file, 'w') as f:
        for ix in range(NX):
            f.write(f'{xpos[ix]:.6g} {moho_flat[ix]:.5f}\n')
    print(f'  Written: {moho_file} (flat Moho at {MOHO_BASE:.1f} km)')

    sed_flat = np.full(NX, SED_BASE)
    sed_file = os.path.join(WORK_DIR, 'init_sed_bottom.dat')
    with open(sed_file, 'w') as f:
        for ix in range(NX):
            f.write(f'{xpos[ix]:.6g} {sed_flat[ix]:.5f}\n')
    print(f'  Written: {sed_file} (flat sed bottom at {SED_BASE:.1f} km)')

    fig, ax = plt.subplots(figsize=(16, 6))
    im = plot_velocity_model(ax, xpos, topo, ZPOS_REL, vel, 'Initial Model (1D + Flat Sed Bottom)',
                             interfaces=[
                                 (xpos, moho_flat, 'r--', 'red', 'Moho (flat)'),
                                 (xpos, sed_flat, 'b--', 'blue', 'Sed bottom (flat)')
                             ])
    plt.colorbar(im, ax=ax, label='Vp (km/s)', shrink=0.8)
    ax.legend(loc='lower right')
    ax.set_xlim(X_MIN, X_MAX)
    plt.tight_layout()
    outfile = os.path.join(WORK_DIR, 'step1_init_model.pdf')
    plt.savefig(outfile, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'  Plot: {outfile}')


# ============================================================
# Step 2: Generate true model (checkerboard + undulating sed bottom)
# ============================================================
def step2():
    print('\n' + '='*60)
    print('Step 2: Generate true model (checkerboard + undulating sed bottom)')
    print('='*60)

    xpos = np.linspace(X_MIN, X_MAX, NX)
    topo = get_topo(xpos)

    # Moho is flat in this test (not inverted)
    moho_true = np.full(NX, MOHO_BASE)
    # Sediment bottom is undulating
    sed_true = SED_BASE + SED_AMP * np.sin(2 * np.pi * xpos / SED_WAVELENGTH + 1.0)

    vel = np.zeros((NZ, NX))
    for ix in range(NX):
        for iz in range(NZ):
            z = ZPOS_REL[iz]
            if z < Z_SED_BOTTOM:
                v_base = V_SEAFLOOR + (V_SED_BOTTOM - V_SEAFLOOR) * z / Z_SED_BOTTOM
            elif z < Z_UPPER_CRUST:
                v_base = V_SED_BOTTOM + (V_UPPER_CRUST - V_SED_BOTTOM) * (z - Z_SED_BOTTOM) / (Z_UPPER_CRUST - Z_SED_BOTTOM)
            elif z < Z_LOWER_CRUST:
                v_base = V_UPPER_CRUST + (V_LOWER_CRUST - V_UPPER_CRUST) * (z - Z_UPPER_CRUST) / (Z_LOWER_CRUST - Z_UPPER_CRUST)
            else:
                v_base = V_MOHO + (V_BOTTOM - V_MOHO) * (z - Z_LOWER_CRUST) / (Z_BOTTOM - Z_LOWER_CRUST)

            z_abs = topo[ix] + z
            if z > 0.01:
                cb_sign = np.sin(2 * np.pi * xpos[ix] / CB_DX) * np.sin(2 * np.pi * z_abs / CB_DZ)
                vel[iz, ix] = v_base * (1.0 + CB_AMP * cb_sign)
            else:
                vel[iz, ix] = v_base

    smesh_file = os.path.join(WORK_DIR, 'true_vp.smesh')
    write_smesh(smesh_file, NX, NZ, VWATER, VAIR, xpos, topo, ZPOS_REL, vel)
    print(f'  Written: {smesh_file}')

    moho_file = os.path.join(WORK_DIR, 'true_moho.dat')
    with open(moho_file, 'w') as f:
        for ix in range(NX):
            f.write(f'{xpos[ix]:.6g} {moho_true[ix]:.5f}\n')
    print(f'  Written: {moho_file} (Moho flat at {MOHO_BASE:.1f} km)')

    sed_file = os.path.join(WORK_DIR, 'true_sed_bottom.dat')
    with open(sed_file, 'w') as f:
        for ix in range(NX):
            f.write(f'{xpos[ix]:.6g} {sed_true[ix]:.5f}\n')
    print(f'  Written: {sed_file} (Sed bottom {SED_BASE:.1f} +/- {SED_AMP} km)')

    nx, nz, vwater, x_smesh, topo_r, z_rel, vel_init = read_smesh(os.path.join(WORK_DIR, 'init_vp.smesh'))
    solid = vel_init >= 2.0
    pert_pct = (vel - vel_init) / vel_init * 100
    print(f'  Velocity perturbation: {pert_pct[solid].min():.2f}% to {pert_pct[solid].max():.2f}%')

    fig, axes = plt.subplots(2, 1, figsize=(16, 10))

    ax = axes[0]
    im = plot_velocity_model(ax, xpos, topo, ZPOS_REL, vel, 'True Model (Checkerboard + Undulating Sed Bottom)',
                             interfaces=[
                                 (xpos, moho_true, 'r--', 'red', 'True Moho (flat)'),
                                 (xpos, sed_true, 'b--', 'blue', 'True Sed bottom')
                             ])
    plt.colorbar(im, ax=ax, label='Vp (km/s)', shrink=0.8)
    ax.legend(loc='lower right')
    ax.set_xlim(X_MIN, X_MAX)

    ax = axes[1]
    x_grid, z_grid, pert_grid = build_vel_grid(xpos, topo, ZPOS_REL, pert_pct)
    # 修改：使用 RdBu_r，与 step5_comparison.pdf 中反演结果色阶一致
    im = ax.pcolormesh(x_grid, z_grid, pert_grid, cmap='RdBu_r', vmin=-5, vmax=5, shading='auto')
    ax.fill_between(xpos, 0, topo, color='lightblue', alpha=0.3)
    ax.plot(xpos, topo, 'k-', linewidth=2)
    ax.plot(xpos, moho_true, 'r--', linewidth=2, label='True Moho')
    ax.plot(xpos, sed_true, 'b--', linewidth=2, label='True Sed bottom')
    plt.colorbar(im, ax=ax, label='Perturbation (%)', shrink=0.8)
    ax.set_ylim(topo.max() + ZPOS_REL[-1], 0)
    ax.set_xlim(X_MIN, X_MAX)
    ax.legend(loc='lower right')
    ax.set_title(f'Velocity Perturbation (checkerboard +/-{CB_AMP*100:.0f}%)')
    ax.set_xlabel('Distance (km)')
    ax.set_ylabel('Depth (km)')

    plt.tight_layout()
    outfile = os.path.join(WORK_DIR, 'step2_true_model.pdf')
    plt.savefig(outfile, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'  Plot: {outfile}')


# ============================================================
# Step 3: Generate observation geometry (Pg + PsedP)
# ============================================================
def step3():
    print('\n' + '='*60)
    print('Step 3: Generate observation geometry (Pg + PsedP)')
    print('='*60)

    all_src = []
    all_rcv = []

    for isrc, obs_x in enumerate(OBS_X):
        obs_z = get_topo(obs_x)  # source on seafloor

        pg_rcv = []
        psed_rcv = []

        # Pg: offset 3-50 km
        rcv_x_min = max(X_MIN + 0.5, obs_x - MAX_OFFSET_PG)
        rcv_x_max = min(X_MAX - 0.5, obs_x + MAX_OFFSET_PG)
        for rx in np.arange(rcv_x_min, rcv_x_max + 0.01, RCV_SPACING):
            offset = abs(rx - obs_x)
            if MIN_OFFSET_PG <= offset <= MAX_OFFSET_PG:
                pg_rcv.append((rx, RCV_Z, 1, 0.0, 0.0))  # ray_type=1 (Pg)

        # PsedP: offset 3-15 km
        rcv_x_min = max(X_MIN + 0.5, obs_x - MAX_OFFSET_PSED)
        rcv_x_max = min(X_MAX - 0.5, obs_x + MAX_OFFSET_PSED)
        for rx in np.arange(rcv_x_min, rcv_x_max + 0.01, RCV_SPACING):
            offset = abs(rx - obs_x)
            if MIN_OFFSET_PSED <= offset <= MAX_OFFSET_PSED:
                psed_rcv.append((rx, RCV_Z, 2, 0.0, 0.0))  # ray_type=2 (PsedP)

        rcv_list = pg_rcv + psed_rcv
        all_src.append((obs_x, obs_z, len(rcv_list)))
        all_rcv.append(rcv_list)

    geom_file = os.path.join(WORK_DIR, 'geom.dat')
    with open(geom_file, 'w') as f:
        f.write(f'{len(all_src)}\n')
        for isrc, (sx, sz, nrcv) in enumerate(all_src):
            f.write(f'S  {sx:.3f}  {sz:.5f}  {nrcv}\n')
            for rx, rz, rt, wt, we in all_rcv[isrc]:
                f.write(f'R  {rx:.3f}  {rz:.3f}  {rt}  {wt:.3f}  {we:.3f}\n')

    total_rcv = sum(len(r) for r in all_rcv)
    n_pg = sum(sum(1 for _, _, rt, _, _ in r if rt == 1) for r in all_rcv)
    n_psed = sum(sum(1 for _, _, rt, _, _ in r if rt == 2) for r in all_rcv)
    print(f'  Written: {geom_file}')
    print(f'  {len(all_src)} sources, {total_rcv} receivers (Pg: {n_pg}, PsedP: {n_psed})')

    fig, ax = plt.subplots(figsize=(16, 6))
    nx, nz, vwater, x_smesh, topo_arr, z_rel, vel = read_smesh(os.path.join(WORK_DIR, 'init_vp.smesh'))
    im = plot_velocity_model(ax, x_smesh, topo_arr, z_rel, vel, 'Observation Geometry (Pg + PsedP)')

    PHASE_COLORS_GEOM = {1: 'red', 2: 'green'}
    PHASE_NAMES_GEOM = {1: 'Pg', 2: 'PsedP'}

    for isrc, (sx, sz, nrcv) in enumerate(all_src):
        ax.plot(sx, sz, 'k*', markersize=14, zorder=5)
        for rx, rz, rt, _, _ in all_rcv[isrc]:
            ax.plot(rx, rz, 'v', color=PHASE_COLORS_GEOM[rt], markersize=2, zorder=4, alpha=0.5)

    legend_elements = [
        Line2D([0], [0], marker='v', color=PHASE_COLORS_GEOM[pid], linestyle='None',
               markersize=8, label=f'{PHASE_NAMES_GEOM[pid]} ({cnt})')
        for pid, cnt in [(1, n_pg), (2, n_psed)]
    ]
    ax.legend(handles=legend_elements, loc='lower right', fontsize=16)
    plt.colorbar(im, ax=ax, label='Vp (km/s)', shrink=0.8)
    ax.set_xlim(X_MIN, X_MAX)
    plt.tight_layout()
    outfile = os.path.join(WORK_DIR, 'step3_geometry.pdf')
    plt.savefig(outfile, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'  Plot: {outfile}')


# ============================================================
# Step 4: Forward modeling
# ============================================================
def step4():
    print('\n' + '='*60)
    print('Step 4: Forward modeling (true model)')
    print('='*60)

    fmt2d_in = os.path.join(WORK_DIR, 'fmt2d.in')
    with open(fmt2d_in, 'w') as f:
        f.write('true_vp.smesh\n')
        f.write('true_moho.dat\n')
        f.write('true_sed_bottom.dat\n')
        f.write('geom.dat\n')
        f.write('T\n')
        f.write('arrivals.dat\n')
        f.write('obs_times.dat\n')
        f.write('rays_src\n')
    print(f'  Written: {fmt2d_in}')

    result = run_cmd(f'{FMT2D_BIN} < fmt2d.in', cwd=WORK_DIR)
    print(f'  fmt2d exit code: {result.returncode}')

    arr_file = os.path.join(WORK_DIR, 'arrivals.dat')
    if os.path.exists(arr_file):
        arr = np.loadtxt(arr_file)
        print(f'  Arrivals: {len(arr)} rays')
        for p in [1, 2]:
            pmask = arr[:, 2] == p
            pname = {1: 'Pg', 2: 'PsedP'}.get(p, f'P{p}')
            if np.any(pmask):
                print(f'    {pname}: {np.sum(pmask)} rays, t_total range: {arr[pmask,3].min():.3f} - {arr[pmask,3].max():.3f} s')

    # Plot
    nx, nz, vwater, x_smesh, topo, z_rel, vel = read_smesh(os.path.join(WORK_DIR, 'true_vp.smesh'))

    sed_x, sed_z = read_refl(os.path.join(WORK_DIR, 'true_sed_bottom.dat'))
    sed_abs = np.interp(x_smesh, sed_x, sed_z)
    moho_x, moho_z = read_refl(os.path.join(WORK_DIR, 'true_moho.dat'))
    moho_abs = np.interp(x_smesh, moho_x, moho_z)

    sources = []
    receivers_per_src = []
    with open(os.path.join(WORK_DIR, 'geom.dat'), 'r') as f:
        nobs = int(f.readline().strip())
        for _ in range(nobs):
            parts = f.readline().split()
            sx, sz = float(parts[1]), float(parts[2])
            nrec = int(parts[3])
            sources.append((sx, sz))
            src_recs = []
            for _ in range(nrec):
                rparts = f.readline().split()
                rx, rz = float(rparts[1]), float(rparts[2])
                src_recs.append((rx, rz))
            receivers_per_src.append(src_recs)

    PHASE_NAMES = {1: 'Pg', 2: 'PsedP', 3: 'PmP'}
    PHASE_COLORS = {1: '#e41a1c', 2: '#4daf4a', 3: '#377eb8'}
    PHASE_MARKERS = {1: 'o', 2: '^', 3: 'D'}
    SRC_COLORS = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00',
                  '#a65628', '#f781bf', '#66c2a5', '#fc8d62', '#8da0cb']

    def read_rays_with_info(fname):
        rays, ray_info = [], []
        with open(fname, 'r') as f:
            while True:
                header = f.readline()
                if not header: break
                header = header.strip()
                if not header or not header.startswith('>'): continue
                parts = header[1:].split()
                if len(parts) < 4: continue
                try:
                    src_id = int(parts[0])
                    rcv_id = int(parts[1])
                    phase_id = int(parts[2])
                    travel_time = float(parts[3])
                except (ValueError, IndexError): continue
                ray_path = []
                while True:
                    pos = f.tell()
                    line = f.readline()
                    if not line: break
                    ls = line.strip()
                    if ls.startswith('>'):
                        f.seek(pos)
                        break
                    if not ls: continue
                    vals = ls.split()
                    if len(vals) >= 2:
                        try:
                            ray_path.append((float(vals[0]), float(vals[1])))
                        except ValueError:
                            f.seek(pos)
                            break
                if len(ray_path) > 0:
                    rays.append(np.array(ray_path))
                    ray_info.append((src_id, rcv_id, phase_id, travel_time))
        return rays, ray_info

    all_rays, all_info = [], []
    for src_file in sorted(os.listdir(WORK_DIR)):
        if src_file.startswith('rays_src_') and src_file.endswith('.dat'):
            rays, info = read_rays_with_info(os.path.join(WORK_DIR, src_file))
            all_rays.extend(rays)
            all_info.extend(info)

    phase_counts = {}
    for info in all_info:
        pid = info[2]
        phase_counts[pid] = phase_counts.get(pid, 0) + 1
    print(f'  Ray phases: {dict((PHASE_NAMES.get(k, f"P{k}"), v) for k, v in phase_counts.items())}')

    if len(phase_counts) == 0:
        print('  WARNING: No ray data found, skipping step4 plot')
        return

    # Build velocity mesh with water layer
    z_max_abs = topo.max() + z_rel.max()
    n_water = 10
    iz_max = nz
    n_total_z = n_water + iz_max

    vel_full = np.zeros((nx, n_total_z))
    z_mesh = np.zeros((nx, n_total_z))
    x_mesh = np.zeros((nx, n_total_z))
    for ix in range(nx):
        for iw in range(n_water):
            x_mesh[ix, iw] = x_smesh[ix]
            z_mesh[ix, iw] = topo[ix] * iw / (n_water - 1)
            vel_full[ix, iw] = vwater
        for iz in range(iz_max):
            idx = n_water + iz
            x_mesh[ix, idx] = x_smesh[ix]
            z_mesh[ix, idx] = topo[ix] + z_rel[iz]
            vel_full[ix, idx] = vel[iz, ix]

    tt_data = {}
    for idx_ray in range(len(all_info)):
        src_id, rcv_id, phase_id, travel_time = all_info[idx_ray]
        key = (src_id, phase_id)
        if key not in tt_data:
            tt_data[key] = ([], [])
        ray = all_rays[idx_ray]
        rx = ray[0, 0]
        tt_data[key][0].append(rx)
        tt_data[key][1].append(travel_time)

    phase_ids = sorted(phase_counts.keys())
    n_phases = len(phase_ids)

    fig, axes = plt.subplots(n_phases + 1, 1, figsize=(16, 4 * (n_phases + 1)))

    for i, pid in enumerate(phase_ids):
        ax = axes[i]
        pcm = ax.pcolormesh(x_mesh, z_mesh, vel_full, cmap='jet', vmin=1.5, vmax=8.5, shading='gouraud')
        ax.plot(x_smesh, topo, 'k-', linewidth=2)
        ax.plot(x_smesh, sed_abs, color='brown', linewidth=2, linestyle='--')
        ax.plot(x_smesh, moho_abs, 'k--', linewidth=2)

        for ray, info in zip(all_rays, all_info):
            if info[2] != pid: continue
            src_id = info[0]
            c = SRC_COLORS[(src_id - 1) % len(SRC_COLORS)]
            ax.plot(ray[:, 0], ray[:, 1], '-', color=c, linewidth=0.6, alpha=0.5)

        for j, (sx, sz) in enumerate(sources):
            ax.plot(sx, sz, 'k*', markersize=14, zorder=10)

        ax.set_ylabel('Depth (km)')
        pname = PHASE_NAMES.get(pid, f'P{pid}')
        ax.set_title(f'{pname} ({phase_counts[pid]} rays)')
        ax.invert_yaxis()
        ax.set_xlim(X_MIN, X_MAX)
        ax.set_ylim(z_max_abs, 0)
        plt.colorbar(pcm, ax=ax, label='Vp (km/s)', shrink=0.8)

    ax_tt = axes[n_phases]
    for pid in phase_ids:
        for src_id in range(1, len(sources) + 1):
            key = (src_id, pid)
            if key not in tt_data: continue
            rec_xs, times = tt_data[key]
            rec_xs = np.array(rec_xs)
            times = np.array(times)
            sort_idx = np.argsort(rec_xs)
            ax_tt.plot(rec_xs[sort_idx], times[sort_idx], PHASE_MARKERS.get(pid, 'o'),
                       color=PHASE_COLORS.get(pid, 'black'), markersize=4, alpha=0.6)

    legend_elements = [
        Line2D([0], [0], marker=PHASE_MARKERS.get(pid, 'o'), color=PHASE_COLORS.get(pid, 'black'),
               linestyle='None', markersize=8,
               label=f'{PHASE_NAMES.get(pid, f"P{pid}")} ({phase_counts.get(pid, 0)} rays)')
        for pid in phase_ids
    ]
    ax_tt.legend(handles=legend_elements, loc='upper left', fontsize=16, ncol=3)
    ax_tt.set_xlabel('Distance (km)')
    ax_tt.set_ylabel('Travel Time (s)')
    ax_tt.set_title('Travel Time Curves')
    ax_tt.grid(True, alpha=0.3)
    ax_tt.set_xlim(X_MIN, X_MAX)
    ax_tt.invert_yaxis()

    for j, (sx, sz) in enumerate(sources):
        ax_tt.axvline(x=sx, color=SRC_COLORS[j % len(SRC_COLORS)], linestyle=':', alpha=0.3, linewidth=1.5)

    fig.suptitle('FMT2D Forward Modeling - True Model (Pg + PsedP)', fontsize=22, fontweight='bold', y=0.99)
    fig.tight_layout(rect=[0, 0, 1, 0.98])

    # 修改：使走时子图宽度与前面的速度子图（经 colorbar 压缩后）一致
    pos_vel = axes[0].get_position()
    pos_tt = ax_tt.get_position()
    ax_tt.set_position([pos_vel.x0, pos_tt.y0, pos_vel.width, pos_tt.height])

    outfile = os.path.join(WORK_DIR, 'step4_forward.pdf')
    fig.savefig(outfile, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'  Plot: {outfile}')


# ============================================================
# Step 5: Run 3 inversions and compare
# ============================================================
def step5():
    print('\n' + '='*60)
    print('Step 5: Run 3 inversions and compare')
    print('='*60)

    # Read forward modeling results
    obs_file = os.path.join(WORK_DIR, 'obs_times.dat')
    if not os.path.exists(obs_file):
        print('  ERROR: obs_times.dat not found. Run step4 first.')
        return
    obs_data = np.loadtxt(obs_file, skiprows=1)

    # Read initial model info
    nx, nz_smesh, vwater, x_smesh, topo, z_rel, vel_init = read_smesh(os.path.join(WORK_DIR, 'init_vp.smesh'))
    _, _, _, _, _, _, vel_true = read_smesh(os.path.join(WORK_DIR, 'true_vp.smesh'))

    moho_init_x, moho_init_z = read_refl(os.path.join(WORK_DIR, 'init_moho.dat'))
    moho_true_x, moho_true_z = read_refl(os.path.join(WORK_DIR, 'true_moho.dat'))
    sed_init_x, sed_init_z = read_refl(os.path.join(WORK_DIR, 'init_sed_bottom.dat'))
    sed_true_x, sed_true_z = read_refl(os.path.join(WORK_DIR, 'true_sed_bottom.dat'))

    # Interpolate to smesh grid
    moho_init_z_interp = np.interp(x_smesh, moho_init_x, moho_init_z)
    moho_true_z_interp = np.interp(x_smesh, moho_true_x, moho_true_z)
    sed_init_z_interp = np.interp(x_smesh, sed_init_x, sed_init_z)
    sed_true_z_interp = np.interp(x_smesh, sed_true_x, sed_true_z)

    results = {}
    exp_list = list(EXPERIMENTS.keys())

    for exp_name in exp_list:
        exp = EXPERIMENTS[exp_name]
        exp_dir = os.path.join(WORK_DIR, exp['dir'])
        os.makedirs(exp_dir, exist_ok=True)

        print(f'\n--- Running inversion: {exp["label"]} ---')

        # Clean old files
        for f in os.listdir(exp_dir):
            fp = os.path.join(exp_dir, f)
            if os.path.isfile(fp):
                os.remove(fp)

        # Copy input files
        import shutil
        for fname in ['init_vp.smesh', 'init_moho.dat', 'init_sed_bottom.dat', 'geom.dat']:
            shutil.copy2(os.path.join(WORK_DIR, fname), os.path.join(exp_dir, fname))

        # Prepare geom_inv.dat with travel times for this experiment's phases
        obs_lookup = {}
        for row in obs_data:
            key = (int(row[0]), int(row[1]), int(row[2]))
            obs_lookup[key] = row[3]

        inv_phases = set(exp['phases'])
        default_unc = 0.05
        geom_orig = os.path.join(exp_dir, 'geom.dat')
        geom_inv = os.path.join(exp_dir, 'geom_inv.dat')
        n_valid = 0
        with open(geom_orig, 'r') as fin, open(geom_inv, 'w') as fout:
            nobs = int(fin.readline().strip())
            fout.write(f'{nobs}\n')
            ircv_global = 0
            for iobs in range(1, nobs + 1):
                line = fin.readline()
                parts = line.split()
                sx, sz, nrec = float(parts[1]), float(parts[2]), int(parts[3])
                fout.write(f'S  {sx:.3f}  {sz:.5f}  {nrec}\n')
                for _ in range(nrec):
                    rline = fin.readline()
                    rparts = rline.split()
                    rx, rz, rt = float(rparts[1]), float(rparts[2]), int(rparts[3])
                    ircv_global += 1
                    key = (ircv_global, iobs, rt)
                    if rt in inv_phases and key in obs_lookup:
                        t_obs = obs_lookup[key]
                        fout.write(f'R  {rx:.3f}  {rz:.3f}  {rt}  {t_obs:.6f}  {default_unc:.6f}\n')
                        n_valid += 1
                    else:
                        fout.write(f'R  {rx:.3f}  {rz:.3f}  {rt}  -1.000000  {default_unc:.6f}\n')
        print(f'  geom_inv.dat: {n_valid} valid observations')

        # Write frechet.in
        frechet_file = os.path.join(exp_dir, 'frechet.in')
        with open(frechet_file, 'w') as f:
            f.write('1\n')     # do_frechet = 1
            f.write('1\n')     # n_inv_vgrid = 1
            f.write('1\n')     # vgrid index = 1
            f.write('1\n')     # vtype = 1 (Vp)
            if exp['inv_iface']:
                f.write('1\n')  # n_inv_iface = 1
                f.write('3\n')  # iface index = 3 (sediment bottom)
            else:
                f.write('0\n')  # n_inv_iface = 0
        print(f'  Written: {frechet_file} (inv_iface={exp["inv_iface"]})')

        # Write invert2d.in
        inv_in = os.path.join(exp_dir, 'invert2d.in')
        with open(inv_in, 'w') as f:
            f.write('init_vp.smesh\n')
            f.write('init_moho.dat\n')
            f.write('init_sed_bottom.dat\n')
            f.write('geom_inv.dat\n')
            f.write(f'{INV_NITER}\n')
            f.write(f'{INV_LAM_DAMP}\n')
            f.write(f'{INV_LAM_SMOOTH}\n')
            f.write(f'{INV_CRIT_CHI}\n')
            f.write(f'{INV_LSQR_ITER}\n')
            f.write(f'{INV_TARGET_PERT}\n')
            f.write(f'{INV_LAM_DAMP_IFACE}\n')
            f.write(f'{INV_LAM_SMOOTH_IFACE}\n')
            f.write(f'{INV_REFL_WEIGHT}\n')
            f.write(f'{INV_PHASE_WT_PMP}\n')
            f.write(f'1.0\n')  # phase_wt_pgss (unused)
            f.write(f'{INV_PHASE_WT_PSED}\n')
            f.write(f'{INV_CORR_X}\n')
            f.write(f'{INV_CORR_Z}\n')
            f.write(f'{INV_SMOOTH_TYPE} {INV_TARGET_DD} {INV_IFACE_STEP_RATIO} {INV_MAX_REL_PERT}\n')
            f.write(f'{INV_ALT_MODE}\n')
        print(f'  Written: {inv_in}')

        # Run invert2d
        result = run_cmd(f'{INVERT2D_BIN}', cwd=exp_dir)
        print(f'  invert2d exit code: {result.returncode}')

        log_file = os.path.join(exp_dir, 'invert2d.log')
        with open(log_file, 'w') as f:
            f.write(result.stdout)

        results[exp_name] = exp_dir

    # ============================================================
    # Read and compare results
    # ============================================================
    corr_all = {}
    corr_shallow = {}
    rms_diff = {}
    sed_inv_z_dict = {}

    for exp_name, exp_dir in results.items():
        model_files = sorted([f for f in os.listdir(exp_dir) if f.startswith('model_iter') and f.endswith('.dat')],
                             key=lambda f: int(f.replace('model_iter','').replace('.dat','')))
        if not model_files:
            print(f'  WARNING: No model_iter files found for {exp_name}!')
            continue

        best_model = model_files[-1]
        print(f'  {EXPERIMENTS[exp_name]["label"]}: best model = {best_model}')

        # Read inverted model
        def read_model_iter_joint(fname, nx):
            with open(fname, 'r') as f:
                lines = f.readlines()
            vel_tokens = []
            iface_data = {}
            section = None
            current_iface = None
            current_tokens = []
            for line in lines:
                stripped = line.strip()
                if stripped.startswith('Iteration'):
                    continue
                elif stripped.startswith('vgrid'):
                    section = 'vgrid'
                    continue
                elif stripped.startswith('iface'):
                    if section == 'iface' and current_iface and current_tokens:
                        iface_data[current_iface] = np.array([float(t) for t in current_tokens])
                    section = 'iface'
                    current_iface = stripped
                    current_tokens = []
                    continue
                if section == 'vgrid':
                    vel_tokens.extend(stripped.split())
                elif section == 'iface':
                    current_tokens.extend(stripped.split())
            if section == 'iface' and current_iface and current_tokens:
                iface_data[current_iface] = np.array([float(t) for t in current_tokens])

            nz_inv = len(vel_tokens) // nx
            vel = np.zeros((nz_inv, nx))
            for ix in range(nx):
                for iz in range(nz_inv):
                    vel[iz, ix] = float(vel_tokens[ix*nz_inv + iz])

            sed_z = iface_data.get('iface 3', None)
            return vel, sed_z, nz_inv

        vel_inv, sed_inv_z, nz_inv = read_model_iter_joint(
            os.path.join(exp_dir, best_model), nx)
        z_inv = np.linspace(0, 23.0, nz_inv)

        if sed_inv_z is None or len(sed_inv_z) == 0:
            sed_inv_file = os.path.join(exp_dir, 'iface_3.refl')
            if os.path.exists(sed_inv_file):
                _, sed_inv_z = read_refl(sed_inv_file)
            else:
                sed_inv_z = None

        sed_inv_z_dict[exp_name] = sed_inv_z

        # Build initial model on inversion grid
        vel_init_inv = np.zeros((nz_inv, nx))
        for ix in range(nx):
            for iz in range(nz_inv):
                z_abs = z_inv[iz]
                if z_abs < topo[ix]:
                    vel_init_inv[iz, ix] = vwater
                else:
                    z_r = z_abs - topo[ix]
                    if z_r > z_rel[-1]:
                        vel_init_inv[iz, ix] = vel_init[-1, ix]
                    else:
                        iz0 = np.searchsorted(z_rel, z_r) - 1
                        iz0 = max(0, min(iz0, len(z_rel)-2))
                        frac = (z_r - z_rel[iz0]) / (z_rel[iz0+1] - z_rel[iz0]) if z_rel[iz0+1] != z_rel[iz0] else 0
                        vel_init_inv[iz, ix] = vel_init[iz0, ix] * (1-frac) + vel_init[iz0+1, ix] * frac

        # Build true model on inversion grid
        vel_true_inv = np.zeros((nz_inv, nx))
        for ix in range(nx):
            for iz in range(nz_inv):
                z_abs = z_inv[iz]
                if z_abs < topo[ix]:
                    vel_true_inv[iz, ix] = vwater
                else:
                    z_r = z_abs - topo[ix]
                    if z_r > z_rel[-1]:
                        vel_true_inv[iz, ix] = vel_true[-1, ix]
                    else:
                        iz0 = np.searchsorted(z_rel, z_r) - 1
                        iz0 = max(0, min(iz0, len(z_rel)-2))
                        frac = (z_r - z_rel[iz0]) / (z_rel[iz0+1] - z_rel[iz0]) if z_rel[iz0+1] != z_rel[iz0] else 0
                        vel_true_inv[iz, ix] = vel_true[iz0, ix] * (1-frac) + vel_true[iz0+1, ix] * frac

        # Compute perturbations
        pert_true = np.full((nz_inv, nx), np.nan)
        pert_inv = np.full((nz_inv, nx), np.nan)
        for ix in range(nx):
            for iz in range(nz_inv):
                if vel_init_inv[iz, ix] >= 2.0:
                    pert_true[iz, ix] = (vel_true_inv[iz, ix] - vel_init_inv[iz, ix]) / vel_init_inv[iz, ix] * 100
                    pert_inv[iz, ix] = (vel_inv[iz, ix] - vel_init_inv[iz, ix]) / vel_init_inv[iz, ix] * 100

        # Correlation
        solid = vel_init_inv >= 2.0
        valid = ~np.isnan(pert_true) & ~np.isnan(pert_inv)
        mask = solid & valid
        corr_all[exp_name] = np.corrcoef(pert_true[mask], pert_inv[mask])[0, 1] if mask.any() else 0.0

        # Shallow correlation (above sed bottom + 1 km)
        shallow_mask = mask.copy()
        for ix in range(nx):
            for iz in range(nz_inv):
                if z_inv[iz] > SED_BASE + 1.0:
                    shallow_mask[iz, ix] = False
        corr_shallow[exp_name] = np.corrcoef(pert_true[shallow_mask], pert_inv[shallow_mask])[0, 1] if shallow_mask.any() else 0.0

        # RMS Diff
        pert_diff = pert_inv - pert_true
        diff_mask = ~np.isnan(pert_diff) & solid
        rms_diff[exp_name] = np.sqrt(np.mean(pert_diff[diff_mask]**2)) if diff_mask.any() else 0.0

        print(f'  {EXPERIMENTS[exp_name]["label"]}: Corr={corr_all[exp_name]:.4f}, '
              f'Shallow Corr={corr_shallow[exp_name]:.4f}, RMS Diff={rms_diff[exp_name]:.3f}%')

    # ============================================================
    # Figure 1: Perturbation recovery comparison
    # ============================================================
    n_exp = len(exp_list)
    n_rows = n_exp + 1  # true + experiments
    fig, axes = plt.subplots(n_rows, 1, figsize=(16, 4 * n_rows))

    # True perturbation
    ax = axes[0]
    # Re-read for plotting
    vel_init_inv_plot = np.zeros((nz_inv, nx))
    vel_true_inv_plot = np.zeros((nz_inv, nx))
    for ix in range(nx):
        for iz in range(nz_inv):
            z_abs = z_inv[iz]
            if z_abs < topo[ix]:
                vel_init_inv_plot[iz, ix] = vwater
                vel_true_inv_plot[iz, ix] = vwater
            else:
                z_r = z_abs - topo[ix]
                if z_r > z_rel[-1]:
                    vel_init_inv_plot[iz, ix] = vel_init[-1, ix]
                    vel_true_inv_plot[iz, ix] = vel_true[-1, ix]
                else:
                    iz0 = np.searchsorted(z_rel, z_r) - 1
                    iz0 = max(0, min(iz0, len(z_rel)-2))
                    frac = (z_r - z_rel[iz0]) / (z_rel[iz0+1] - z_rel[iz0]) if z_rel[iz0+1] != z_rel[iz0] else 0
                    vel_init_inv_plot[iz, ix] = vel_init[iz0, ix] * (1-frac) + vel_init[iz0+1, ix] * frac
                    vel_true_inv_plot[iz, ix] = vel_true[iz0, ix] * (1-frac) + vel_true[iz0+1, ix] * frac

    pert_true_plot = np.full((nz_inv, nx), np.nan)
    for ix in range(nx):
        for iz in range(nz_inv):
            if vel_init_inv_plot[iz, ix] >= 2.0:
                pert_true_plot[iz, ix] = (vel_true_inv_plot[iz, ix] - vel_init_inv_plot[iz, ix]) / vel_init_inv_plot[iz, ix] * 100

    x_mesh_inv = np.zeros((nx, nz_inv))
    z_mesh_inv = np.zeros((nx, nz_inv))
    for ix in range(nx):
        for iz in range(nz_inv):
            x_mesh_inv[ix, iz] = x_smesh[ix]
            z_mesh_inv[ix, iz] = z_inv[iz]

    im = ax.pcolormesh(x_mesh_inv, z_mesh_inv, pert_true_plot.T, cmap='RdBu_r', vmin=-5, vmax=5, shading='auto')
    ax.fill_between(x_smesh, 0, topo, color='lightblue', alpha=0.3)
    ax.plot(x_smesh, topo, 'k-', linewidth=2)
    ax.plot(x_smesh, moho_true_z_interp, 'r--', linewidth=2, label='True Moho')
    ax.plot(x_smesh, sed_true_z_interp, 'b--', linewidth=2, label='True Sed')
    plt.colorbar(im, ax=ax, label='Perturbation (%)', shrink=0.8)
    ax.set_ylim(z_inv[-1], 0)
    ax.set_xlim(X_MIN, X_MAX)
    ax.set_title('True Perturbation (Checkerboard)')
    ax.set_ylabel('Depth (km)')
    ax.legend(fontsize=16)

    # Inverted perturbations
    for idx, exp_name in enumerate(exp_list):
        ax = axes[idx + 1]
        exp_dir = results[exp_name]
        model_files = sorted([f for f in os.listdir(exp_dir) if f.startswith('model_iter') and f.endswith('.dat')],
                             key=lambda f: int(f.replace('model_iter','').replace('.dat','')))
        if not model_files:
            continue
        vel_inv, sed_inv_z_local, nz_inv_local = read_model_iter_joint(
            os.path.join(exp_dir, model_files[-1]), nx)

        pert_inv_plot = np.full((nz_inv, nx), np.nan)
        for ix in range(nx):
            for iz in range(nz_inv):
                if vel_init_inv_plot[iz, ix] >= 2.0:
                    pert_inv_plot[iz, ix] = (vel_inv[iz, ix] - vel_init_inv_plot[iz, ix]) / vel_init_inv_plot[iz, ix] * 100

        im = ax.pcolormesh(x_mesh_inv, z_mesh_inv, pert_inv_plot.T, cmap='RdBu_r', vmin=-5, vmax=5, shading='auto')
        ax.fill_between(x_smesh, 0, topo, color='lightblue', alpha=0.3)
        ax.plot(x_smesh, topo, 'k-', linewidth=2)
        ax.plot(x_smesh, moho_true_z_interp, 'r--', linewidth=2, label='True Moho')
        ax.plot(x_smesh, sed_true_z_interp, 'b--', linewidth=2, label='True Sed')
        if sed_inv_z_dict.get(exp_name) is not None:
            ax.plot(x_smesh, sed_inv_z_dict[exp_name], 'g-', linewidth=2, label='Inv Sed')
        plt.colorbar(im, ax=ax, label='Perturbation (%)', shrink=0.8)
        ax.set_ylim(z_inv[-1], 0)
        ax.set_xlim(X_MIN, X_MAX)
        label = EXPERIMENTS[exp_name]['label']
        ax.set_title(f'{label} — Corr={corr_all.get(exp_name,0):.3f}, Shallow={corr_shallow.get(exp_name,0):.3f}')
        ax.set_ylabel('Depth (km)')
        if idx == len(exp_list) - 1:
            ax.set_xlabel('Distance (km)')
        ax.legend(fontsize=16, loc='lower right')

    fig.suptitle('Velocity Perturbation Recovery (Pg + PsedP)', fontsize=22, fontweight='bold')
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    outfile = os.path.join(WORK_DIR, 'step5_comparison.pdf')
    fig.savefig(outfile, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'  Plot: {outfile}')

    # ============================================================
    # Figure 2: Sediment bottom interface comparison
    # ============================================================
    fig2, ax2 = plt.subplots(figsize=(14, 5))

    ax2.plot(x_smesh, sed_true_z_interp, 'r-', linewidth=2, label='True Sed Bottom')
    ax2.plot(x_smesh, sed_init_z_interp, 'b--', linewidth=2, label='Initial Sed Bottom (flat)')
    sed_colors = {'pg_vel': 'gray', 'pg_psed_vel': 'orange', 'pg_psed_vel_sed': 'green'}
    for exp_name in exp_list:
        sed_z = sed_inv_z_dict.get(exp_name)
        if sed_z is not None:
            ax2.plot(x_smesh, sed_z, '-', color=sed_colors.get(exp_name, 'black'),
                     linewidth=2, label=f'{EXPERIMENTS[exp_name]["label"]}')

    sed_err_init = np.sqrt(np.mean((sed_init_z_interp - sed_true_z_interp)**2))
    ax2.set_title(f'Sediment Bottom Recovery (Init RMSE={sed_err_init:.3f}km)', fontsize=20, fontweight='bold')
    ax2.set_xlabel('X (km)')
    ax2.set_ylabel('Depth (km)')
    ax2.invert_yaxis()
    ax2.legend(fontsize=16)
    ax2.grid(True, alpha=0.3)
    fig2.tight_layout()
    outfile2 = os.path.join(WORK_DIR, 'step5_sed_recovery.pdf')
    fig2.savefig(outfile2, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'  Plot: {outfile2}')

    # ============================================================
    # Figure 3: Bar comparison
    # ============================================================
    fig3, axes3 = plt.subplots(1, 3, figsize=(15, 5))

    labels = [EXPERIMENTS[e]['label'] for e in exp_list]
    x_pos = np.arange(len(exp_list))
    bar_colors = ['#e41a1c', '#ff7f00', '#4daf4a']

    # Correlation
    vals = [corr_all.get(e, 0) for e in exp_list]
    axes3[0].bar(x_pos, vals, color=bar_colors[:len(exp_list)], edgecolor='black', linewidth=0.5)
    axes3[0].set_xticks(x_pos)
    axes3[0].set_xticklabels(labels, fontsize=14, rotation=15, ha='right')
    axes3[0].set_ylabel('Correlation', fontsize=18)
    axes3[0].set_title('Overall Correlation', fontsize=20)
    axes3[0].set_ylim(0, 1)
    for i, v in enumerate(vals):
        axes3[0].text(i, v + 0.02, f'{v:.3f}', ha='center', fontsize=16)

    # Shallow correlation
    vals = [corr_shallow.get(e, 0) for e in exp_list]
    axes3[1].bar(x_pos, vals, color=bar_colors[:len(exp_list)], edgecolor='black', linewidth=0.5)
    axes3[1].set_xticks(x_pos)
    axes3[1].set_xticklabels(labels, fontsize=14, rotation=15, ha='right')
    axes3[1].set_ylabel('Correlation', fontsize=18)
    axes3[1].set_title('Shallow Correlation', fontsize=20)
    axes3[1].set_ylim(0, 1)
    for i, v in enumerate(vals):
        axes3[1].text(i, v + 0.02, f'{v:.3f}', ha='center', fontsize=16)

    # RMS Diff
    vals = [rms_diff.get(e, 0) for e in exp_list]
    axes3[2].bar(x_pos, vals, color=bar_colors[:len(exp_list)], edgecolor='black', linewidth=0.5)
    axes3[2].set_xticks(x_pos)
    axes3[2].set_xticklabels(labels, fontsize=14, rotation=15, ha='right')
    axes3[2].set_ylabel('RMS Diff (%)', fontsize=18)
    axes3[2].set_title('RMS Difference', fontsize=20)
    for i, v in enumerate(vals):
        axes3[2].text(i, v + 0.01, f'{v:.3f}', ha='center', fontsize=16)

    fig3.suptitle('Inversion Quality Comparison', fontsize=22, fontweight='bold')
    fig3.tight_layout(rect=[0, 0, 1, 0.95])
    outfile3 = os.path.join(WORK_DIR, 'step5_bar_comparison.pdf')
    fig3.savefig(outfile3, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'  Plot: {outfile3}')

    # ============================================================
    # Figure 4: Convergence
    # ============================================================
    fig4, ax4 = plt.subplots(figsize=(10, 6))
    conv_colors = {'pg_vel': '#e41a1c', 'pg_psed_vel': '#ff7f00', 'pg_psed_vel_sed': '#4daf4a'}

    for exp_name in exp_list:
        exp_dir = results[exp_name]
        log_file = os.path.join(exp_dir, 'invert_log_latest.txt')
        if not os.path.exists(log_file):
            continue
        iters, rms_vals = [], []
        with open(log_file, 'r') as f:
            for line in f:
                if not line.strip().startswith('Iteration'): continue
                try:
                    idx_rms = line.index('RMS=')
                    idx_rms_new = line.index('RMS_new=')
                    after_rms = line[idx_rms+4:idx_rms_new].strip()
                    rms = float(after_rms.split()[0])
                    it = int(line.split()[1])
                    iters.append(it)
                    rms_vals.append(rms)
                except (ValueError, IndexError):
                    continue
        if iters:
            ax4.plot(iters, rms_vals, 'o-', color=conv_colors.get(exp_name, 'black'),
                     markersize=8, linewidth=2, label=EXPERIMENTS[exp_name]['label'])

    ax4.set_xlabel('Iteration', fontsize=20)
    ax4.set_ylabel('RMS Travel Time Residual (s)', fontsize=20)
    ax4.set_title('Inversion Convergence', fontsize=22, fontweight='bold')
    ax4.legend(fontsize=16)
    ax4.grid(True, alpha=0.3)
    ax4.set_yscale('log')
    fig4.tight_layout()
    outfile4 = os.path.join(WORK_DIR, 'step5_convergence.pdf')
    fig4.savefig(outfile4, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'  Plot: {outfile4}')

    # ============================================================
    # Figure 5: Ray density comparison
    # ============================================================
    n_cols = min(3, n_exp)
    n_rows_rd = (n_exp + n_cols - 1) // n_cols
    fig5, axes5 = plt.subplots(n_rows_rd, n_cols, figsize=(7 * n_cols, 5 * n_rows_rd))
    if n_exp == 1:
        axes5 = np.array([[axes5]])
    elif n_rows_rd == 1:
        axes5 = axes5[np.newaxis, :]

    for idx, (exp_name, exp_dir) in enumerate(results.items()):
        row = idx // n_cols
        col = idx % n_cols
        ax = axes5[row, col]
        ray_dens_file = os.path.join(exp_dir, 'ray_density_final.dat')
        if not os.path.exists(ray_dens_file):
            ax.text(0.5, 0.5, 'No ray density data', transform=ax.transAxes, ha='center')
            continue

        with open(ray_dens_file, 'r') as f:
            tokens = f.read().split()
        rd_nx = int(tokens[0])
        rd_nz = int(tokens[1])
        rd_x = np.array([float(tokens[2 + i]) for i in range(rd_nx)])
        rd_z = np.array([float(tokens[2 + rd_nx + i]) for i in range(rd_nz)])
        ray_density = np.zeros((rd_nz, rd_nx))
        tidx = 2 + rd_nx + rd_nz
        while tidx + 2 < len(tokens):
            ix = int(tokens[tidx])
            iz = int(tokens[tidx + 1])
            cnt = float(tokens[tidx + 2])
            if 1 <= ix <= rd_nx and 1 <= iz <= rd_nz:
                ray_density[iz - 1, ix - 1] = cnt
            tidx += 3

        # Mask water
        for ix_rd in range(rd_nx):
            x_val = rd_x[ix_rd]
            topo_val = np.interp(x_val, x_smesh, topo)
            for iz_rd in range(rd_nz):
                if rd_z[iz_rd] < topo_val:
                    ray_density[iz_rd, ix_rd] = 0.0

        rd_plot = ray_density.copy()
        rd_plot[rd_plot <= 0] = np.nan
        rd_max = np.nanmax(rd_plot) if np.any(~np.isnan(rd_plot)) else 1.0
        im = ax.pcolormesh(rd_x, rd_z, rd_plot, cmap='magma', vmin=0, vmax=rd_max, shading='nearest')
        ax.fill_between(x_smesh, 0, topo, color='lightblue', alpha=0.6, zorder=2)
        ax.plot(x_smesh, topo, 'k-', linewidth=2, zorder=3)
        ax.plot(x_smesh, moho_true_z_interp, 'r--', linewidth=2, label='True Moho', zorder=3)
        ax.plot(x_smesh, sed_true_z_interp, 'b--', linewidth=2, label='True Sed', zorder=3)
        plt.colorbar(im, ax=ax, label='Ray Density', shrink=0.8)
        ax.set_ylim(max(moho_true_z_interp) + 3.0, 0)
        ax.set_xlim(X_MIN, X_MAX)
        ax.set_title(f'{EXPERIMENTS[exp_name]["label"]} — Ray Density', fontsize=18)
        ax.set_xlabel('Distance (km)', fontsize=16)
        if col == 0:
            ax.set_ylabel('Depth (km)', fontsize=16)
        ax.legend(fontsize=14, loc='lower right')

    # Hide unused subplots
    for idx in range(n_exp, n_rows_rd * n_cols):
        row = idx // n_cols
        col = idx % n_cols
        axes5[row, col].set_visible(False)

    fig5.suptitle('Ray Density Comparison', fontsize=22, fontweight='bold')
    fig5.tight_layout(rect=[0, 0, 1, 0.95])
    outfile5 = os.path.join(WORK_DIR, 'step5_ray_density.pdf')
    fig5.savefig(outfile5, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'  Plot: {outfile5}')

    # ============================================================
    # Summary
    # ============================================================
    print('\n' + '='*60)
    print('SUMMARY: Pg + PsedP Joint Inversion Assessment')
    print('='*60)
    print(f'{"Experiment":<25} {"Corr":>8} {"Shallow":>8} {"RMS%":>8} {"Sed RMSE":>10}')
    print('-'*59)
    for exp_name in exp_list:
        label = EXPERIMENTS[exp_name]['label']
        sed_z = sed_inv_z_dict.get(exp_name)
        if sed_z is not None:
            sed_rmse = np.sqrt(np.mean((sed_z - sed_true_z_interp)**2))
            sed_rmse_str = f'{sed_rmse:.3f}km'
        else:
            sed_rmse_str = 'N/A'
        print(f'{label:<25} {corr_all.get(exp_name,0):>8.4f} {corr_shallow.get(exp_name,0):>8.4f} {rms_diff.get(exp_name,0):>8.3f} {sed_rmse_str:>10}')
    print()


# ============================================================
# Main
# ============================================================
def cleanup_workspace():
    """Remove old output files before a fresh run."""
    patterns = [
        'model_iter*.dat', 'iface_*.refl', 'residuals_*.dat',
        'ray_density_*.dat', 'arrivals.dat', 'obs_times.dat',
        'obs_times_inv.dat', 'geom_inv.dat', 'invert_log_*.txt', 'invert2d.log',
        'vgrid_*.smesh', 'rays_src_*.dat', 'res_*.dat',
        'entry_points.dat',
    ]
    import glob as _glob
    for pat in patterns:
        for f in _glob.glob(os.path.join(WORK_DIR, pat)):
            try:
                os.remove(f)
            except OSError:
                pass


if __name__ == '__main__':
    step = sys.argv[1] if len(sys.argv) > 1 else 'all'

    os.makedirs(WORK_DIR, exist_ok=True)

    if step == 'all':
        print('Cleaning up old output files...')
        cleanup_workspace()

    if step == 'all' or step == '1':
        step1()
    if step == 'all' or step == '2':
        step2()
    if step == 'all' or step == '3':
        step3()
    if step == 'all' or step == '4':
        step4()
    if step == 'all' or step == '5':
        step5()

    print('\nDone!')

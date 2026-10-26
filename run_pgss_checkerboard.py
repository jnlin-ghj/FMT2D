#!/usr/bin/env python3
"""
FMT2D Pg_ss Phase Value Assessment — Checkerboard Test
========================================================
Evaluate the inversion value of the Pg_ss phase (sea-surface + seafloor
double multiple) by comparing velocity-only inversions with different
phase combinations:

  A) Pg only
  B) Pg + Pg_ss
  C) Pg + Pg_ss2
  D) Pg + Pg_ss + Pg_ss2

Steps:
  1: Generate initial velocity model (1D + true interfaces)
  2: Generate true model (checkerboard + undulating interfaces)
  3: Generate observation geometry (Pg + Pg_ss + Pg_ss2)
  4: Forward modeling with true model
  5: Run 4 inversions and compare

Usage: python run_pgss_checkerboard.py [step]
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
WORK_DIR = os.path.join(SCRIPT_DIR, 'example_inv_pgss')
FMT2D_BIN = os.path.join(SCRIPT_DIR, 'bin', 'fmt2d')
INVERT2D_BIN = os.path.join(SCRIPT_DIR, 'bin', 'invert2d')

# Model parameters (same as run_vel_checkerboard.py)
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

ZPOS_REL = np.linspace(0, Z_BOTTOM, NZ)

# True Moho
MOHO_BASE = SEAFLOOR_BASE + Z_LOWER_CRUST
MOHO_AMP = 0.72
MOHO_WAVELENGTH = 40.0

# True sediment bottom
SED_BASE = SEAFLOOR_BASE + Z_SED_BOTTOM
SED_AMP = 0.3
SED_WAVELENGTH = 25.0

# Checkerboard parameters
CB_DX = 18.0
CB_DZ = 5.0
CB_AMP = 0.05  # 5%

# Observation geometry — more OBS to get good Pg_ss coverage
OBS_X = np.array([5, 13, 21, 29, 37, 45, 53, 61, 69, 77, 85, 93, 101, 109], dtype=float)
RCV_Z = 0.005
MIN_OFFSET_PG = 3.0
MAX_OFFSET_PG = 50.0
MIN_OFFSET_PGSS = 3.0
MAX_OFFSET_PGSS = 50.0
RCV_SPACING = 0.3

# Inversion parameters
INV_NITER = 15
INV_LAM_DAMP = 0.15
INV_LAM_SMOOTH = 0.3
INV_CRIT_CHI = 0.5
INV_LSQR_ITER = 1000
INV_TARGET_PERT = 0.02
INV_LAM_DAMP_IFACE = 0.0
INV_LAM_SMOOTH_IFACE = 0.0
INV_REFL_WEIGHT = 2.0
INV_PHASE_WT_PMP = 8.0
INV_PHASE_WT_PGSS = 1.0
INV_PHASE_WT_PSED = 5.0
INV_CORR_X = 18.0
INV_CORR_Z = 5.0
INV_SMOOTH_TYPE = 1
INV_TARGET_DD = 0.2
INV_IFACE_STEP_RATIO = 1.0
INV_MAX_REL_PERT = 0.20
INV_ALT_MODE = 0

# Inversion experiments — compare Pg only with multi-phase combinations
# (Pg_ss only and Pg_ss2 only are omitted as they are rarely used alone)
EXPERIMENTS = {
    'pg_only':   {'phases': [1],       'label': 'Pg only'},
    'pg_pgss':   {'phases': [1, 4],    'label': 'Pg + Pg_ss'},
    'pg_pgss2':  {'phases': [1, 5],    'label': 'Pg + Pg_ss2'},
    'all':       {'phases': [1, 4, 5], 'label': 'Pg + Pg_ss + Pg_ss2'},
}

PHASE_NAMES = {1: 'Pg', 2: 'PsedP', 3: 'PmP', 4: 'Pg_ss', 5: 'Pg_ss2'}
PHASE_COLORS = {1: '#e41a1c', 2: '#4daf4a', 3: '#377eb8', 4: '#ff7f00', 5: '#984ea3'}
PHASE_MARKERS = {1: 'o', 2: '^', 3: 'D', 4: 's', 5: 'p'}

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
            ax.plot(iface_x, iface_z, style, linewidth=2, label=label)
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

def read_rays_with_info(fname):
    rays, ray_info = [], []
    with open(fname, 'r') as f:
        while True:
            header = f.readline()
            if not header:
                break
            header = header.strip()
            if not header or not header.startswith('>'):
                continue
            parts = header[1:].split()
            if len(parts) < 4:
                continue
            try:
                src_id = int(parts[0])
                rcv_id = int(parts[1])
                phase_id = int(parts[2])
                travel_time = float(parts[3])
            except (ValueError, IndexError):
                continue
            ray_path = []
            while True:
                pos = f.tell()
                line = f.readline()
                if not line:
                    break
                ls = line.strip()
                if ls.startswith('>'):
                    f.seek(pos)
                    break
                if not ls:
                    continue
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


# ============================================================
# Step 1: Generate initial velocity model
# ============================================================
def step1():
    print('\n' + '='*60)
    print('Step 1: Generate initial velocity model (1D + True Interfaces)')
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

    moho_true = MOHO_BASE + MOHO_AMP * np.sin(2 * np.pi * xpos / MOHO_WAVELENGTH)
    moho_file = os.path.join(WORK_DIR, 'init_moho.dat')
    with open(moho_file, 'w') as f:
        for ix in range(NX):
            f.write(f'{xpos[ix]:.6g} {moho_true[ix]:.5f}\n')
    print(f'  Written: {moho_file}')

    sed_true = SED_BASE + SED_AMP * np.sin(2 * np.pi * xpos / SED_WAVELENGTH + 1.0)
    sed_file = os.path.join(WORK_DIR, 'init_sed_bottom.dat')
    with open(sed_file, 'w') as f:
        for ix in range(NX):
            f.write(f'{xpos[ix]:.6g} {sed_true[ix]:.5f}\n')
    print(f'  Written: {sed_file}')

    fig, ax = plt.subplots(figsize=(16, 6))
    im = plot_velocity_model(ax, xpos, topo, ZPOS_REL, vel, 'Initial Model (1D + True Interfaces)',
                             interfaces=[
                                 (xpos, moho_true, 'r--', 'red', 'Moho (true)'),
                                 (xpos, sed_true, 'b--', 'blue', 'Sed bottom (true)')
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
# Step 2: Generate true model (checkerboard)
# ============================================================
def step2():
    print('\n' + '='*60)
    print('Step 2: Generate true model (checkerboard + undulating Moho)')
    print('='*60)

    xpos = np.linspace(X_MIN, X_MAX, NX)
    topo = get_topo(xpos)

    moho_true = MOHO_BASE + MOHO_AMP * np.sin(2 * np.pi * xpos / MOHO_WAVELENGTH)
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

    sed_file = os.path.join(WORK_DIR, 'true_sed_bottom.dat')
    with open(sed_file, 'w') as f:
        for ix in range(NX):
            f.write(f'{xpos[ix]:.6g} {sed_true[ix]:.5f}\n')

    nx, nz, vwater, x_smesh, topo_r, z_rel, vel_init = read_smesh(os.path.join(WORK_DIR, 'init_vp.smesh'))
    solid = vel_init >= 2.0
    pert_pct = (vel - vel_init) / vel_init * 100
    print(f'  Velocity perturbation: {pert_pct[solid].min():.2f}% to {pert_pct[solid].max():.2f}%')

    fig, axes = plt.subplots(2, 1, figsize=(16, 10))

    ax = axes[0]
    im = plot_velocity_model(ax, xpos, topo, ZPOS_REL, vel, 'True Model (Checkerboard + Undulating Moho)',
                             interfaces=[
                                 (xpos, moho_true, 'r--', 'red', 'True Moho'),
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
# Step 3: Generate observation geometry (Pg + Pg_ss + Pg_ss2)
# ============================================================
def step3():
    print('\n' + '='*60)
    print('Step 3: Generate observation geometry (Pg + Pg_ss + Pg_ss2)')
    print('='*60)

    all_src = []
    all_rcv = []

    for isrc, obs_x in enumerate(OBS_X):
        obs_z = get_topo(obs_x)

        pg_rcv = []
        pgss_rcv = []
        pgss2_rcv = []

        # Pg: offset 3-50 km
        rcv_x_min = max(X_MIN + 0.5, obs_x - MAX_OFFSET_PG)
        rcv_x_max = min(X_MAX - 0.5, obs_x + MAX_OFFSET_PG)
        for rx in np.arange(rcv_x_min, rcv_x_max + 0.01, RCV_SPACING):
            offset = abs(rx - obs_x)
            if MIN_OFFSET_PG <= offset <= MAX_OFFSET_PG:
                pg_rcv.append((rx, RCV_Z, 1, 0.0, 0.0))  # ray_type=1 (Pg)

        # Pg_ss: offset 3-50 km
        rcv_x_min = max(X_MIN + 0.5, obs_x - MAX_OFFSET_PGSS)
        rcv_x_max = min(X_MAX - 0.5, obs_x + MAX_OFFSET_PGSS)
        for rx in np.arange(rcv_x_min, rcv_x_max + 0.01, RCV_SPACING):
            offset = abs(rx - obs_x)
            if MIN_OFFSET_PGSS <= offset <= MAX_OFFSET_PGSS:
                pgss_rcv.append((rx, RCV_Z, 4, 0.0, 0.0))  # ray_type=4 (Pg_ss)

        # Pg_ss2: offset 3-50 km (same range as Pg_ss)
        for rx in np.arange(rcv_x_min, rcv_x_max + 0.01, RCV_SPACING):
            offset = abs(rx - obs_x)
            if MIN_OFFSET_PGSS <= offset <= MAX_OFFSET_PGSS:
                pgss2_rcv.append((rx, RCV_Z, 5, 0.0, 0.0))  # ray_type=5 (Pg_ss2)

        rcv_list = pg_rcv + pgss_rcv + pgss2_rcv
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
    n_pgss = sum(sum(1 for _, _, rt, _, _ in r if rt == 4) for r in all_rcv)
    n_pgss2 = sum(sum(1 for _, _, rt, _, _ in r if rt == 5) for r in all_rcv)
    print(f'  Written: {geom_file}')
    print(f'  {len(all_src)} sources, {total_rcv} receivers (Pg: {n_pg}, Pg_ss: {n_pgss}, Pg_ss2: {n_pgss2})')

    fig, ax = plt.subplots(figsize=(16, 6))
    nx, nz, vwater, x_smesh, topo_arr, z_rel, vel = read_smesh(os.path.join(WORK_DIR, 'init_vp.smesh'))
    im = plot_velocity_model(ax, x_smesh, topo_arr, z_rel, vel, 'Observation Geometry (Pg + Pg_ss + Pg_ss2)')

    PHASE_COLORS_GEOM = {1: 'red', 4: 'orange', 5: 'purple'}
    PHASE_NAMES_GEOM = {1: 'Pg', 4: 'Pg_ss', 5: 'Pg_ss2'}

    for isrc, (sx, sz, nrcv) in enumerate(all_src):
        ax.plot(sx, sz, 'k*', markersize=14, zorder=5)
        for rx, rz, rt, _, _ in all_rcv[isrc]:
            ax.plot(rx, rz, 'v', color=PHASE_COLORS_GEOM[rt], markersize=2, zorder=4, alpha=0.5)

    legend_elements = [
        Line2D([0], [0], marker='v', color=PHASE_COLORS_GEOM[pid], linestyle='None',
               markersize=8, label=f'{PHASE_NAMES_GEOM[pid]} ({cnt})')
        for pid, cnt in [(1, n_pg), (4, n_pgss), (5, n_pgss2)]
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
        for p in [1, 4, 5]:
            pmask = arr[:, 2] == p
            pname = PHASE_NAMES.get(p, f'P{p}')
            if np.any(pmask):
                print(f'    {pname}: {np.sum(pmask)} rays, t_total range: {arr[pmask,3].min():.3f} - {arr[pmask,3].max():.3f} s')

    # Plot ray paths
    nx, nz, vwater, x_smesh, topo, z_rel, vel = read_smesh(os.path.join(WORK_DIR, 'true_vp.smesh'))
    sed_x, sed_z = read_refl(os.path.join(WORK_DIR, 'true_sed_bottom.dat'))
    sed_abs = np.interp(x_smesh, sed_x, sed_z)
    moho_x, moho_z = read_refl(os.path.join(WORK_DIR, 'true_moho.dat'))
    moho_abs = np.interp(x_smesh, moho_x, moho_z)

    sources = []
    with open(os.path.join(WORK_DIR, 'geom.dat'), 'r') as f:
        nobs = int(f.readline().strip())
        for _ in range(nobs):
            parts = f.readline().split()
            sx, sz = float(parts[1]), float(parts[2])
            nrec = int(parts[3])
            sources.append((sx, sz))
            for _ in range(nrec):
                f.readline()

    SRC_COLORS = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00',
                  '#a65628', '#f781bf', '#66c2a5', '#fc8d62', '#8da0cb']

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
                       color=PHASE_COLORS.get(pid, 'gray'), markersize=4, alpha=0.6)

    legend_elements = [
        Line2D([0], [0], marker=PHASE_MARKERS.get(pid, 'o'), color=PHASE_COLORS.get(pid, 'gray'),
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

    fig.suptitle('FMT2D Forward Modeling - True Model (Pg + Pg_ss)', fontsize=22, fontweight='bold', y=0.99)
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
# Step 5: Run inversions and compare
# ============================================================
def run_single_inversion(exp_name, phase_ids_inv):
    """Run a single velocity-only inversion with specified phases."""
    exp_dir = os.path.join(WORK_DIR, exp_name)
    os.makedirs(exp_dir, exist_ok=True)

    # Clean old output
    import glob as _glob
    for pat in ['model_iter*.dat', 'iface_*.refl', 'residuals_*.dat',
                'ray_density_*.dat', 'invert_log_*.txt', 'invert2d.log',
                'vgrid_*.smesh', 'rays_src_*.dat', 'res_*.dat', 'entry_points.dat']:
        for f in _glob.glob(os.path.join(exp_dir, pat)):
            try:
                os.remove(f)
            except OSError:
                pass

    # Prepare geom_inv.dat with travel times for specified phases
    obs_file = os.path.join(WORK_DIR, 'obs_times.dat')
    obs_data = np.loadtxt(obs_file, skiprows=1)

    obs_lookup = {}
    for row in obs_data:
        key = (int(row[0]), int(row[1]), int(row[2]))
        obs_lookup[key] = row[3]

    inv_phases = set(phase_ids_inv)
    default_unc = 0.05

    # Symlink shared files into experiment directory
    shared_files = ['init_vp.smesh', 'init_moho.dat', 'init_sed_bottom.dat', 'geom.dat']
    for fname in shared_files:
        src = os.path.join(WORK_DIR, fname)
        dst = os.path.join(exp_dir, fname)
        if os.path.islink(dst):
            os.remove(dst)
        elif os.path.exists(dst):
            os.remove(dst)
        os.symlink(src, dst)

    # Create geom_inv.dat from geom.dat + obs_times.dat
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
        f.write('0\n')     # n_inv_iface = 0

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
        f.write(f'{INV_PHASE_WT_PGSS}\n')
        f.write(f'{INV_PHASE_WT_PSED}\n')
        f.write(f'{INV_CORR_X}\n')
        f.write(f'{INV_CORR_Z}\n')
        f.write(f'{INV_SMOOTH_TYPE} {INV_TARGET_DD} {INV_IFACE_STEP_RATIO} {INV_MAX_REL_PERT}\n')
        f.write(f'{INV_ALT_MODE}\n')

    # Run invert2d
    result = run_cmd(f'{INVERT2D_BIN}', cwd=exp_dir)
    print(f'  invert2d exit code: {result.returncode}')

    log_file = os.path.join(exp_dir, 'invert2d.log')
    with open(log_file, 'w') as f:
        f.write(result.stdout)

    return exp_dir


def read_model_iter(fname, nx):
    """Read model_iter file and return velocity perturbation grid."""
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
    return vel, nz_inv


def compute_perturbation(vel, vel_init, topo, z_rel, z_inv, nx, nz_inv, vel_on_inv_grid=False):
    """Compute perturbation on inversion grid."""
    def interp_to_inv_grid(v, topo, z_rel, z_inv, nx, nz_inv):
        v_inv = np.zeros((nz_inv, nx))
        for ix in range(nx):
            for iz in range(nz_inv):
                z_abs = z_inv[iz]
                if z_abs < topo[ix]:
                    v_inv[iz, ix] = 1.5
                else:
                    z_r = z_abs - topo[ix]
                    if z_r > z_rel[-1]:
                        v_inv[iz, ix] = v[-1, ix]
                    else:
                        iz0 = np.searchsorted(z_rel, z_r) - 1
                        iz0 = max(0, min(iz0, len(z_rel)-2))
                        frac = (z_r - z_rel[iz0]) / (z_rel[iz0+1] - z_rel[iz0]) if z_rel[iz0+1] != z_rel[iz0] else 0
                        v_inv[iz, ix] = v[iz0, ix] * (1-frac) + v[iz0+1, ix] * frac
        return v_inv

    if vel_on_inv_grid:
        vel_inv = vel
    else:
        vel_inv = interp_to_inv_grid(vel, topo, z_rel, z_inv, nx, nz_inv)
    vel_init_inv = interp_to_inv_grid(vel_init, topo, z_rel, z_inv, nx, nz_inv)

    pert = np.full((nz_inv, nx), np.nan)
    for ix in range(nx):
        for iz in range(nz_inv):
            if vel_init_inv[iz, ix] >= 2.0:
                pert[iz, ix] = (vel_inv[iz, ix] - vel_init_inv[iz, ix]) / vel_init_inv[iz, ix] * 100
    return pert, vel_init_inv


def step5():
    print('\n' + '='*60)
    print('Step 5: Run 4 inversions and compare')
    print('='*60)

    results = {}
    for exp_name, exp_cfg in EXPERIMENTS.items():
        print(f'\n--- Running inversion: {exp_cfg["label"]} ---')
        exp_dir = run_single_inversion(exp_name, exp_cfg['phases'])
        results[exp_name] = exp_dir

    # ============================================================
    # Read all results
    # ============================================================
    nx, nz_smesh, vwater, x_smesh, topo, z_rel, vel_init = read_smesh(
        os.path.join(WORK_DIR, 'init_vp.smesh'))
    _, _, _, _, _, _, vel_true = read_smesh(os.path.join(WORK_DIR, 'true_vp.smesh'))

    moho_true_x, moho_true_z = read_refl(os.path.join(WORK_DIR, 'true_moho.dat'))
    sed_true_x, sed_true_z = read_refl(os.path.join(WORK_DIR, 'true_sed_bottom.dat'))

    first_exp = list(results.values())[0]
    model_files = sorted([f for f in os.listdir(first_exp) if f.startswith('model_iter') and f.endswith('.dat')],
                         key=lambda f: int(f.replace('model_iter','').replace('.dat','')))
    if not model_files:
        print('  ERROR: No model_iter files found!')
        return

    vel_sample, nz_inv = read_model_iter(os.path.join(first_exp, model_files[-1]), nx)
    z_inv = np.linspace(0, 23.0, nz_inv)

    pert_true, vel_init_inv = compute_perturbation(vel_true, vel_init, topo, z_rel, z_inv, nx, nz_inv)

    pert_inv = {}
    corr_all = {}
    corr_shallow = {}
    rms_diff = {}

    for exp_name, exp_dir in results.items():
        model_files = sorted([f for f in os.listdir(exp_dir) if f.startswith('model_iter') and f.endswith('.dat')],
                             key=lambda f: int(f.replace('model_iter','').replace('.dat','')))
        if not model_files:
            print(f'  WARNING: No model_iter files for {exp_name}')
            continue

        vel_inv, nz_inv_exp = read_model_iter(os.path.join(exp_dir, model_files[-1]), nx)
        pert, _ = compute_perturbation(vel_inv, vel_init, topo, z_rel, z_inv, nx, nz_inv, vel_on_inv_grid=True)
        pert_inv[exp_name] = pert

        solid = vel_init_inv >= 2.0
        valid = ~np.isnan(pert_true) & ~np.isnan(pert)
        mask = solid & valid
        if mask.any():
            corr_all[exp_name] = np.corrcoef(pert_true[mask], pert[mask])[0, 1]
        else:
            corr_all[exp_name] = 0.0

        shallow_mask = mask.copy()
        for ix in range(nx):
            for iz in range(nz_inv):
                if z_inv[iz] > MOHO_BASE:
                    shallow_mask[iz, ix] = False
        if shallow_mask.any():
            corr_shallow[exp_name] = np.corrcoef(pert_true[shallow_mask], pert[shallow_mask])[0, 1]
        else:
            corr_shallow[exp_name] = 0.0

        pert_diff = pert - pert_true
        diff_mask = ~np.isnan(pert_diff) & solid
        if diff_mask.any():
            rms_diff[exp_name] = np.sqrt(np.mean(pert_diff[diff_mask]**2))
        else:
            rms_diff[exp_name] = 0.0

        print(f'  {EXPERIMENTS[exp_name]["label"]}: Corr={corr_all[exp_name]:.4f}, '
              f'Shallow Corr={corr_shallow[exp_name]:.4f}, RMS Diff={rms_diff[exp_name]:.3f}%')

    # ============================================================
    # Figure 1: Comparison of perturbation recovery (vertical layout)
    # ============================================================
    x_mesh_inv = np.zeros((nx, nz_inv))
    z_mesh_inv = np.zeros((nx, nz_inv))
    for ix in range(nx):
        for iz in range(nz_inv):
            x_mesh_inv[ix, iz] = x_smesh[ix]
            z_mesh_inv[ix, iz] = z_inv[iz]

    exp_list = list(EXPERIMENTS.keys())
    n_exp = len(exp_list)
    n_panels = n_exp + 1  # +1 for true perturbation

    fig, axes = plt.subplots(n_panels, 1, figsize=(16, 4.5 * n_panels))

    # True perturbation
    ax = axes[0]
    im = ax.pcolormesh(x_mesh_inv, z_mesh_inv, pert_true.T, cmap='RdBu_r', vmin=-5, vmax=5, shading='auto')
    ax.fill_between(x_smesh, 0, topo, color='lightblue', alpha=0.3)
    ax.plot(x_smesh, topo, 'k-', linewidth=2)
    ax.plot(x_smesh, moho_true_z, 'r--', linewidth=2, label='True Moho')
    ax.plot(x_smesh, sed_true_z, 'b--', linewidth=2, label='True Sed')
    plt.colorbar(im, ax=ax, label='Perturbation (%)', shrink=0.8)
    ax.set_ylim(z_inv[-1], 0)
    ax.set_xlim(X_MIN, X_MAX)
    ax.set_title('True Perturbation (Checkerboard +/-5%)')
    ax.set_ylabel('Depth (km)')
    ax.legend(fontsize=16)

    # Each experiment
    for idx, exp_name in enumerate(exp_list):
        ax = axes[idx + 1]
        if exp_name not in pert_inv:
            ax.text(0.5, 0.5, 'No data', transform=ax.transAxes, ha='center')
            continue
        im = ax.pcolormesh(x_mesh_inv, z_mesh_inv, pert_inv[exp_name].T, cmap='RdBu_r', vmin=-5, vmax=5, shading='auto')
        ax.fill_between(x_smesh, 0, topo, color='lightblue', alpha=0.3)
        ax.plot(x_smesh, topo, 'k-', linewidth=2)
        ax.plot(x_smesh, moho_true_z, 'r--', linewidth=2)
        ax.plot(x_smesh, sed_true_z, 'b--', linewidth=2)
        plt.colorbar(im, ax=ax, label='Perturbation (%)', shrink=0.8)
        ax.set_ylim(z_inv[-1], 0)
        ax.set_xlim(X_MIN, X_MAX)
        label = EXPERIMENTS[exp_name]['label']
        ax.set_title(f'{label} — Corr={corr_all.get(exp_name,0):.3f}, Shallow={corr_shallow.get(exp_name,0):.3f}, RMS={rms_diff.get(exp_name,0):.2f}%')
        ax.set_ylabel('Depth (km)')
        if idx == n_exp - 1:
            ax.set_xlabel('Distance (km)')

    fig.suptitle('Pg_ss Phase Value Assessment — Velocity Perturbation Recovery', fontsize=22, fontweight='bold')
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    outfile = os.path.join(WORK_DIR, 'step5_comparison.pdf')
    fig.savefig(outfile, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'  Plot: {outfile}')

    # ============================================================
    # Figure 2: Perturbation difference maps (vertical layout)
    # ============================================================
    fig2, axes2 = plt.subplots(n_exp, 1, figsize=(16, 4.5 * n_exp))

    for idx, exp_name in enumerate(exp_list):
        ax = axes2[idx]
        if exp_name not in pert_inv:
            ax.text(0.5, 0.5, 'No data', transform=ax.transAxes, ha='center')
            continue
        pert_diff = pert_inv[exp_name] - pert_true
        vmax_diff = max(np.nanmax(np.abs(pert_diff[~np.isnan(pert_diff)])), 1.0) if not np.all(np.isnan(pert_diff)) else 2.0
        im = ax.pcolormesh(x_mesh_inv, z_mesh_inv, pert_diff.T, cmap='RdBu_r',
                           vmin=-vmax_diff, vmax=vmax_diff, shading='auto')
        ax.fill_between(x_smesh, 0, topo, color='lightblue', alpha=0.3)
        ax.plot(x_smesh, topo, 'k-', linewidth=2)
        ax.plot(x_smesh, moho_true_z, 'r--', linewidth=2)
        ax.plot(x_smesh, sed_true_z, 'b--', linewidth=2)
        plt.colorbar(im, ax=ax, label='Diff (%)', shrink=0.8)
        ax.set_ylim(z_inv[-1], 0)
        ax.set_xlim(X_MIN, X_MAX)
        label = EXPERIMENTS[exp_name]['label']
        ax.set_title(f'{label} — RMS Diff = {rms_diff.get(exp_name,0):.3f}%')
        ax.set_ylabel('Depth (km)')
        if idx == n_exp - 1:
            ax.set_xlabel('Distance (km)')

    fig2.suptitle('Perturbation Difference (Recovered - True)', fontsize=22, fontweight='bold')
    fig2.tight_layout(rect=[0, 0, 1, 0.97])
    outfile2 = os.path.join(WORK_DIR, 'step5_difference.pdf')
    fig2.savefig(outfile2, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'  Plot: {outfile2}')

    # ============================================================
    # Figure 3: Convergence curves
    # ============================================================
    fig3, ax3 = plt.subplots(figsize=(10, 6))
    exp_colors = {'pg_only': '#e41a1c', 'pg_pgss': '#ff7f00',
                  'pg_pgss2': '#984ea3', 'all': '#4daf4a'}

    for exp_name, exp_dir in results.items():
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
            ax3.plot(iters, rms_vals, 'o-', color=exp_colors.get(exp_name, 'gray'),
                     markersize=8, linewidth=2, label=EXPERIMENTS[exp_name]['label'])

    ax3.set_xlabel('Iteration', fontsize=20)
    ax3.set_ylabel('RMS Travel Time Residual (s)', fontsize=20)
    ax3.set_title('Inversion Convergence Comparison', fontsize=22, fontweight='bold')
    ax3.legend(fontsize=16)
    ax3.grid(True, alpha=0.3)
    ax3.set_yscale('log')
    fig3.tight_layout()
    outfile3 = os.path.join(WORK_DIR, 'step5_convergence.pdf')
    fig3.savefig(outfile3, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'  Plot: {outfile3}')

    # ============================================================
    # Summary
    # ============================================================
    print('\n' + '='*60)
    print('SUMMARY: Pg_ss Phase Value Assessment')
    print('='*60)
    print(f'{"Experiment":<25} {"Corr":>8} {"Shallow":>8} {"RMS%":>8}')
    print('-'*50)
    for exp_name in exp_list:
        label = EXPERIMENTS[exp_name]['label']
        print(f'{label:<25} {corr_all.get(exp_name,0):>8.4f} {corr_shallow.get(exp_name,0):>8.4f} {rms_diff.get(exp_name,0):>8.3f}')
    print()


# ============================================================
# Main
# ============================================================
def cleanup_workspace():
    """Remove old output files before a fresh run."""
    import glob as _glob
    patterns = [
        'model_iter*.dat', 'iface_*.refl', 'residuals_*.dat',
        'ray_density_*.dat', 'arrivals.dat', 'obs_times.dat',
        'obs_times_inv.dat', 'geom_inv.dat', 'invert_log_*.txt', 'invert2d.log',
        'vgrid_*.smesh', 'rays_src_*.dat', 'res_*.dat',
        'entry_points.dat',
    ]
    for pat in patterns:
        for f in _glob.glob(os.path.join(WORK_DIR, pat)):
            try:
                os.remove(f)
            except OSError:
                pass
    # Also clean experiment subdirectories
    for exp_name in EXPERIMENTS:
        exp_dir = os.path.join(WORK_DIR, exp_name)
        if os.path.isdir(exp_dir):
            import shutil
            shutil.rmtree(exp_dir)


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

# DPLB — Deformable Particle Lattice Boltzmann

A hybrid Lattice Boltzmann Method (LBM) and Molecular Dynamics (MD) simulation code for modeling the flow of soft particles such as polymers and biological cells in confined geometries.

**ASLB** (the internal name) extends the [Susp3D](http://www.che.ufl.edu/ladd/) particle suspension code (Tony Ladd, 2003) with deformable particle and polymer capabilities.

## Features

- **D3Q19 Lattice Boltzmann** fluid solver with BGK collision in mode space
- **Deformable particles** — icosahedral surface meshes (12 to 2562 beads) with spring networks, bending rigidity, and volume/area constraints; suitable for modeling red blood cells and capsules
- **Polymers** — bead-spring chains with FENE, Worm-Like Chain (WLC), or harmonic bonds
- **Deformable tubes/channels** — soft sinusoidal tubes, slits, and trapezoidal channels
- **Rigid colloids** — spheres and walls with lubrication corrections
- **Fluid-solid coupling** — immersed boundary method with bounce-back no-slip condition and implicit velocity update
- **Many-body lubrication** — asymptotic Happel-Brenner formulas with cluster-based implicit solver
- **Excluded volume interactions** — hard-sphere, WCA (truncated Lennard-Jones), and Gaussian soft repulsion
- **Thermal fluctuations** — fluctuating LBM and Langevin thermostat for Brownian dynamics
- **Parallelization** — MPI domain decomposition (x-direction) and OpenMP threading

## Requirements

- C compiler (GCC, ICC, etc.)
- MPI library (OpenMPI, MPICH, Intel MPI)
- Intel MKL (for the VSL random number generator)
- OpenMP support

## Building

Compile all `.c` files in `MASTER/` and link against MPI and MKL:

```bash
cd MASTER
mpicc -O2 -fopenmp -o dplb *.c -lmkl_rt -lm
```

Adjust compiler flags and library paths for your system. For Intel compilers:

```bash
mpiicc -O2 -qopenmp -o dplb *.c -mkl -lm
```

## Input Files

The simulation reads input from the `init/` directory relative to the working directory:

| File | Description |
|------|-------------|
| `init/task_file.dat` | Task number (integer) for labeling output files |
| `init/input_file.dat` | Global simulation parameters (see below) |
| `init/c_inp.<task>` | Particle configuration: count, domain size, positions, radii |
| `init/p_inp.dat` | Deformable particle / polymer parameters |
| `init/alloc_file.dat` | Processor allocation |

### Global parameters (`input_file.dat`)

```
num_cycle  num_step  num_modes  t_lbe
wall_flag  backflow_flag  add_noise
rho  ux  uy  uz
mass_fac  vel_fac
lub_cut  del_hy
f_sph_x  f_sph_y  f_sph_z
f_ext_x  f_ext_y  f_ext_z
u_wall_x  u_wall_y  u_wall_z
tau  tau_v  tau_g              (fluid 1 relaxation times)
tau  tau_v  tau_g              (fluid 2 relaxation times)
lub_N  lub_T  lub_R
seed
```

Key parameters:
- `wall_flag`: 0 = fully periodic, 1 = y-walls, 2 = y+z walls, 3 = x+y+z walls
- `backflow_flag`: 0 = off, 1 = on, 2 = on + fluid at rest
- `tau`: BGK relaxation time (controls kinematic viscosity: nu = (tau - 0.5) * cs^2)
- `add_noise`: 0 = deterministic, >0 = fluctuating hydrodynamics

### Deformable particle parameters (`p_inp.dat`)

Includes: number of particles, particle types, triangulation level, spring type (0=FENE, 1=WLC, 2=harmonic), spring constants, bending rigidity (`k_bend`), volume constraint (`k_V`), area constraint (`k_A`), excluded volume type, friction coefficient, temperature (`kT`), Verlet timestep, and output frequencies.

## Running

```bash
# Serial
./dplb [work_dir]

# Parallel (e.g., 4 MPI ranks)
mpirun -np 4 ./dplb [work_dir]
```

If `work_dir` is omitted, the current directory is used. The code expects `init/` and `data/` subdirectories under the working directory.

## Output Files

Output is written to the `data/` directory:

| File | Description |
|------|-------------|
| `data/properties.dat` | Echoed input parameters |
| `data/p_out.<task>` | Particle trajectories (position, velocity, displacement) |
| `data/m_out.<task>` | Monomer positions and deformable particle properties |
| `data/chk_p.<task>` | Particle checkpoint (restart configuration) |
| `data/chk_f.<task>` | Fluid field checkpoint |
| `data/final.config` | Final particle configuration |
| `data/modes.<task>` | Fluid mode output |

## Usage

### Directory Setup

Create a working directory with the required subdirectories:

```bash
mkdir -p myrun/init myrun/data
```

### Step 1: Create the Task File

```bash
# init/task_file.dat — a single integer identifying this run
echo "1" > myrun/init/task_file.dat
```

### Step 2: Create the Global Input File

`init/input_file.dat` controls the fluid solver and overall simulation. Each data line is preceded by a label line (read and discarded). Example for Poiseuille flow in a slit channel:

```
num_cycle num_step num_modes t_lbe
100       1000     0         1.0
wall_flag backflow_flag add_noise
1         1              0
rho       ux        uy        uz
1.0       0.0       0.0       0.0
mass_fac  vel_fac
1.0       1.0
lub_cut   del_hy
0.5       0.0
f_sph_x   f_sph_y   f_sph_z
0.0       0.0       0.0
f_ext_x   f_ext_y   f_ext_z
1.0e-5    0.0       0.0
u_wall_x  u_wall_y  u_wall_z
0.0       0.0       0.0
tau       tau_v     tau_g
1.0       1.0       1.0
tau       tau_v     tau_g
1.0       1.0       1.0
lub_N     lub_T     lub_R
2.0       2.0       2.0
seed
12345
```

Parameter reference:

| Parameter | Description |
|-----------|-------------|
| `num_cycle` | Number of output cycles |
| `num_step` | LBM timesteps per cycle (total steps = num_cycle x num_step) |
| `num_modes` | Number of fluid modes to output (0 = none) |
| `t_lbe` | LBM timestep duration (MD steps per LBM step = t_lbe / dt) |
| `wall_flag` | Wall geometry: 0 = periodic, 1 = slit (y-walls), 2 = rectangular channel (y+z walls), 3 = box (all walls) |
| `backflow_flag` | Backflow correction: 0 = off, 1 = on, 2 = on + fluid at rest, 3/4/5 = single-axis only (x/y/z) |
| `add_noise` | Thermal fluctuations in fluid: 0 = off, 1 = on |
| `rho` | Fluid density ratio (reference density = 36 in lattice units) |
| `ux, uy, uz` | Initial fluid velocity |
| `mass_fac` | Particle mass scaling factor |
| `vel_fac` | Particle velocity scaling factor |
| `lub_cut` | Lubrication force cutoff distance |
| `del_hy` | Hydrodynamic layer thickness |
| `f_sph` | Particle body acceleration (x, y, z) |
| `f_ext` | External fluid body force density (x, y, z) |
| `u_wall` | Wall velocity for shear/channel flow |
| `tau` | BGK relaxation time; kinematic viscosity nu = (tau - 0.5)/3 |
| `tau_v` | Viscous relaxation time |
| `tau_g` | Ghost relaxation time |
| `lub_N, lub_T, lub_R` | Lubrication ranges (normal, tangential, rotational) |
| `seed` | Random number seed (overridden by system clock) |

### Step 3: Create the Particle Configuration File

`init/c_inp.<task>` (e.g., `init/c_inp.1`) defines the domain and rigid particles. The first line sets the number of rigid colloids and the domain dimensions:

```
num_sph  max_x  max_y  max_z
0        64     32     32
```

For each rigid colloid, add two lines:

```
x   y   z   ex  ey  ez          (position and orientation)
ux  uy  uz  wx  wy  wz          (linear and angular velocity)
shape_flag  mass_flag  move_flag  rho  r_a  r_b
```

Where:
- `shape_flag`: 0 = sphere, 1 = ellipsoid, 2 = open cylinder, 3 = closed cylinder, 4 = capped cylinder
- `mass_flag`: 0 = infinite mass (fixed), 1 = finite mass
- `move_flag`: 0 = stationary, 1 = mobile, 2 = wall-attached
- `rho`: density ratio (particle/fluid)
- `r_a, r_b`: semi-major and semi-minor axes

Example with one sphere at the domain center:

```
1  64  32  32
16.0  16.0  16.0  0.0  0.0  0.0
0.0   0.0   0.0   0.0  0.0  0.0
0  1  1  1.0  5.0  5.0
```

### Step 4: Create the Deformable Particle / Polymer Input File

`init/p_inp.dat` defines deformable spheres, polymers, and tubes. The file has three sections read sequentially:

**Section 1 — Deformable spheres:**

```
Ntype0 Ntype1 nlevel0 nlevel1 max_x max_y max_z       (header line skipped)
2      0      1       0       64    32    32
spring ev_type verlet initconfig                       (header line skipped)
0      1       1      0
H0 H1 Q0 Q1 sigma nks kb0 kb1 kV0 kV1 kA0 kA1        (header line skipped)
100.0 0.0 1.5 1.0 1.0 1.0 0.01 0.0 50.0 0.0 50.0 0.0
evcutoff eps fric dt monmass kT                        (header line skipped)
2.5 1.0 5.0 0.005 1.0 0.0
f_ext_x f_ext_y f_ext_z                                (header line skipped)
0.0 0.0 0.0
write_time write_config write_fluid                    (header line skipped)
1000 10000 100000
```

| Parameter | Description |
|-----------|-------------|
| `Ntype0, Ntype1` | Number of type-0 and type-1 deformable spheres |
| `nlevel` | Triangulation level: 0 = 12 beads, 1 = 42, 2 = 162, -1 = 162 (template), -2 = 642 (RBC template) |
| `spring` | Bond type: 0 = FENE, 1 = WLC, 2 = harmonic |
| `ev_type` | Excluded volume: 0 = hard-sphere, 1 = WCA, 2 = Gaussian |
| `verlet` | Integration scheme: 0 = forward Euler, 1 = velocity Verlet (recommended), 2 = leapfrog |
| `initconfig` | Initial placement: 0 = random, 3 = read from file |
| `H_fene` | Spring constant |
| `Q_fene` | Maximum spring extension ratio |
| `k_bend` | Bending rigidity |
| `k_V` | Volume constraint strength |
| `k_A` | Area constraint strength |
| `evcutoff` | Excluded volume cutoff distance |
| `eps` | Lennard-Jones energy parameter |
| `fric` | Stokes friction coefficient for a point force |
| `dt` | MD timestep (must satisfy t_lbe/dt = integer) |
| `monmass` | Monomer mass |
| `kT` | Thermal energy (0 = athermal) |
| `write_time` | Output properties every N steps |
| `write_config` | Output configuration every N steps |
| `write_fluid` | Output fluid field every N steps |

**Section 2 — Polymer chains** (appended after sphere section):

```
Ntype0 Ntype1 Nbeads0 Nbeads1                          (header line skipped)
0      0      0       0
spring ev_type initconfig                              (header line skipped)
0      1       1
H0 H1 Q0 Q1 sigma nks kb0 kb1                         (header line skipped)
30.0 0.0 1.3 1.0 1.0 1.0 0.0 0.0
evcutoff eps fric dt monmass kT                        (header line skipped)
2.5 1.0 5.0 0.005 1.0 0.001
f_ext_x f_ext_y f_ext_z                                (header line skipped)
0.0 0.0 0.0
```

| Parameter | Description |
|-----------|-------------|
| `Ntype0, Ntype1` | Number of type-0 and type-1 polymer chains |
| `Nbeads0, Nbeads1` | Beads per chain for each type |
| `initconfig` | 1 = random walk, 3 = read from file |

**Section 3 — Tube/channel** (appended after polymer section):

```
NDP spring ev_type initconfig                          (header lines skipped)
0   0      0       1
radius H_fene Q_fene k_bend k_V k_A                   (header line skipped)
10.0 100.0 1.5 0.01 0.0 0.0
evcutoff fric monmass kT eps                           (header line skipped)
2.5 5.0 1.0 0.0 1.0
ramp rperiod offset                                    (header line skipped)
2.0 32.0 0
```

| Parameter | Description |
|-----------|-------------|
| `NDP` | Number of tube segments (0 = no tube) |
| `initconfig` | 1 = soft sinusoidal tube, 2 = slit with contraction, 3 = curved slit, 5 = cosine slit, 6 = trapezoid |
| `radius` | Tube radius (or slit half-height) |
| `ramp` | Amplitude of radius variation (config 1) or contraction angle in degrees (config 2) |
| `rperiod` | Wavelength of radius variation (config 1) or central flat length (config 2) |
| `offset` | Initial offset region for tube geometry |

### Step 5: Processor Allocation (optional)

`init/alloc_file.dat` assigns x-slices to each MPI rank. One integer per line, one line per processor. The sum must equal `max_x`. If omitted, slices are divided evenly.

```
16
16
16
16
```

### Step 6: Run the Simulation

```bash
# Serial
./dplb myrun

# Parallel (4 MPI ranks)
mpirun -np 4 ./dplb myrun
```

### Step 7: Examine Output

All output goes to `myrun/data/`:

| File | Contents | Frequency |
|------|----------|-----------|
| `properties.dat` | Echoed parameters, derived quantities (viscosity, fluctuation variance, etc.) | Once at startup |
| `p_out.<task>` | Per-cycle: rigid particle coordinates, velocities, forces, wall forces | Every cycle |
| `m_out.<task>` | Per-cycle: fluid mass/momenta profiles, Reynolds/viscous/particle/collisional stress | Every cycle |
| `data/mon_*.dat` | Monomer positions, velocities, and per-DP properties (Rg, volume, area, asphericity) | Every `write_time` steps |
| `data/stress_*.dat` | Particle-fluid stress tensor | Every `write_time` steps |
| `data/config_*.dat` | Full bead configurations for restart | Every `write_config` steps |
| `data/fluid_*.dat` | 3D fluid velocity/density field | Every `write_fluid` steps |
| `data/nodemap_*.dat` | 3D node map (fluid/solid/boundary) | Every `write_fluid` steps |
| `chk_p.<task>` | Particle checkpoint for restart | End of each cycle |
| `chk_f.<task>` | Fluid field checkpoint for restart | End of each cycle |
| `final.config` | Final particle configuration | End of simulation |

### Example: RBC in Poiseuille Flow

Set up a single red blood cell (RBC) template flowing through a slit channel:

```bash
mkdir -p rbc_flow/init rbc_flow/data

# Task
echo "1" > rbc_flow/init/task_file.dat

# Domain: 128 x 32 x 32 with y-walls, no rigid colloids
cat > rbc_flow/init/c_inp.1 << 'EOF'
0  128  32  32
EOF

# Global parameters: slit channel, pressure-driven flow
cat > rbc_flow/init/input_file.dat << 'EOF'
num_cycle num_step num_modes t_lbe
50        2000     0         1.0
wall_flag backflow_flag add_noise
1         1              0
rho       ux        uy        uz
1.0       0.0       0.0       0.0
mass_fac  vel_fac
1.0       1.0
lub_cut   del_hy
0.5       0.0
f_sph_x   f_sph_y   f_sph_z
0.0       0.0       0.0
f_ext_x   f_ext_y   f_ext_z
5.0e-6    0.0       0.0
u_wall_x  u_wall_y  u_wall_z
0.0       0.0       0.0
tau       tau_v     tau_g
1.0       1.0       1.0
tau       tau_v     tau_g
1.0       1.0       1.0
lub_N     lub_T     lub_R
2.0       2.0       2.0
seed
12345
EOF
```

Then create `rbc_flow/init/p_inp.dat` with 1 type-0 RBC (nlevel=-2 for the 642-bead template), appropriate spring/bending/volume/area constants, zero polymers, and zero tubes. Run with:

```bash
mpirun -np 4 ./dplb rbc_flow
```

### Wall Configurations

The `wall_flag` parameter selects the confinement geometry:

| wall_flag | Geometry | Walls | Typical flow |
|-----------|----------|-------|--------------|
| 0 | Fully periodic | None | Bulk suspension |
| 1 | Slit channel | y-, y+ | Poiseuille, Couette (set `u_wall.x`) |
| 2 | Rectangular channel | y-, y+, z-, z+ | Duct flow, elongational (set `u_wall.y`) |
| 3 | Closed box | All 6 faces | Sedimentation, confined systems |

### Flow Driving Mechanisms

| Method | Parameters | Use case |
|--------|-----------|----------|
| Pressure-driven (Poiseuille) | `f_ext` != 0, `u_wall` = 0 | Channel/pipe flow |
| Shear (Couette) | `u_wall.x` != 0, `wall_flag` = 1 | Rheology, shear flow |
| Gravity/sedimentation | `f_sph` != 0 with `backflow_flag` = 1 | Settling particles |
| Elongational | `u_wall.y` != 0, `wall_flag` = 2, `max_y` = `max_z` | Extensional flow |

## Code Structure

```
MASTER/
  main.c              Entry point
  driver.c             Initialization, main time loop, I/O
  header.h             Global constants, includes
  struct.h             Data structures (object, monomer, face, DP, DP_param)
  lbe.h                D3Q19 lattice vectors, eigenvector matrices
  macro.h              Utility macros (periodic images, range checks)
  func_proto.h         Function prototypes

  initlbe.c            LBM eigenvector setup and equilibrium initialization
  lbe.c                Streaming, boundary exchange, periodic/wall BCs
  lbe_update.c         Full LBM cycle: boundary mapping + collision + streaming

  objects_init.c        Rigid particle and wall initialization
  objects_map.c         Shape discriminator functions (sphere, cylinder)
  init_sphere.c         Deformable sphere (icosahedron) mesh generation
  init_polymer.c        Polymer chain initialization
  init_cylinder.c       Deformable tube/channel geometry

  bnodes.c             Particle-fluid force exchange (bounce-back)
  bnodes_init.c         Boundary node identification and friction tensors
  bnodes_dp.c           Deformable particle boundary mapping
  bnodes_tube.c         Tube/channel boundary conditions

  get_forces.c          Force computation: springs, EV, drag, bending, constraints
  implicit_force.c      Implicit velocity update (removes CFL restriction)
  lub.c                Lubrication force tensors (sphere-sphere, sphere-wall)
  cluster_force.c       Many-body lubrication in particle clusters

  update.c             Orchestrates LBM + MD time stepping
  verlet_update.c       Velocity Verlet integration for deformable particles
  velcs_update.c        Rigid particle velocity update with lubrication

  clusters.c           Cluster identification (neighbor graph labeling)
  cluster_update.c      Implicit cluster velocity solve
  n_list.c             Cell-linked neighbor lists
  hs3d.c               Event-driven hard-sphere collision detection

  msg_mpi.c            MPI communication (domain exchange, global reductions)
  msg_ser.c            Serial communication stubs
  global_sums.c         MPI global reduction routines

  output.c             Data output (trajectories, fluid fields, stress)
  utils.c              File naming, error handling
  ran_num.c            L'Ecuyer long-period random number generator
  cj_grad.c            Conjugate gradient / matrix utilities
  modes_write.c         Fluid mode output
```

## Key Physical Models

### Lattice Boltzmann
- D3Q19 velocity model, cs^2 = 1/3
- BGK collision in moment space with tunable relaxation times
- Optional thermal fluctuations (fluctuating LBM)

### Deformable Particles
- Icosahedral triangulation with configurable refinement levels (12 to 2562 vertices)
- FENE springs: F = -H r / (1 - (r/Q)^2)
- Bending resistance between adjacent triangles
- Global volume constraint: F_V = -k_V (V - V0) / V0
- Global area constraint: F_A = -k_A (A - A0) / A0

### Polymer Chains
- Bead-spring chains with FENE, WLC, or harmonic bonds
- Excluded volume (WCA, hard-sphere, or Gaussian)
- Hydrodynamic coupling via Stokes drag from LBM fluid

### Time Integration
- Velocity Verlet (multiple variants) for molecular dynamics
- Implicit fluid-particle coupling for numerical stability

## License

GNU General Public License v2.0 or later. See the license headers in source files.

## Author

Yeng-Long Chen (2019), based on Susp3D by Tony Ladd (2003).

# Preparation of Basic VASP Input Files

## Required Files:

- **INCAR** (Instructs VASP what to calculate and how to calculate)
- **KPOINTS** (Contains K-point information)
- **POSCAR** (Describes the structure of the model, i.e., the positions of atoms in xyz coordinates)
- **POTCAR** (Contains the pseudopotentials for the atoms in the system, including information about atomic nuclei and electrons)

---

## 1. Key Concepts

- **Rubbish in, Rubbish out!**: The program only performs calculations; it is up to the user to ensure the inputs are correct. There are three main types of errors:

  - **Model error**: Mistakes in modeling, mainly related to `POSCAR`.
  - **Calculation parameters**: Mistakes in `INCAR`, `KPOINTS`, or `POTCAR`.
  - **Submission script or command errors**.

- **Check parameters repeatedly** on the official VASP website.

---

## 2. VASP Input Files

### 2.1 INCAR

#### Purpose:

Instructs VASP what to calculate and how to calculate it.

#### Preparation Principle:

Keep it as simple as possible. Do not include what you don’t understand.

#### Example of INCAR File:

```ini
SYSTEM = O atom # This refers to Oxygen (O), not the digit zero
ISMEAR = 0      # This is zero, not O (Oxygen); use 0 for molecules or atoms
SIGMA = 0.01    # Use 0.01 for molecules or atoms
```

#### Detailed Explanation:

- **SYSTEM**: Describes the task of the calculation (purpose, system).
- **ISMEAR**:
  - Different values correspond to different smearing methods.
  - For molecules or atoms, use 0.
- **SIGMA**:
  - Related to ISMEAR.
  - When ISMEAR = -5, SIGMA can be ignored.
  - For metals: ISMEAR = 1 or 0, for non-metals: ISMEAR = 0, generally use SIGMA = 0.01.
  - For gases or atomic systems (molecule/atom in a box), use ISMEAR = 0 and SIGMA = 0.01.

#### Testing Criterion:

Check that the average `entropy T*S` per atom in the `OUTCAR` file is less than 1-2 meV using the command:

```bash
grep 'entropy T' OUTCAR
```

Divide the energy by the number of atoms and compare with 0.001 eV. If it’s smaller, SIGMA is fine; if larger, reduce SIGMA and test again.

---

### 2.2 KPOINTS

#### Purpose:

Determines the precision of the calculation and influences the required calculation time. Higher precision means longer calculation time, while lower precision shortens the time.

#### Example of KPOINTS File:

```
K-POINTS  # First line can contain anything but must not be empty
0         # Zero, auto-generate the grid
Gamma     # Gamma-centered grid
1 1 1     # 1x1x1 grid
0 0 0     # S1 S2 S3 (usually kept at 0 0 0)
```

#### Detailed Explanation:

- **Gamma-centered** grids are used in most cases.
- **Monkhorst-Pack (M-centered)** grids shift the grid by 1/(2N) in each direction, which can lead to mismatches between the grid and the crystal symmetry, especially for hexagonal systems. Therefore, use Gamma-centered grids for hexagonal structures.

---

### 2.3 POSCAR

#### Purpose:

Contains the structural information of the model for calculation.

#### Example of POSCAR File (for an O atom in a box):

```
O atom in a box
1.0  # Universal scaling factor
8.0 0.0 0.0  # Lattice vector a(1)
0.0 8.0 0.0  # Lattice vector a(2)
0.0 0.0 8.0  # Lattice vector a(3)
O           # Oxygen element (O)
1           # Number of atoms
Cartesian   # Cartesian coordinate system
0 0 0       # Atom position in Cartesian coordinates
```

#### Detailed Explanation:

- **First line**: Description of the system (can be anything).
- **Second line**: Scaling factor, usually set to 1.0.
- **Third to Fifth lines**: Lattice vectors of the system.
- **Sixth line**: Element name, in this case, Oxygen (O).
- **Seventh line**: Number of atoms, here it’s 1 atom.
- **Eighth line**: Specifies the coordinate system. "Cartesian" means the positions are in Cartesian coordinates; "Direct" means fractional coordinates.
- **Ninth line**: Atom positions. In this example, the O atom is placed at the origin (0.0, 0.0, 0.0).

---

### 2.4 POTCAR

#### Purpose:

Contains pseudopotential information about each atom in the system, including the atomic nucleus and electron configuration.

#### Example of POTCAR:

```
PAW_PBE Fe 06Sep2000
8.00000000000000000
parameters from PSCTR are:
VRHFIN =Fe: d7 s1
LEXCH = PE
EATOM = 594.4687 eV, 43.6922 Ry
...
```

#### Key Parameters:

- **VRHFIN**: Electron configuration of the atom.
- **LEXCH**: Specifies the exchange-correlation functional, here GGA-PBE.
- **TITEL**: Describes the element and the date of the pseudopotential.
- **ZVAL**: Number of valence electrons, critical for Bader charge analysis.
- **ENMAX**: Default cutoff energy, related to the `ENCUT` parameter in `INCAR`.

#### Linux Commands for Checking POTCAR:

- View the element in POTCAR:

  ```bash
  grep TIT POTCAR
  ```

- View the cutoff energy:

  ```bash
  grep ENMAX POTCAR
  ```

- View the number of valence electrons:
  ```bash
  grep ZVAL POTCAR
  ```

---

## References:

- [VASP File Structure](http://cms.mpi.univie.ac.at/vasp/guide/node50.html)
- [INCAR Documentation](http://cms.mpi.univie.ac.at/vasp/guide/node91.html)
- [KPOINTS Documentation](https://cms.mpi.univie.ac.at/vasp/vasp/Automatic_k_mesh_generation.html)
- [POSCAR Documentation](http://cms.mpi.univie.ac.at/vasp/guide/node59.html)
- [POTCAR Documentation](http://cms.mpi.univie.ac.at/vasp/vasp/Recommended_PAW_potentials_DFT_calculations_using_vasp_5_2.html)

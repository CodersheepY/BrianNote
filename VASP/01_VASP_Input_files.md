# Preparation of Basic VASP Input Files

**Required files:**

- `INCAR` (tells VASP what to calculate and how to calculate)
- `KPOINTS` (contains K-point information for the calculation)
- `POSCAR` (describes the structure of the model in xyz coordinates)
- `POTCAR` (describes the basis sets of atoms in the system)

---

## 1. Key Concepts

- **Rubbish in, Rubbish out!**: The program only computes; it's up to the user to ensure the inputs are correct. Errors usually fall into three categories:

  - **Model errors:** Errors in modeling (usually in `POSCAR`)
  - **Calculation parameters:** Errors in `INCAR`, `KPOINTS`, or `POTCAR`
  - **Task submission script or command errors**

- **Check parameters repeatedly** from the official VASP site.

---

## 2. VASP Input Files

### 2.1 INCAR

#### Purpose:

Tells VASP what to calculate and how to calculate it.

#### Preparation Principle:

Keep it as simple as possible. Do not include what you don't understand.

#### Sample INCAR File:

```ini
SYSTEM = O atom # The element is oxygen, not the digit zero
ISMEAR = 0      # Use 0 for molecules or atoms, not the letter O
SIGMA = 0.01    # For molecules or atoms, use 0.01
```

#### Explanation:

- **SYSTEM**: Description of the task (objective, system)
- **ISMEAR**:
  - Different values correspond to different smearing methods.
  - Use `0` for molecules or atoms.
- **SIGMA**:
  - Value related to ISMEAR.
  - For metals: ISMEAR = 1 or 0, non-metals: ISMEAR = 0. Use `SIGMA = 0.01` for non-metals or atoms.
  - Check: Use `grep 'entropy T' OUTCAR` to find entropy. Ensure that `T*S` divided by the number of atoms is less than 1-2 meV. Adjust `SIGMA` accordingly.

### 2.2 KPOINTS

#### Purpose:

Determines the precision of the calculation, which affects the time needed. High precision requires longer time, lower precision requires shorter time.

#### Sample KPOINTS File:

```
K-POINTS  # First line can be anything but must not be empty
0         # Use automatic generation
Gamma     # Gamma-centered grid
1 1 1     # 1*1*1 grid
0 0 0     # S1 S2 S3 (typically kept at 0 0 0)
```

#### Explanation:

- **Gamma-centered** is usually recommended for molecules, atoms, and hexagonal systems.
- **M-centered** (Monkhorst-Pack) is another grid option but can mismatch with hexagonal symmetry.

### 2.3 POSCAR

#### Purpose:

Contains structural information of the model.

#### Sample POSCAR File (for O atom in a box):

```
O atom in a box
1.0  # Universal scaling factor
8.0 0.0 0.0  # Lattice vector a(1)
0.0 8.0 0.0  # Lattice vector a(2)
0.0 0.0 8.0  # Lattice vector a(3)
O           # Element (Oxygen, not zero)
1           # Number of atoms
Cartesian   # Coordinate system (use Direct for fractional)
0 0 0       # Atomic position in Cartesian coordinates
```

#### Explanation:

- The coordinate system can be Cartesian or fractional. For Cartesian, use "Cartesian" or "C"; for fractional, use "Direct" or "D".
- Coordinates of atoms are placed after the coordinate system.

### 2.4 POTCAR

#### Purpose:

Contains basis set information for atoms in the system, including core radius, cutoff energies, and electron configurations.

#### Example Output of POTCAR:

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

- **VRHFIN**: Electron configuration of the element.
- **LEXCH**: Specifies GGA-PBE functional.
- **TITEL**: Element and potential creation date.
- **ZVAL**: Number of valence electrons (important for Bader charge analysis).
- **ENMAX**: Default cutoff energy, related to ENCUT in INCAR.

#### Common Linux Commands to Check POTCAR:

- View the element in POTCAR: `grep TIT POTCAR`
- View the cutoff energy: `grep ENMAX POTCAR`
- View the valence electron number: `grep ZVAL POTCAR`

---

## References:

- [VASP Documentation - File Structure](http://cms.mpi.univie.ac.at/vasp/guide/node50.html)
- [INCAR Documentation](http://cms.mpi.univie.ac.at/vasp/guide/node91.html)
- [KPOINTS Documentation](https://cms.mpi.univie.ac.at/vasp/vasp/Automatic_k_mesh_generation.html)
- [POSCAR Documentation](http://cms.mpi.univie.ac.at/vasp/guide/node59.html)
- [POTCAR Documentation](http://cms.mpi.univie.ac.at/vasp/vasp/Recommended_PAW_potentials_DFT_calculations_using_vasp_5_2.html)

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

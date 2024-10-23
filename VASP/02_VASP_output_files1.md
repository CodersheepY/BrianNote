# VASP Output Files

## 1. Difference Between `CONTCAR` and `POSCAR`

In VASP, `POSCAR` is used to store the structural information of the model. After the structure optimization is completed, VASP generates a new structure file, `CONTCAR`, which stores the optimized structure. For a single atom, even though there is no structural optimization, VASP will still generate a `CONTCAR`. If no optimization occurs, the `CONTCAR` and `POSCAR` will be identical. However, the contents of these two files may or may not differ.

### Example Comparison of `POSCAR` and `CONTCAR` (Oxygen Atom):

- **Indentation**: In both `CONTCAR` and `POSCAR`, the only allowed separator between elements on each line is spaces.
- **Coordinate System**: In `POSCAR`, the Cartesian coordinate system (denoted as `Cartesian` or `C`) can be used, while `CONTCAR` usually defaults to fractional coordinates (`Direct`).
  - Insight: To convert a system from Cartesian to Direct, you can run a single-point calculation.
- **Cartesian vs Direct**: Both systems are allowed in `POSCAR`, but Cartesian is more intuitive.
- **Origin Atom**: In this example, since the O atom is located at the origin (0, 0, 0), there's no apparent difference between Cartesian and Direct coordinates.
- **Extra Line in `CONTCAR`**: `CONTCAR` contains an additional line, which shows the movement of atoms in the xyz directions. For our example (a single atom), there’s only one extra line, all filled with zeroes (0 0 0).
  - If you’re performing structural optimization, there will be as many additional lines as there are atoms in the system. In molecular dynamics, these lines will represent the velocities of the atoms. When using the dimer method to calculate transition states, these lines relate to the vibrational modes of the transition state.
  - Summary: If you are not running molecular dynamics or transition state calculations, these extra lines are all zeros and can be safely removed without affecting the calculation.

---

## 2. Role of `CONTCAR`

- **Same as `POSCAR`**: `CONTCAR` stores the structure of the optimized model, which may differ from the manually written `POSCAR`.
- **Two Uses for `CONTCAR`**:
  1. After a calculation is completed, `CONTCAR` contains the structure from the last optimization step. You can open this file with visualization software to check the accuracy of the calculation.
  2. If the calculation is interrupted (due to server issues, power outages, or accidental file deletions), `CONTCAR` can be used to resume the calculation. (I don’t fully understand how to continue a calculation yet, so this is left as a "to be explored" topic.)

---

## 3. Two Optimization Processes in VASP

- **Electronic Structure Optimization**: For a fixed geometry, the Schrödinger equation is solved iteratively to find the energy minimum of the system. This process is called electronic structure optimization (also known as an "electronic step").
- **Geometrical Structure Optimization**: Based on the forces acting on the atoms from the optimized electronic structure, the atomic positions are adjusted iteratively until the system reaches a minimum on the potential energy surface (also known as an "ionic step").
  - **Example**: For an O atom, only electronic structure optimization is needed because its geometry doesn’t change.

---

## 4. OSZICAR

- **Purpose**: The iterative process of solving the electronic structure is recorded in the `OSZICAR` file. This includes both electronic and geometrical optimization steps.
  - **Official Explanation**: "Information about convergence speed and about the current step is written to stdout and to the `OSZICAR` file. Always keep a copy of the `OSZICAR` file, it might give important information." (Source: [OSZICAR](https://cms.mpi.univie.ac.at/wiki/index.php/OSZICAR))

#### Example Output of OSZICAR:

```bash
  N       E                     dE             d eps       ncg     rms          rms(c)
DAV:  1   0.324969965196E+02   0.32497E+02   -0.10270E+03    48   0.977E+01
DAV:  2   0.501749892771E+00  -0.31995E+02   -0.31995E+02    72   0.202E+01
DAV:  3  -0.182605770767E-01  -0.52001E+00   -0.50521E+00    48   0.521E+00
DAV:  4  -0.203547758465E-01  -0.20942E-02   -0.20860E-02    96   0.333E-01
DAV:  5  -0.203547873947E-01  -0.11548E-07   -0.11210E-07    48   0.844E-04  0.307E-01
DAV:  6  -0.213726161828E-01  -0.10178E-02   -0.17884E-03    48   0.111E-01  0.155E-01
DAV:  7  -0.214708381542E-01  -0.98222E-04   -0.23522E-04    48   0.459E-02
 1 F= -0.21470838E-01  E0= -0.13757722E-01  d E = -0.154262E-01
```

#### Explanation of Each Term:

1. **N**: Number of electronic structure iteration steps (called "electronic steps").
2. **E**: Energy of the current electronic step.
3. **dE**: Difference in energy between the current and previous electronic steps.
4. **d eps**: Change in the band structure energy.
5. **ncg**: Number of evaluations of the Hamiltonian acting on a wavefunction.
6. **rms**: Norm of the residuum of the trial wavefunctions (i.e., the approximate error).
7. **rms (c)**: Difference between input and output charge density.

---

### Blocked Davidson Algorithm (DAV):

- **DAV**: Represents the Blocked Davidson algorithm, a self-consistent algorithm for solving electronic structure iteratively.

  - Other algorithms include **RMM** (Residual Minimization Scheme) and **CG** (Conjugate Gradient Algorithm).
  - To set a specific algorithm in VASP, use the `ALGO` parameter in `INCAR`. For most cases, `ALGO = Fast` suffices (official reference: [ALGO](https://cms.mpi.univie.ac.at/wiki/index.php/ALGO)).

- **Why does DAV appear in OSZICAR even if it wasn’t set in INCAR?**
  - VASP has many default parameters to avoid the calculation getting stuck due to missing settings. If no specific algorithm is set in INCAR, VASP defaults to using DAV (equivalent to `ALGO = N`).

---

### The Last Line in OSZICAR:

1. **F**: Number of ionic steps (geometry optimization steps). In this example, there’s only one step.
2. **F =**: Total energy of the system, equivalent to `free energy TOTEN` in the `OUTCAR`.
3. **E0**: The energy corresponding to `energy (sigma->0)` in the `OUTCAR`.

---

## 5. Role of OSZICAR

- The `OSZICAR` file records the optimization process of the entire system. Though `OUTCAR` also contains detailed information, `OSZICAR` provides a more straightforward view of the energy changes during optimization.
- To get the system’s energy from `OSZICAR`, use the following command:
  ```bash
  grep E0 OSZICAR
  ```

---

## Summary:

1. **Optimization in VASP**: Understanding the basics of electronic and geometrical structure optimization.
2. **OSZICAR Fields**: Familiarity with the different terms in `OSZICAR`.
3. **ALGO Parameter in INCAR**: `ALGO = Fast` can satisfy most calculation needs.
4. **Default Parameters in VASP**: Knowing that VASP has defaults to prevent calculation failure due to missing settings.
5. **Retrieving Energy from OSZICAR**: Knowing how to get the system's energy from `OSZICAR`.

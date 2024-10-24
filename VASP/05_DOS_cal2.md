# DOS Calculations Are Done in Two Steps

The VASP official guide mentions two methods for calculating DOS. One is to perform a self-consistent calculation (commonly referred to as static or SCF calculation), and then process the `DOSCAR` and `vasprun.xml` files, as shown in the example:

### Prerequisites:

1. High-quality DOS requires a fine K-point mesh. However, setting too many K points can become computationally expensive, as discussed previously regarding the relationship between K-points and computational time. When computational resources were limited, performing DOS calculations in two steps was a common way to solve this problem.
2. Another reason is related to band structure calculations. As mentioned in an online explanation: during band structure calculations, the K-points are positioned along high-symmetry lines in reciprocal space, where self-consistent calculations are not possible. Refer to this [article](http://blog.sciencenet.cn/blog-567091-675253.html/). For band structure calculations, self-consistent calculations are a necessary step.
3. Even when the number of K-points is increased, charge density and effective potential converge quickly, meaning that variations in K-points have little effect on the convergence of the charge density.

### Analysis:

For now, setting band structure calculations aside and focusing on points 1 and 3: After structure optimization, when calculating DOS:

- **Step 1**: Use a smaller number of K-points to perform a single-point calculation and generate the `CHGCAR` file.
- **Step 2**: Read the `CHGCAR` file from the first step (`ICHARG=11`).

This approach avoids the computational burden of using a high-density K-point mesh directly. The two-step process for DOS calculations is mainly about **saving time**. Performing DOS calculations in two steps is **not mandatory**. A high-density K-point mesh can be used for a one-step calculation. However, for **band structure calculations**, a two-step process is required.

Given modern computational resources, performing a one-step DOS calculation with a high K-point mesh is entirely feasible.

## LDOS and PDOS

- **LORBIT = 10**: Decomposes the density of states for each atom and its `spd` orbitals, referred to as the **Local DOS (LDOS)**.
- **LORBIT = 11**: On top of what LORBIT=10 provides, further decomposes the orbitals into `px`, `py`, `pz`, referred to as the **Projected DOS (PDOS)** or **Partial DOS (PDOS)**.

Using **LORBIT = 11** provides more detailed information.

## WAVECAR

If a `WAVECAR` file exists, VASP will read it; if it doesn't, VASP won’t read it. Understanding the **ISTART** parameter is important:

- If a `WAVECAR` file was saved in a previous calculation, and ISTART is not set, VASP will read it by default.
- If there is no `WAVECAR` file, even if **ISTART=1** or **ISTART=2** is set, VASP will not find the file and will continue the calculation without error.

### Controlling `WAVECAR` Output:

1. Control `WAVECAR` output using the **LWAVE** parameter.
2. Reading `WAVECAR` can significantly reduce self-consistent calculation time. However, `WAVECAR` files are often very large (hundreds of megabytes or even several gigabytes).
3. If a `WAVECAR` file was saved in previous steps (during optimization), it can be read during the DOS calculation (whether it's a one-step or two-step process), which will speed up the calculations.

## Summary:

After structure optimization:

- **One-step DOS calculation** requires the following parameters:

  1. **ISMEAR = -5**
  2. **LORBIT = 11**
  3. High-density K-point mesh

- **Two-step DOS calculation** requires the following parameters:
  1. **Step 1**:
     - **ISMEAR = -5**
     - **LCHARG = .TRUE.**
     - A slightly lower density K-point mesh
  2. **Step 2**:
     - **ISMEAR = -5**
     - **ICHARG = 11**
     - **LORBIT = 11**
     - High-density K-point mesh

If the `WAVECAR` file is saved during structure optimization, it can be read during the DOS calculation, allowing for a one-step calculation.

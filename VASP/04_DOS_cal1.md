# VASP Calculation: DOS and Band Structure Calculation for Materials

## Basics: How to Fit or Optimize a Stable Crystal Cell Structure and Calculate the Corresponding DOS Information

### 1. KPOINTS

#### 1.1 Number of K-Points

Compared to structural optimization, more K-points are required when calculating DOS because the more K-points used, the higher the quality of the DOS plot. As stated on the official VASP website:

> A high-quality DOS usually requires very fine K-meshes.

#### 1.2 Choosing the Number of K-Points

The more K-points, the better. But how should we choose the number of K-points?

Remember the empirical rule for selecting K-points that we previously discussed? That rule can be considered the standard for selecting K-points in normal calculations. For DOS calculations, the configuration needs to be upgraded.

Generally, the following formula can be used as a guideline for choosing K-points:

> \(K \times a = 45\)

This should be sufficient in most cases, and you can use this rule to select your K-points.

### 2. NEDOS

The NEDOS parameter also plays a crucial role in the quality of the DOS plot. For example, if our DOS energy range (the x-axis of the DOS plot) is: [-10 eV, 10 eV], VASP will, by default, divide this energy range into 301 points and then plot the graph. The 301 points correspond to the default NEDOS value.

If the NEDOS value is large enough, the DOS range will be more finely divided, improving accuracy.

Typically, a value of:

> NEDOS = 3000

is sufficient. Increasing the value too much won't be meaningful, but remember:

- The larger the NEDOS, the larger the DOSCAR and vasprun.xml files will be, taking up more storage space.
- If the DOS plot has many sharp peaks, you can try increasing NEDOS to smooth them out.

For more information, refer to [VASP Wiki on NEDOS](https://cms.mpi.univie.ac.at/wiki/index.php/NEDOS).

### 3. ISMEAR (Part 1)

According to the official VASP Wiki: [ISMEAR VASP Documentation](https://cms.mpi.univie.ac.at/wiki/index.php/ISMEAR)

> For the calculation of the total energy in bulk materials, we recommend the tetrahedron method with Blöchl corrections (ISMEAR=-5). This method also provides a smooth, nice electronic density of states (DOS).

This means that when ISMEAR = -5 (Blöchl correction tetrahedron method) is used, we can obtain a very smooth DOS plot.

#### Notes:

##### 3.1 Number of K-Points

When setting ISMEAR = -5, if the number of K-points is less than or equal to 4, an error will occur during the calculation, and the following error message will be displayed:

> VERY BAD NEWS! internal error in subroutine IBZKPT: Tetrahedron method fails for NKPT < 4. NKPT = 1

This is a common error. The official documentation states that the tetrahedron method is not applicable if less than 3 K-points are used:

> "The tetrahedron method is not applicable if fewer than three k-points are used."

If there are not enough K-points, using ISMEAR = -5 will result in an error. The solution is simple: increase the number of K-points and then use ISMEAR = -5 (this is a direct but recommended solution). However, if increasing the K-points still causes errors, you might encounter this:

> WARNING: DENTET: can't reach specified precision
> Number of Electrons is NELECT = ...

This is due to insufficient precision in the Fermi level with the tetrahedron method, causing an inconsistency between the number of electrons integrated to the Fermi surface and the number of valence electrons in the system. For more details, see [VASP Forum on DENTET Error](http://cms.mpi.univie.ac.at/vasp-forum/viewtopic.php?t=416).

##### 3.2 Suitable Systems:

ISMEAR = -5 is suitable for DOS calculations in all systems. However, it should not be used for structural optimization in metallic systems because the tetrahedron method cannot properly handle electron occupancy at the Fermi level, leading to a percentage error in the calculated forces. For metal structural optimizations, use ISMEAR ≥ 0.

For semiconductors and insulators, ISMEAR > 0 is not allowed, so ISMEAR <= 0 must be used.

Additionally, ISMEAR = -5 is not suitable for band structure calculations for any system.

When using ISMEAR = -5, the value of SIGMA does not affect the results. But if you're unsure, set SIGMA = 0.01.

##### 3.3 Summary:

For DOS calculations, as long as the number of K-points is no less than 3, using ISMEAR = -5 is a good choice.

### 4. ISMEAR (Part 2)

If the system is large and you can only use gamma points for the calculation, ISMEAR = -5 will surely fail. But if the server's resources are insufficient and you can't increase the number of K-points, what can you do?

For any system (even with fewer than 4 K-points), you can use:

> ISMEAR = 0; SIGMA = 0.01

This will yield ideal results in most cases. When using the Gaussian smearing method (ISMEAR=0), the value of SIGMA should be tested to ensure the entropy term \(T \times S\), when averaged per atom, is less than 0.001 eV (1 meV). If you don’t want to test, use a small value like SIGMA = 0.01.

For metallic systems, you can also use:

> ISMEAR = 1; SIGMA = 0.01.

If SIGMA is too large, the calculated energy may be incorrect. The smaller SIGMA is, the more accurate the calculation, but the longer it takes. A SIGMA value of 0.01 is already quite small, so there's no need to make it even smaller. For metallic systems, using the MP method (ISMEAR=1..N), SIGMA = 0.10 should suffice. The official recommendation is SIGMA = 0.20.

For more information, refer to the [VASP Guide on Smearing Methods](http://cms.mpi.univie.ac.at/vasp/guide/node124.html).

After determining KPOINTS, it's a good idea to test different values of SIGMA. The principle is as follows: the value of SIGMA should be as large as possible while keeping the average entropy energy **T\*S** per atom below **1 meV**. This ensures both accuracy and faster convergence.

### 5.Summary

- Settings for **KPOINTS** and **NEDOS**.
- Use **ISMEAR = -5** for DOS calculations.
- If the number of K points is too small, or if there are hardware limitations, use:
  - **ISMEAR = 0** for all systems.
  - For metallic systems, **ISMEAR = 1..N**, and the recommended SIGMA is **0.20**. However, values between **0.01 and 0.10** are also safe.
- Test the value of **SIGMA** to ensure the entropy energy **T\*S** per atom is less than **1 meV**.
- For non-DOS calculations:
  - For metallic systems, do not use **ISMEAR = -5**. Use **ISMEAR = 1** instead.
  - For semiconductors and insulators, **ISMEAR must not be greater than 0**.
  - **ISMEAR = 0** is a safe choice for all systems, but the value of **SIGMA** should be tested.

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

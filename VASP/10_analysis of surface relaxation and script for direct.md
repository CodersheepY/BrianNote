## Analysis of Surface Relaxation Calculation Results and Direct Coordinate Conversion Script

For most metallic systems, surface relaxation generally occurs in the vertical direction. When a metal is split along the surface direction, the coordination number of surface atoms changes.

We can imagine: the binding ability between the surface atoms and the layer below will be stronger, thereby lowering the energy of the system. This is reflected in two aspects: the interlayer distance and energy. Example: Exercise 44 on the Cu(111) surface.

### 1. Changes Before and After Structural Optimization

First, observe the changes in the structure. The structure of the POSCAR file is as follows (the freshly cut surface):

```
Cu\(1\1\1)
1.00000000000000
2.5717000960999998 0.0000000000000000 0.0000000000000000
-1.2858500480999999 2.2271576141999998 0.0000000000000000
0.0000000000000000 0.0000000000000000 21.2994003296000010
Cu
4
Selective
Cartesian
+0.0000000000 +0.0000000000 +0.0000000000 F F F
-0.0000128765 +1.4847792201 +2.0999078998 F F F
+1.2858629618 +0.7423784587 +4.1996028482 T T T
+0.0000000000 +0.0000000000 +6.2995107693 T T T
```

The structure after optimization (cooled down and converted to Cartesian coordinates from the CONTCAR file):

```
Cu\(1\1\1)
1.00000000000000
2.5717000960999998 0.0000000000000000 0.0000000000000000
-1.2858500480999999 2.2271576141999998 0.0000000000000000
0.0000000000000000 0.0000000000000000 21.2994003296000010
Cu
4
Selective dynamics
Cartesian
+0.0000000000 +0.0000000000 +0.0000000000 F F F
-0.0000128765 +1.4847792201 +2.0999078998 F F F
+1.2858542192 +0.7423835062 +4.1890255872 T T T
+0.0000182442 -0.0000105333 +6.2592676588 T T T
```

The z-coordinate of the surface atom decreased from 6.2995 Å to 6.2593 Å, indicating that the surface atom contracted toward the bulk.

The interlayer distance between the first and second layers (in Å):

```
6.2995 - 4.1996 = 2.0999 (POSCAR)
6.2593 - 4.1890 = 2.0703 (CONTCAR)
```

The change is:

\[
\frac{(2.0999 - 2.0703)}{2.0999} \times 100\% = 1.4\%
\]

**Note**: Here, we did not distinguish between positive and negative signs. Generally, a negative value indicates contraction towards the bulk.

### 2. Relaxation Energy

How do we know the energy change before and after relaxation?

Step 1: Perform a single-point calculation on the freshly cut surface and obtain an energy.

Step 2: Subtract the energy from the optimized structure from the single-point calculation result.

In practice, the first step can be skipped because VASP calculates the initial structure during optimization, which corresponds to the energy of the first ionic step.

```bash
grep ' without' OUTCAR | awk '{print $7}'
```

Results:

```
-13.96892338
-13.96932426
-13.96974724
-13.97030280
-13.97016827
-13.97041579
-13.97047584
-13.97063103
-13.97082922
-13.97086603
```

Thus, the energy change (in eV) before and after relaxation is:

\[
-13.97086603 - (-13.96892338) = -0.00194265
\]

The energy is negative, indicating that the relaxation process is exothermic. This suggests that the freshly cut surface is unstable, and after the surface atom contracts toward the bulk, the system's energy decreases, making it more stable.

**Note**: When performing this calculation, always check if the electronic steps of the first ionic step have converged. For example, VASP's default setting for the number of electronic steps in a single ionic step is 60 (NELM=60). If the default value is used and the 60 steps have not converged, the energy from this step cannot be used.

In such cases, you need to:

Increase `NELM` to 100 or a larger number and redo the single-point calculation for the unoptimized structure. If convergence is difficult, try adjusting the `ALGO` parameter.

---

### 3. `dire2cart.py`

The script `dire2cart.py` converts Direct coordinates to Cartesian coordinates. I write in 10_dire2cart.py

The conversion of the CONTCAR file above was done using this script. Usage example:

```bash
ls
CONTCAR INCAR KPOINTS OSZICAR OUTCAR POSCAR POTCAR sub4 vasprun.xml
dire2cart.py POSCAR

###################################
# for VASP 5.x or higher versions #
###################################

## This POSCAR has Direct Coordinations, Conversion is starting....
## POSCAR with Cartesian Coordinates is named as POSCAR_C
```

**Note**:

If the current POSCAR or CONTCAR is already in Cartesian coordinates, the script will automatically detect it and terminate the conversion. The converted file will be named `AAA_C`, where `AAA` is the file you want to convert (POSCAR or CONTCAR).

---

### 4. Extended Exercise

Download the PDF documentation on surface-related calculations from the official VASP website: _VASP Tutorial: A bit of surface science_.

**Supplement**:
In materials simulation, especially when using software like VASP, atomic positions in a crystal structure can be represented in two main coordinate systems: **Direct (fractional coordinates)** and **Cartesian coordinates**. Converting Direct coordinates to Cartesian coordinates is a common step that translates the relative positions of atoms into actual physical positions. This step is significant for the following reasons:

#### 1. **Direct Coordinates** (Fractional Coordinates)

Direct coordinates represent the position of atoms relative to the size of the unit cell (i.e., lattice constants). They are expressed as fractions of the lattice vectors, for example:

\[
\mathbf{r}\_{\text{direct}} = (x, y, z)
\]

where \(x\), \(y\), and \(z\) are numbers between 0 and 1, indicating the relative position of the atom within the unit cell.

#### 2. **Cartesian Coordinates**

Cartesian coordinates represent the actual physical positions of atoms in three-dimensional space, typically in Ångströms or nanometers. In crystallography, Cartesian coordinates are derived from the linear combination of lattice vectors, expressed as:

\[
\mathbf{r}\_{\text{cartesian}} = x \mathbf{a}\_1 + y \mathbf{a}\_2 + z \mathbf{a}\_3
\]

where \(\mathbf{a}\_1\), \(\mathbf{a}\_2\), and \(\mathbf{a}\_3\) are the basis vectors of the unit cell, and \(x\), \(y\), and \(z\) are the fractional coordinates in the Direct system.

#### 3. **Importance of the Conversion**

In material simulations, converting Direct coordinates to Cartesian coordinates serves the following purposes:

- **Physical Location**: Direct coordinates only represent the relative position within the unit cell and do not reflect the actual physical location of atoms. Cartesian coordinates provide the actual position in three-dimensional space, which is essential for combining with other physical quantities (e.g., potential energy, forces).

- **Geometric Analysis**: Some calculations or visualization tools require the actual absolute positions of atoms, thus necessitating the conversion to Cartesian coordinates. For example, in VASP, quantities like forces and velocities are computed based on Cartesian coordinates.

- **Simulation Efficiency**: In molecular dynamics, energy optimization, or other simulations, Cartesian coordinates are often required to compute geometrical information like interatomic distances and angles because these quantities are more intuitive in the Cartesian system.

#### 4. **Conversion Formula**

The formula for converting Direct coordinates to Cartesian coordinates is:

\[
\mathbf{r}_{\text{cartesian}} = \mathbf{r}_{\text{direct}} \cdot \mathbf{A}
\]

where \(\mathbf{A}\) is the matrix of lattice vectors:

\[
\mathbf{A} = \begin{bmatrix}
a_1 & a_2 & a_3 \\
b_1 & b_2 & b_3 \\
c_1 & c_2 & c_3
\end{bmatrix}
\]

This means multiplying the Direct coordinates \((x, y, z)\) by the lattice vector matrix \(\mathbf{A}\) gives the Cartesian coordinates \((x*{\text{cart}}, y*{\text{cart}}, z\_{\text{cart}})\). This formula ensures that the relative positions in the crystal lattice are transformed into their absolute positions in space.

For example, given the lattice vectors:

\[
\mathbf{A} = \begin{bmatrix}
a_1 & a_2 & a_3 \\
b_1 & b_2 & b_3 \\
c_1 & c_2 & c_3
\end{bmatrix}
\]

The conversion for each atom's position in Direct coordinates \((x, y, z)\) is calculated as:

\[
\mathbf{r}\_{\text{cartesian}} = x \mathbf{a}\_1 + y \mathbf{a}\_2 + z \mathbf{a}\_3
\]

Thus, this conversion allows researchers to move between relative and absolute positions in the unit cell for calculations that require Cartesian coordinates.

---

### 5. **Conclusion**

In summary, converting Direct coordinates to Cartesian coordinates is crucial in materials simulation, especially when analyzing physical quantities such as energy, forces, and interatomic distances. Direct coordinates are suitable for describing positions within a unit cell, while Cartesian coordinates are necessary for calculating physical properties in three-dimensional space.

By utilizing conversion scripts like `dire2cart.py`, the process becomes automated, allowing for smooth transitions between these two coordinate systems during your simulations.

Other resources related to surface science and VASP calculations, including tutorials, can be found on the [NERSC VASP Workshop](http://www.nersc.gov/users/training/events/3-day-vasp-workshop/) page.

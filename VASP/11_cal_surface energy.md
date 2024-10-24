## Surface Energy Calculation

### Definition of Surface Energy

Surface energy is defined as the energy consumed in breaking chemical bonds between molecules when creating a surface. In solid-state physics theory, surface atoms possess more energy compared to atoms within the bulk material. Therefore, according to the principle of minimum energy, atoms will spontaneously move towards the bulk rather than the surface. Another definition of surface energy is the extra energy of the surface relative to the bulk material. Breaking the internal chemical bonds of a solid material to divide it into smaller pieces requires energy. If this process is reversible, the energy required to divide the material into smaller pieces is equal to the energy added to the surface of the smaller pieces. However, in reality, only surfaces freshly formed in a vacuum adhere to this energy conservation principle. This is because newly formed surfaces are highly unstable and reduce surface energy through surface atom rearrangement, interaction with other atoms or molecules, or adsorption from the surroundings.

---

### 2. VASP Calculation

Notes from learning VASP’s official surface energy calculation process.

First, browse the formula and analyze the meaning of each term:

- **Erel**: Relaxation energy, the energy released when the surface relaxes to stability after being freshly cut (as mentioned in section 09). The first point in the diagram represents the energy of the freshly cut surface (Esurf), and the last point represents the energy after optimization is complete.
- **σ**: In _Hand-on Session III_, this denotes surface energy.
- **Esurf**: The energy of the freshly cut slab, which can either be calculated through a single-point calculation or by taking the energy from the first ionic step during optimization.
- **Natoms**: The number of atoms in the slab.
- **Ebulk**: The energy per atom in the bulk structure.

---

### Notes:

**Note 1:**

Regarding the final point above, **Ebulk** does not refer to the total energy of the bulk but to the energy per atom in the bulk, i.e., the bulk energy divided by the number of atoms. Let's review the previous optimization of the unit cell:

1. **Obtain the stable lattice constant** using two methods:
   - ISIF + Large ENCUT: By calculating the energy of the unit cell with different lattice constants, the Birch-Murnaghan equation is used for fitting.
2. After obtaining the stable lattice constant, if ISIF=3 is used, copy the CONTCAR file to POSCAR, and then set ISIF=2 in the INCAR file along with normal ENCUT to perform a single-point calculation to determine the unit cell energy.
3. Divide the calculated unit cell energy by the number of atoms in the unit cell to obtain **Ebulk**. For example, if the Cu unit cell contains 4 atoms, the total energy of the system after calculation is divided by 4.

---

**Note 2:**

There is some doubt about the surface energy presented in _Hand-on Session III_ because surface energy, by definition, should have units of energy per area. However, the diagram only shows energy without specifying the area. Therefore, the "surface energy" in the diagram refers to pure energy.

1. **σunrel**: This represents the work required to directly cut the bulk into two pieces. Since the bulk is split into two surfaces, the formula divides by 2. This calculation ignores surface relaxation and focuses solely on the process of splitting the bulk into two surfaces. The energy is the freshly cut slab's energy: -25.560 eV.
2. **Total energy change** from bulk to the cut surface until reaching a stable state.

---

**Note 3:**

**Erel** is not divided by 2 because only the top two layers of atoms are relaxed, i.e., only one surface is optimized.

Thus, the energy change from the bulk to the stable surface has been fully calculated.

As mentioned earlier, surface energy includes area, so how do we calculate the area?

The area is the surface area of the slab being calculated.

Now, let's look at the slab model. If the area is represented by **A**, then:

\[
A = 2.57170 \times 2.22716 - 0.000 \times (-1.28585)
\]

---

**Note 4:**

The area is calculated as:

\[
A = x_1 \times y_2 - y_1 \times x_2
\]

In _Hand-on Session III_, the area of the Ni(100) surface is:

\[
A = (0.5 \times 0.5 - 0.5 \times (-0.5)) = 0.5 \, \text{Å}^2
\]

We also need to multiply by the scaling factor, so:

\[
A = 0.5 \times 3.53 \times 3.53 = 6.23045 \, \text{Å}^2
\]

Finally, the surface energy is:

\[
\sigma = \frac{0.71 \, \text{eV}}{6.23045 \, \text{Å}^2} = 0.114 \, \text{eV/Å}^2
\]

In general, the unit of surface energy is **J/m²**.

The conversion factor is:

\[
1 \, \text{eV/Å}^2 = 16.02 \, \text{J/m}^2
\]

After conversion, the surface energy is:

\[
\sigma = 0.114 \times 16.02 \, \text{J/m}^2 = 1.826 \, \text{J/m}^2
\]

---

**Note 5:**

When obtaining **Esurf**, it is crucial to check the following:

1. Has the single-point calculation converged?
2. Has the first ionic step converged?

## Surface Energy Calculation in Hand-On-Session-III

### 1. Explanation of Dividing by 2 in the Calculation

When performing a single-point calculation, there are two surfaces included. This is because the bulk is split into two halves. Let’s first examine the slab model:

**Problem**: There are two surfaces in the structure shown in the figure because our slab model is periodic in the z-direction. Even though we fix the bottom two layers of atoms to simulate the bulk, they still count as surfaces. See the illustration below:

So, why isn't **Erel** divided by 2?

This is because only one surface is optimized.

If both surfaces need to be optimized, **Erel** should be divided by 2.

If we want to optimize both the top and bottom surfaces, this leads to two types of slab models: symmetric and asymmetric models. The previous models were all asymmetric slab models. Symmetric slab models are as shown below:

In this model, the three middle layers of atoms are fixed to simulate the bulk, while the top and bottom three layers of atoms are relaxed to simulate the surfaces. In this case, when calculating surface energy, the value of **Erel** should be divided by 2. The simplified formula for the calculation is:

\[
\sigma = \frac{E*{\text{rel_surf}} - N*{\text{atoms}}E\_{\text{bulk}}}{2A}
\]

Where:

- **Erel_surf** is the energy of the optimized symmetric slab model.
- **A** is the surface area of a single surface.
- **Natoms** is the number of atoms in the slab model.
- **Ebulk** is the bulk energy per atom.

**Note**: Compare this formula with the one at the bottom of page 96 in the reference book, and carefully read the explanation of the terms in the formula on page 97.

---

### 2. Differences Between Symmetric and Asymmetric Slab Models

1. **Symmetric slab** models have more atoms and are longer in the z-direction, requiring more vacuum layers (top and bottom), which inevitably leads to larger computational costs. Especially for those with limited computational resources, this model may not be ideal.

2. **Asymmetric models**, due to their asymmetry, will introduce a dipole moment during the calculation. This dipole moment is caused by the model itself and is not the same as the molecular dipole moment we commonly refer to. Particularly when calculating surface adsorption or reactions, the interaction between the dipole moments of the slabs can affect the results. However, this is not a major issue since most computational software can control the calculation parameters to minimize or eliminate this effect.

3. **In VASP**, this issue can be addressed by setting:
   - `LDIPOL = .TRUE.` (to enable dipole correction),
   - `IDIPOL = 1, 2, or 3`, where 1, 2, and 3 correspond to corrections along the x, y, and z directions, respectively. Setting `IDIPOL = 4` applies corrections in all directions.

---

### 3. Parameters Influencing Surface Energy Calculation

From the formula, it is clear that the energy of the slab and the energy of the bulk are the main factors:

1. **Slab Energy**:

   - **Number of layers in the slab**: Refer to the results in the reference book.
   - **Surface size of the slab**: Generally, a p(1x1) surface is sufficient. However, you can also compare the difference between p(1x1) and p(2x2). It is important to note that when changing the surface size, the **KPOINTS** should also be adjusted accordingly to ensure the results are comparable.
   - **Vacuum layer thickness**: Consider testing different thicknesses.

2. **Bulk Energy**: We’ve already discussed how to calculate bulk energy in the previous section. Some students asked: How many atoms should be included in the bulk model for the calculation?

   - **Personal opinion**: The calculation can be done using either a unit cell or a primitive cell. Of course, you can also expand the unit cell. However, when comparing the energy per atom in models of different sizes, make sure to pay attention to the selection of **KPOINTS**. If the **KPOINT** density is inconsistent, even small differences in the energy per atom may occur, which can be confusing. You can test this if you're interested.

---

### 4. Reference

For further reading on surface energy, see _Theoretical Surface Science: A Microscopic Perspective_ by Axel Groß. The book introduces surface energy concepts on page 68. The author’s group homepage can be found at:

[Institute of Theoretical Chemistry, University of Ulm](https://www.uni-ulm.de/en/nawi/institute-of-theoretical-chemistry/).

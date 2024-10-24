## Slab Model

### Mechanism Research on Surface Reactions

The first step is to construct a model.

1. Calculate the structure of the bulk and obtain stable lattice constants.
2. Based on this, construct the slab model.

### 1. What is the Slab Model?

Reference: _Density Functional Theory: A Practical Introduction_, Chapter 4: DFT Calculations for Surfaces of Solids.

In heterogeneous catalysis, reactions occur on the surface of catalysts, which is at the interface between the gas (or liquid) and solid phases. To simulate the surface of a solid, a structure with periodicity in two dimensions is needed to mimic the surface of the solid, while the third dimension is non-periodic to simulate the gas or liquid phase. However, it is more common to add periodicity in the third dimension as well, with a free space (vacuum layer) between the periodic structures to separate them. It is important to ensure that there is no interaction between the two periodic structures in the third dimension, meaning the electron density in one periodic unit along the vacuum layer direction should approach zero, avoiding any influence on the neighboring structures. The steps to achieve this are as follows:

1. **Increase the vacuum layer thickness**. The thickness must be appropriate—too small is obviously not acceptable, but too large is also unsuitable. For example, one of the factors affecting calculation speed is the thickness of the vacuum layer. A thicker vacuum layer implies a larger model size, which in turn slows down the calculation. Generally, for surface reaction calculations, a vacuum thickness of 15 Å is sufficient. However, for calculations sensitive to the vacuum layer, such as work function, special care is needed.

2. **In VASP**, the following parameters need to be included: `LDIPOL = .TRUE.` and `IDIPOL = 3` (indicating the z-direction) to eliminate the dipole moment caused by the asymmetric slab surface. [More details](https://cms.mpi.univie.ac.at/wiki/index.php/LDIPOL).

For surface structures, the following points need attention:

1. **Size of the surface in the xy direction**:

   - A) This affects the coverage of adsorbed species on the surface.
   - B) It affects the system size and calculation time.
   - C) Different sizes require selecting corresponding k-points.

2. **Different crystal planes, such as (111), (100), (110)**:

   - A) It depends on the research direction.
   - B) Different crystal planes have different surface energies.
   - C) Different crystal planes also exhibit different surface structures and reactivity.

3. **Number of layers in the surface structure**:

   - A) More layers mean more atoms, increasing the system size in the z-direction and slowing down the calculation.
   - B) The number of layers that need relaxation during calculations.
   - C) The number of layers can affect properties such as surface energy.
   - D) Different crystal planes require different numbers of layers. Generally, more open surfaces need more layers.
   - E) Choose the appropriate number of layers based on the property being studied.

4. There are two types of slab models: symmetric (top and bottom surfaces are identical) and asymmetric. Symmetric structures usually require more layers, leading to a larger system. Asymmetric structures are smaller but can be affected by dipole moments. In such cases, it is necessary to add `LDIPOL` and `IDIPOL` to remove the dipole moment.

---

### 2. Building the Surface Model

The primary tool for cutting surfaces is _Material Studio_. Example: Cu(111).

### 2.1 Bulk Calculation

1. Build the Cu unit cell model and export the POSCAR file for VASP calculation using VESTA. If unclear, review the previous content related to bulk calculations.

2. **Preparation of the INCAR file**:
   - A) Since Cu is a metal, set `ISMEAR=1` and `SIGMA=0.1`.
   - B) Set `ENCUT=700`. Since the Cu unit cell contains very few atoms, the energy cutoff can be set higher.
   - C) `EDIFF` and `EDIFFG` control the convergence of electronic steps and the precision of ion relaxation.
   - D) `PREC=High`.

```ini
SYSTEM = Cu-Bulk calculation RELAX
ISTART = 0
ICHARG = 2
ENCUT = 700
EDIFF = 1e-6
IBRION = 1
NSW = 200
EDIFFG = -0.01
ISIF = 3
POTIM = 0.1
ISMEAR = 1
SIGMA = 0.1
PREC = H
LREAL = Auto
ALGO = Fast
NCORE = 4
```

- E) Here, `LREAL` is set to Auto to unify with the subsequent surface energy calculation. When calculating, choose the `LREAL` value based on the number of atoms and subsequent calculations. Do not compare energies calculated with `LREAL=.FALSE.` and `LREAL=ON/.TRUE.`. See the note below.
- F) **Note**: For bulk material calculations, to compute other properties later, it is often required to use the same `ENCUT`, `ENAUG`, `PREC`, `LREAL`, and `ROPT` throughout the entire calculation. You can refer to the manual, section 8.3: "What kind of 'technical' errors do exist" for an overview.

3. **KPOINTS**:
   With `a=b=c=3.614715 Å`, based on previous experience, a `13×13×13` k-point grid (Gamma-centered) is set.

4. **Prepare the POTCAR file** for the calculation. Refer to the previous related script for Ex30.

5. The final calculated lattice constant for Cu is 3.637015 Å, while the experimental value is 3.61515 Å.

---

### 2.2 Notes on Structure Optimization

When optimizing the structure, the following points need attention (Source: 7.6.2 Accurate bulk relaxations with internal parameters).

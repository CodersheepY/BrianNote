## Adsorption Energy Calculation

I have completed the optimization of the adsorption of an oxygen atom (O) on the Cu(111) surface. Let’s summarize and count the energies we have obtained from the previous calculations:

1. The energy of the Cu(111) slab.
2. The optimized energy of the system with O adsorbed on the slab.
3. The energy of the O₂ molecule.
4. The energy of an oxygen atom in the gas phase.

With these energies, we can calculate the adsorption energy of oxygen (O) on the p(1x1) Cu(111) surface, which is the main content of this section.

---

### 1. What is Adsorption Energy?

Adsorption energy, as the term suggests, is the energy released when a molecule or atom is adsorbed onto a surface from the gas phase. When calculating the adsorption energy, there are two key points to keep in mind: the initial state and the final state.

1. **Final State**: The final state refers to the optimized configuration of the O atom adsorbed on the surface. We denote this as `slab+O`, and its energy is labeled as **E(slab+O)**.

2. **Initial State**: For calculating the adsorption energy of oxygen, the initial state includes both an oxygen molecule (O₂) in the gas phase and the clean slab surface. Since the possibility of O₂ dissociating into two O atoms and then adsorbing onto the surface in an actual reaction is extremely low, we will use the energy of the O₂ molecule from previous optimization calculations. The energies of the initial structures are denoted as **E(O₂)** and **E(slab)**.

   - Since the final state contains only one oxygen atom, we need to use half the energy of the O₂ molecule for the initial state: **E(O₂)/2**.

3. **Note**: In some literature, people also use the energy of a single O atom in the gas phase to calculate the adsorption energy. In this case, the physical meaning is that the O₂ molecule dissociates into O atoms, and then the O atom is adsorbed. The initial state would then involve dissociated O atoms, and the dissociation energy of the O₂ molecule would not be considered. Since O₂ dissociation is an endothermic reaction, ignoring the dissociation energy will lead to an overestimation of the adsorption energy of the O atom. The energy of a single O atom is labeled as **E(O)**.

4. Therefore, when calculating, it is essential to clearly state the formula you are using. Additionally, when comparing adsorption energies calculated with different methods, be careful not to compare results that use O₂ as the reference with those that use O atoms as the reference.

5. Which of these two methods is correct? **Answer**: Both are correct. Since the dissociation energy of O₂ is a constant, whether you include it (by using O₂ energy as the starting point) or ignore it (by using O atom energy as the starting point), the difference between the results is merely this constant. Refer to the calculation below.

---

### 2. Calculation of Adsorption Energy (Single O Atom)

#### Formula 1:

$E_{\text{ads}} = E_{\text{slab+O}} - \left(E_{\text{slab}} + \frac{E(O_2)}{2}\right)$

#### Formula 2:

$E_{\text{ads}} = E_{\text{slab+O}} - \left(E_{\text{slab}} + E(O)\right)$

Where:

- **E(slab+O)**: The energy of the optimized structure with O adsorbed on the slab.
- **E(slab)**: The energy of the optimized slab.
- **E(O₂)**: The energy of an O₂ molecule in the gas phase.
- **E(O)**: The energy of a single O atom in the gas phase.

After substituting the data:

---

If I use the energy of a single O atom in the gas phase, the dissociation energy of the O₂ molecule (or the binding energy of the O atom) is ignored. The dissociation of the O₂ molecule is an endothermic reaction, and the dissociation energy can be obtained in two ways:

1. Subtracting two adsorption energies.
2. Directly calculating the dissociation energy from the energy of the O₂ molecule and the energy of an O atom.

The binding energy of a single O atom is:

$E_{\text{bind}} = \frac{E(O_2)}{2} - E(O)$

After substituting the values, we obtain:

---

### 3. Analysis:

The calculated adsorption energy \( E\_{\text{ads}}(O) = 1.216 \, \text{eV} \) is a positive value. Why is this the case? Was the calculation wrong? Or is there another physical meaning?

---

### Explanation:

The positive adsorption energy means that the adsorption process is endothermic. In other words, energy is required to adsorb the O atom onto the surface, which can happen if the initial state (O₂ molecule or O atom in the gas phase) is more stable than the final state (O atom adsorbed on the surface). This could be due to factors such as the weak interaction between O and the Cu surface, or possibly because the energy of the adsorbed state is not sufficiently lower than that of the gas-phase O atom or molecule.

This result doesn't necessarily imply an error in the calculation; it might indicate an inherent characteristic of the system being studied.
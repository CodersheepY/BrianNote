Adsorption Energy Calculation

    1.	Definition of Adsorption Energy

The adsorption energy is the energy difference between the initial and final states:
• Initial state: Pure slab and an O₂ molecule or O atom.
• Final state: Slab with adsorbed O atom. 2. Optimization of O₂ Molecule Calculation
When optimizing the O₂ molecule calculation, keep in mind:
• Obtain the bond length of the O₂ molecule from the database as the initial value.
• Set the following parameters:
• ISPIN = 2
• ISMEAR = 0; SIGMA = 0.05
• IBRION = 2; POTIM = 0.1
• Use only the Gamma point for calculations. 3. Optimization of O Atom Calculation
The calculation details are similar to those for the O₂ molecule, but note:
• Since the O atom remains a single O atom regardless of optimization, a single-point calculation is more appropriate.
• The box size should be 13x14x15 Å³ or 13.1x13.2x13.3 Å³. Avoid cubic boxes like 13x13x13. 4. Bulk Structure Optimization
To optimize the bulk structure, perform a lattice sweep and fitting or optimize directly with ISIF = 3. Set a high ENCUT, ideally around 700 eV. 5. Slab Optimization
Use selective dynamics and ensure the correct relationship between the T and F flags in the coordinates. 6. Setting Up the O Adsorption Model
To create a reasonable initial model for O adsorption:
• Use experimental values if available.
• If not, perform a low-precision initial calculation to estimate values. 7. Identifying Calculation Completion
Learn to determine when calculations have finished. Different task types will have unique completion signatures.

Why is the Adsorption Energy of O on the p-(1x1)-Cu(111) Surface Positive?

Example Analogy:

Imagine a subway car with empty seats:

    •	Empty car: You can pick any seat you want.
    •	Some people are seated: Choices are fewer, but you still have options.
    •	Crowded car: Even getting out is a challenge; finding a seat is out of the question.

Returning to the adsorption scenario, imagine replicating the O adsorption structure on the xy-plane. The surface is densely packed with O atoms! In such a cramped space, if you were an O atom, would you want to squeeze in? Would you choose an empty car or a packed one?

As an O atom, I’d definitely prefer the empty car. So, if adsorption must occur on this crowded surface, some external force is needed — like a strong JR station attendant pushing passengers onto the train. This push is analogous to the adsorption energy for the O atom. A positive adsorption energy indicates that additional force is required to adsorb it onto the surface. The more crowded the surface, the larger the positive value.

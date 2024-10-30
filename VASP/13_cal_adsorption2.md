Adsorption Energy Calculation

1. Definition of Adsorption Energy

The adsorption energy is calculated as the energy difference between the initial and final states:

    •	Initial state: A clean slab and an O₂ molecule or an O atom.
    •	Final state: A slab with adsorbed O atom.

2. Optimization of O₂ Molecule

When optimizing the O₂ molecule, consider the following:

    •	Bond Length: Retrieve the O₂ bond length from the database to use as the initial value.
    •	Calculation Settings:
    •	ISPIN = 2
    •	ISMEAR = 0 and SIGMA = 0.05
    •	IBRION = 2 and POTIM = 0.1
    •	k-point: Use only the Gamma point.

3. Optimization of O Atom

For the O atom calculation, follow the same details as above with a few differences:

    •	Single-Point Calculation: Since an O atom remains a single atom after optimization, perform a single-point calculation rather than a full optimization.
    •	Box Size: Use box dimensions of either 13x14x15 Å³ or 13.1x13.2x13.3 Å³. Avoid cubic boxes like 13x13x13.

4. Bulk Structure Optimization

To optimize the bulk structure, use:

    •	Methods: Lattice sweep and fitting or direct optimization with ISIF = 3.
    •	High Cutoff Energy: Set ENCUT to a high value, around 700 eV.

5. Slab Optimization

For slab optimization:

    •	Selective Dynamics: Use selective dynamics and ensure correct T and F flags for the coordinates.

6. Setting Up the O Adsorption Model

To construct a reasonable initial model for O adsorption:

    •	Experimental Values: Use experimental data for initial parameters if available.
    •	Estimation: If no experimental data is available, perform a low-precision preliminary calculation to estimate values.

7. Recognizing Calculation Completion

Learn to recognize the completion of different types of calculations. Each task type has unique characteristics indicating its completion.

Why is the Adsorption Energy of O on the p-(1x1)-Cu(111) Surface Positive?

Analogy

Imagine entering a subway car:

    1.	Empty Car: You have your pick of seats.
    2.	Partially Filled Car: Your options are more limited, but you can still find a seat.
    3.	Crowded Car: Even getting out becomes challenging, let alone finding a seat.

Applying the Analogy

Now, consider the adsorption structure for O atoms replicated in the xy-plane:

    •	High Surface Density of O Atoms: The surface is tightly packed with O atoms. In this crowded environment, if you were an O atom, would you want to squeeze in? You’d likely prefer the “empty car” over the crowded one.

As an O atom, I wouldn’t want to adsorb onto such a crowded surface. However, if adsorption must occur, an external force is required to make it happen — similar to a strong JR station attendant pushing passengers onto a packed train. This push represents the adsorption energy for the O atom.

Interpretation of Positive Adsorption Energy

A positive adsorption energy means additional force is required for the O atom to adsorb onto the surface. The more crowded the surface, the higher the positive value of the adsorption energy.

Factors Influencing Adsorption Energy on Surface Coverage

Surface coverage affects adsorption energy due to several factors, including:

    •	Electron Distribution Shifts in the Bulk: Redistribution of electrons within the bulk material can influence adsorption strength.
    •	Bond Formation with Adsorbate Species: The nature of bonding between the adsorbate and the substrate.
    •	Repulsion and Attraction Among Adsorbate Species: Interactions between adsorbed atoms or molecules can either enhance or weaken adsorption based on their distance and mutual orientation.

For more detailed insights, refer to this ScienceDirect article.
https://www.sciencedirect.com/science/article/abs/pii/S0360056402450134

Here’s the summary in markdown format:

Summary

When surface coverage is high, the adsorption energy is small or even positive, indicating that species on the surface are unstable and need extra force to stay in place.

The surface area is limited, and when O atoms adsorb on the p-(1x1) surface, it’s too crowded for them to be comfortable. However, when the surface expands to p-(2x2), there’s noticeably more space. Although the O atoms still aren’t completely happy, they are significantly less uncomfortable compared to being on the (1x1) surface. This is evident from the adsorption energy, which decreases by about 1 eV.

From another perspective, this can be seen as a battle among O atoms for limited surface resources.

***一 What is Density Functional Theory***
1. **Density**: Electron Density
2. **Functional**: Energy as a function of Electron Density. By minimizing the energy functional, one can obtain the ground-state electron density and the corresponding energy of the system.
3. **Electron Density**: A function of spatial coordinates.
4. **Functional**: A function of a function.
5. **Density Functional Theory**: A method for studying the electronic structure of multi-electron systems (describing the properties and behavior of atoms and molecules) through electron density. The electronic motion equations are no longer the Schrödinger equation (whose wave function represents the solution to a many-body problem that is difficult to solve in practical situations), but are instead solved in terms of the electron density, simplifying the computational process.

***二 Hohenberg-Kohn Theorems***
1. **H-K First Theorem**: For a given electron density *\(n(r)\)*, there exists a unique corresponding ground-state external potential Vext(r)*\( V_{\text{ext}}(r) \)* and ground-state total energy *\( E \)*. Thus, the ground-state energy of the Schrödinger equation can be expressed in terms of the electron density.
<img width="239" alt="image" src="https://github.com/CodersheepY/BrianNote/assets/87512544/7eb08a3f-7883-4105-92f3-7576332c7070">

We can see that the degrees of freedom of the Schrödinger equation have decreased from 3N to 3, greatly reducing the complexity of computations.

Uniqueness proof:まだproof

2. **Hohenberg-Kohn Second Theorem**: The electron density that yields the lowest energy is the ground-state electron density. The second theorem can be regarded as a variational principle of the Hohenberg-Kohn theorem and can be derived through the variational principle of the Schrödinger equation.

***三  Density Functional Theory (DFT)***
According to the Hohenberg-Kohn theorem, the wave function<img width="14" alt="image" src="https://github.com/CodersheepY/BrianNote/assets/87512544/1ca3d589-fc6c-49e7-ac14-c05670cbe2cd">
, energy<img width="15" alt="image" src="https://github.com/CodersheepY/BrianNote/assets/87512544/4dc80cd7-2dc5-4786-8ef3-e7974c825d4a">
, and various other quantum mechanical properties are all functions of the electron density \( n \). The electron density \( n \) can be expressed as:
<img width="456" alt="image" src="https://github.com/CodersheepY/BrianNote/assets/87512544/ed7e40c7-b948-40cd-b088-44c14c9f40c4">

The total energy of the system can be expressed as:
<img width="395" alt="image" src="https://github.com/CodersheepY/BrianNote/assets/87512544/062a5b94-510b-4df5-8533-c3960aa3bb53">
The system now has only three degrees of freedom, which facilitates analysis of larger systems.

Next, the three components of total energy: kinetic energy, external potential energy, and interaction energy, are approximated and simplified.
For the kinetic energy term: Neglecting the electron-electron interactions and assuming that electrons exist in individual single-electron orbitals, the kinetic energy term can be written as:
<img width="375" alt="image" src="https://github.com/CodersheepY/BrianNote/assets/87512544/6986a360-115d-4025-af2f-bc5c6caa172d">
Where<img width="43" alt="image" src="https://github.com/CodersheepY/BrianNote/assets/87512544/8e095ad3-08f4-4cac-a8b1-1fda88c4c869">
represents the \( i \)th single-electron orbital, and \( T_s(n) \) is the sum of kinetic energies of all such non-interacting single-electron orbitals.

External potential energy term: According to the Born-Oppenheimer approximation<img width="175" alt="image" src="https://github.com/CodersheepY/BrianNote/assets/87512544/64223201-c02f-488c-8959-9fa5c11e90ec">

Interaction term: For electron-electron interactions, a Thomas-Fermi model can be employed to derive an intuitively appealing approximation (Note: The proposal of the Thomas-Fermi model itself was not intuitive).
<img width="290" alt="image" src="https://github.com/CodersheepY/BrianNote/assets/87512544/e68d498c-3aeb-40a6-bb0e-3ca485f58932">
It can be observed that the form of <img width="51" alt="image" src="https://github.com/CodersheepY/BrianNote/assets/87512544/edfef052-30b8-421b-b6da-1d5f43bd7969">
closely resembles the classical Coulomb interaction expression. 

Adding these three terms together yields an approximate energy value.
<img width="285" alt="image" src="https://github.com/CodersheepY/BrianNote/assets/87512544/79a0d6bc-f46d-4e41-a939-5f4ae28662d7">

The error introduced by this approximation is referred to as the exchange-correlation energy.
<img width="262" alt="image" src="https://github.com/CodersheepY/BrianNote/assets/87512544/25b85880-71bf-4b99-81c1-2686295d8e50">

The exact total energy of the system is referred to as the "exact total energy.
<img width="328" alt="image" src="https://github.com/CodersheepY/BrianNote/assets/87512544/49fa4c46-f855-4bf0-95ee-6ccd2c8789c2">

The above constitutes the framework of Density Functional Theory (DFT).

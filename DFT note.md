***一 What is Density Functional Theory***
1. **Density**: Electron Density
2. **Functional**: Energy as a function of Electron Density. By minimizing the energy functional, one can obtain the ground-state electron density and the corresponding energy of the system.
3. **Electron Density**: A function of spatial coordinates.
4. **Functional**: A function of a function.
5. **Density Functional Theory**
A method for studying the electronic structure of multi-electron systems (describing the properties and behavior of atoms and molecules) through electron density. The electronic motion equations are no longer the Schrödinger equation (whose wave function represents the solution to a many-body problem that is difficult to solve in practical situations), but are instead solved in terms of the electron density, simplifying the computational process.

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

1. **kinetic energy term**: Neglecting the electron-electron interactions and assuming that electrons exist in individual single-electron orbitals, the kinetic energy term can be written as:
<img width="375" alt="image" src="https://github.com/CodersheepY/BrianNote/assets/87512544/6986a360-115d-4025-af2f-bc5c6caa172d">
Where<img width="43" alt="image" src="https://github.com/CodersheepY/BrianNote/assets/87512544/8e095ad3-08f4-4cac-a8b1-1fda88c4c869">
represents the \( i \)th single-electron orbital, and \( T_s(n) \) is the sum of kinetic energies of all such non-interacting single-electron orbitals.

2. **External potential energy term**: According to the Born-Oppenheimer approximation<img width="175" alt="image" src="https://github.com/CodersheepY/BrianNote/assets/87512544/64223201-c02f-488c-8959-9fa5c11e90ec">

3. **Interaction term**: For electron-electron interactions, a Thomas-Fermi model can be employed to derive an intuitively appealing approximation (Note: The proposal of the Thomas-Fermi model itself was not intuitive).
<img width="290" alt="image" src="https://github.com/CodersheepY/BrianNote/assets/87512544/e68d498c-3aeb-40a6-bb0e-3ca485f58932">

It can be observed that the form of <img width="51" alt="image" src="https://github.com/CodersheepY/BrianNote/assets/87512544/edfef052-30b8-421b-b6da-1d5f43bd7969">
closely resembles the classical Coulomb interaction expression. 

4. Adding these three terms together yields an approximate energy value.
<img width="285" alt="image" src="https://github.com/CodersheepY/BrianNote/assets/87512544/79a0d6bc-f46d-4e41-a939-5f4ae28662d7">

The error introduced by this approximation is referred to as the exchange-correlation energy.
<img width="262" alt="image" src="https://github.com/CodersheepY/BrianNote/assets/87512544/25b85880-71bf-4b99-81c1-2686295d8e50">

The exact total energy of the system is referred to as the "exact total energy.
<img width="328" alt="image" src="https://github.com/CodersheepY/BrianNote/assets/87512544/49fa4c46-f855-4bf0-95ee-6ccd2c8789c2">

The above constitutes the framework of Density Functional Theory (DFT).

***四 Exchange-Correlation Functional***

**The core idea of DFT**: to derive an easily computable but inaccurate result through approximation and then isolate all errors into a separate term for further analysis. Below are various approximation methods used to calculate the exchange-correlation energy.

1. **Local Density Approximation, LDA**
LDA is the simplest exchange-correlation functional and was proposed quite early—almost concurrently with the inception of DFT.
Note: LDA itself is not the name of a functional but rather an approximation method; LDA encompasses many functionals.

**The idea of LDA**: it can be approximated that the exchange-correlation energy at each point only depends on the local electron density at that point. Therefore, the general expression for LDA can be written as:
<img width="167" alt="image" src="https://github.com/CodersheepY/BrianNote/assets/87512544/fc8cd972-0ad4-4647-890a-de1f14201e25">

The Exc part (x term) and the correlation part (c term) of Exc are discussed separately.
<img width="184" alt="image" src="https://github.com/CodersheepY/BrianNote/assets/87512544/93b626cf-4f56-4ed3-a7a0-3d72f4499504">

1）For the exchange part of LDA, the homogeneous/uniform electron gas (HEG/UEG) approximation is commonly used. This means that the exchange energy at a point in space for an atom/molecule is equivalent to that of a uniform electron cloud with the same electron density at that point. The exchange energy of UEG may be provided by the Dirac functional.

2）For the correlation part, different approximations are used for high-density and low-density UEG.
・Under high-density conditions, the correlation energy at a point in space can be expressed as:
<img width="209" alt="image" src="https://github.com/CodersheepY/BrianNote/assets/87512544/9137cc33-e512-4d67-aaa0-33a533ab9898">

In which <img width="16" alt="image" src="https://github.com/CodersheepY/BrianNote/assets/87512544/0e7900c4-a494-4f93-8265-fd7ff6cf161b">
 is the Wigner-Seitz radius, <img width="183" alt="image" src="https://github.com/CodersheepY/BrianNote/assets/87512544/bd9fbeb6-1aee-418c-a6f5-6073bdc7cc5a">
 ; <img width="35" alt="image" src="https://github.com/CodersheepY/BrianNote/assets/87512544/e76fbb92-cd27-4f2e-9883-b15d69e9fd19">
 is the radius of a sphere containing a charge of 1e in the Uniform Electron Gas (UEG). In the formula, <img width="214" alt="image" src="https://github.com/CodersheepY/BrianNote/assets/87512544/bb9b5212-f444-4ad1-aed0-aaf7472503ab">
. 

・Under the low-density condition, the correlation at a point in space can be written as <img width="395" alt="image" src="https://github.com/CodersheepY/BrianNote/assets/87512544/9fe646ef-4ad2-435a-9984-1c8c267dde6e">
 

Apart from UEG, other functionals that perform well include Carpley-Alder (CA), Vosko-Wilk-Nusair (VWN), Perdew-Zunger (PZ), Perdew-Wang (PW), etc."

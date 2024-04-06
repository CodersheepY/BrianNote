***一 What is Density Functional Theory***
1. **Density**: Electron Density
2. **Functional**: Energy as a function of Electron Density. By minimizing the energy functional, one can obtain the ground-state electron density and the corresponding energy of the system.
3. **Electron Density**: A function of spatial coordinates.
4. **Functional**: A function of a function.
5. **Density Functional Theory**: A method for studying the electronic structure of multi-electron systems (describing the properties and behavior of atoms and molecules) through electron density. The electronic motion equations are no longer the Schrödinger equation (whose wave function represents the solution to a many-body problem that is difficult to solve in practical situations), but are instead solved in terms of the electron density, simplifying the computational process.

***二 Hohenberg-Kohn Theorems***
1. **H-K First Theorem**: For a given electron density \( n(r) \), there exists a unique corresponding ground-state external potential \( V_{\text{ext}}(r) \) and ground-state total energy \( E \). Thus, the ground-state energy of the Schrödinger equation can be expressed in terms of the electron density.
<img width="239" alt="image" src="https://github.com/CodersheepY/BrianNote/assets/87512544/7eb08a3f-7883-4105-92f3-7576332c7070">

We can see that the degrees of freedom of the Schrödinger equation have decreased from 3N to 3, greatly reducing the complexity of computations.

Uniqueness proof:まだproof

2. **Hohenberg-Kohn Second Theorem**: The electron density that yields the lowest energy is the ground-state electron density. The second theorem can be regarded as a variational principle of the Hohenberg-Kohn theorem and can be derived through the variational principle of the Schrödinger equation.

***三  Density Functional Theory (DFT)***

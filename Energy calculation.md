## Task1: Comparing Smearing Algorithms for CO Molecule

### Objective

The goal of this task is to compare the energy changes in a CO molecule using two different smearing algorithms specified by the `ISMEAR` parameter in VASP. Specifically, we will compare the results for `ISMEAR = 0` (default) and `ISMEAR = 1` (Methfessel-Paxton smearing).

### Method

1. **Modifying INCAR File**:
   The INCAR file was modified to include the smearing algorithm by setting `ISMEAR` to `0` and `1` respectively.

   ```plaintext
   ISMEAR = 0
   ```

   ```plaintext
   ISMEAR = 1
   ```

2. **Running VASP Calculations**:
   VASP calculations were performed separately for both `ISMEAR` settings.

3. **Extracting Energy Values**:
   The energy values were extracted from the OUTCAR files using the following command:
   ```bash
   grep "energy(sigma->0)" OUTCAR
   ```

### Results

#### ISMEAR = 0 (Default)

For `ISMEAR = 0`, the energy values for different O positions are:

- **O Position 1.1 Å**:
  ```plaintext
  energy without entropy = 72.59187845  energy(sigma->0) = 72.59187845
  energy without entropy =  0.16303939  energy(sigma->0) =  0.16303939
  energy without entropy = -15.55580481  energy(sigma->0) = -15.55580481
  energy without entropy = -15.65303113  energy(sigma->0) = -15.65303113
  energy without entropy = -15.65370544  energy(sigma->0) = -15.65370544
  energy without entropy = -14.86440519  energy(sigma->0) = -14.86440519
  energy without entropy = -14.70778006  energy(sigma->0) = -14.70778006
  energy without entropy = -14.67537914  energy(sigma->0) = -14.67537914
  energy without entropy = -14.67425952  energy(sigma->0) = -14.67425952
  energy without entropy = -14.67514243  energy(sigma->0) = -14.67514243
  energy without entropy = -14.67776843  energy(sigma->0) = -14.67776843
  energy without entropy = -14.67850059  energy(sigma->0) = -14.67850059
  energy without entropy = -14.67945276  energy(sigma->0) = -14.67945276
  energy without entropy = -14.67958062  energy(sigma->0) = -14.67958062
  energy without entropy = -14.67968767  energy(sigma->0) = -14.67968767
  energy without entropy = -14.67969861  energy(sigma->0) = -14.67969861
  ```
  <img width="784" alt="Pasted Graphic 5" src="https://github.com/CodersheepY/BrianNote/assets/87512544/814af292-c587-4ed8-9d6c-af879851691f">

#### ISMEAR = 1 (Methfessel-Paxton Smearing)

For `ISMEAR = 1`, the energy values for different O positions are:

- **O Position 1.1 Å**:
  ```plaintext
  energy without entropy = -14.67970525  energy(sigma->0) = -14.67197471
  energy without entropy = -14.67970526  energy(sigma->0) = -14.67197471
  energy without entropy = -14.67970526  energy(sigma->0) = -14.67197471
  ```
  <img width="587" alt="-14 67970525" src="https://github.com/CodersheepY/BrianNote/assets/87512544/6d6b7488-0ee5-445f-a149-c0bc398c488d">

### Analysis

#### Energy Comparison

The energy values obtained from the calculations indicate that the energy without entropy and energy with sigma are slightly different for the two smearing algorithms:

- For `ISMEAR = 0`, the minimum energy is around `-14.67969861` eV.
- For `ISMEAR = 1`, the minimum energy is around `-14.67197471` eV.

#### Interpretation

- **ISMEAR = 0**:

  - This setting is generally used for molecules and small systems as it uses the Gaussian smearing method, which is suitable for such calculations.
  - The energy values are more stable and accurate for small systems like CO molecule.

- **ISMEAR = 1**:
  - This setting uses the Methfessel-Paxton scheme, which is better suited for metals and large systems with many electrons.
  - The energy values obtained with this setting show minor discrepancies and are not recommended for small systems like molecules.

### Conclusion

For the CO molecule, using the default smearing algorithm (`ISMEAR = 0`) provides more reliable and accurate energy values compared to using the Methfessel-Paxton smearing (`ISMEAR = 1`). The slight differences in energy values highlight the importance of choosing the appropriate smearing method for different types of systems. For molecular systems, the default setting (`ISMEAR = 0`) is preferred.

---

## Task2: Effect of O Position on Energy in CO Molecule

### Objective

The objective of this task is to investigate how the energy of a CO molecule changes when the position of the oxygen (O) atom is varied. Specifically, we want to verify if the energy minimum corresponds to the experimental C-O bond length of approximately 1.12 Å.

### Method

1. **Modifying the POSCAR File**:
   The POSCAR file is modified to change the position of the O atom while keeping the C atom fixed at (0.0, 0.0, 0.0). The following O positions are used: 0.9 Å, 1.0 Å, 1.1 Å, 1.2 Å, and 1.3 Å.

   ```plaintext
   CO molecule in a box
   1.0
   8.0 0.0 0.0
   0.0 8.0 0.0
   0.0 0.0 8.0
   C O
   1 1
   Cartesian
   0.0 0.0 0.0
   0.0 0.0 x
   ```

   Where `x` is varied as 0.9, 1.0, 1.1, 1.2, and 1.3 Å.

2. **Running VASP Calculations**:
   VASP calculations are performed for each O position.

3. **Extracting Energy Values**:
   The energy values are extracted from the OUTCAR files using the following command:
   ```bash
   grep "energy(sigma->0)" OUTCAR
   ```

### Results

#### Energy Values for Different O Positions

```plaintext
O Position (Å)     Energy (eV)
0.9               -12.30588360
1.0               -14.95947326
1.1               -14.67969861
1.2               -14.64170581
1.3               -13.82653313
```

#### Extracted Data

For O Position = 0.9 Å:

```plaintext
energy without entropy = -12.30588360  energy(sigma->0) = -12.30588360
energy without entropy =  -9.39615930  energy(sigma->0) =  -9.39615930
energy without entropy =  -8.61285519  energy(sigma->0) =  -8.61285519
...
energy without entropy =  -8.44945070  energy(sigma->0) =  -8.44945070
```

<img width="585" alt="-8 61285519" src="https://github.com/CodersheepY/BrianNote/assets/87512544/bbc690d2-5d17-4cb8-a402-5d9bcb057b3c">

For O Position = 1.0 Å:

```plaintext
energy without entropy = -14.95947326  energy(sigma->0) = -14.95947326
energy without entropy = -13.24585835  energy(sigma->0) = -13.24585835
energy without entropy = -13.14619441  energy(sigma->0) = -13.14619441
...
energy without entropy = -13.11980634  energy(sigma->0) = -13.11980634
```

<img width="580" alt="-13 11458679" src="https://github.com/CodersheepY/BrianNote/assets/87512544/af1facc7-6d98-4195-9585-53616f6a744d">

For O Position = 1.1 Å:

```plaintext
energy without entropy =  72.59187845  energy(sigma->0) =  72.59187845
energy without entropy =   0.16303939  energy(sigma->0) =   0.16303939
energy without entropy = -15.55580481  energy(sigma->0) = -15.55580481
...
energy without entropy = -14.67969861  energy(sigma->0) = -14.67969861
```

<img width="784" alt="Pasted Graphic 5" src="https://github.com/CodersheepY/BrianNote/assets/87512544/ade5445a-180c-45ae-8bf9-df1a9a680220">

For O Position = 1.2 Å:

```plaintext
energy without entropy = -22.78839068  energy(sigma->0) = -22.78839068
energy without entropy = -17.07771626  energy(sigma->0) = -17.07771626
energy without entropy = -17.17701353  energy(sigma->0) = -17.17701353
...
energy without entropy = -14.64170581  energy(sigma->0) = -14.64170581
```

<img width="587" alt="energy without entropy" src="https://github.com/CodersheepY/BrianNote/assets/87512544/6822fb64-aa19-41c8-a389-0d070e4b2018">

For O Position = 1.3 Å:

```plaintext
energy without entropy = -15.14833601  energy(sigma->0) = -15.14833601
energy without entropy = -14.44715307  energy(sigma->0) = -14.44715307
energy without entropy = -14.22614163  energy(sigma->0) = -14.22614163
...
energy without entropy = -13.82653313  energy(sigma->0) = -13.82653313
```

<img width="586" alt="energy without" src="https://github.com/CodersheepY/BrianNote/assets/87512544/b78b33ff-9d9e-4274-b699-c73028a48a00">

### Analysis

#### Energy Curve Plot

Using Python and Matplotlib, the energy values are plotted against the O positions to visualize the energy changes:

```python
import matplotlib.pyplot as plt

# Data
o_positions = [0.9, 1.0, 1.1, 1.2, 1.3]
energies = [-12.30588360, -14.95947326, -14.67969861, -14.64170581, -13.82653313]

# Plotting
plt.plot(o_positions, energies, marker='o', linestyle='-', color='b')
plt.axvline(x=1.12, color='r', linestyle='--', label='Experimental C-O Bond Length (1.12 Å)')
plt.title('Energy Curve for CO Molecule with Varying O Positions')
plt.xlabel('O Position (Å)')
plt.ylabel('Energy (eV)')
plt.legend()
plt.grid(True)
plt.show()
```

#### Energy Curve Visualization

<img width="854" alt="Pasted Graphic 6" src="https://github.com/CodersheepY/BrianNote/assets/87512544/a66905b5-45a3-4d01-8c0d-c2ae850b502a">

The plot shows that the energy minimum occurs at an O position of 1.1 Å with an energy value of -14.67969861 eV.

### Conclusion

From the analysis, the energy minimum at an O position of 1.1 Å closely matches the experimental C-O bond length of 1.12 Å, validating the VASP calculations and the chosen methodology. This confirms that the C-O bond length in the CO molecule is accurately predicted by the computational model.

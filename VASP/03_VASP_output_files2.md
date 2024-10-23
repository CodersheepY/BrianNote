# OUTCAR File in VASP

The **OSZICAR** file contains simple information about the system's structure optimization and the electronic structure iteration convergence. The **OUTCAR** file also contains this information, but in much greater detail. Additionally, OUTCAR includes the parameters of the calculation (some from INCAR and some default values), the K-point setup, basic information about the pseudopotentials, the VASP version, the system's geometry, electronic properties, and computation times. This makes OUTCAR a comprehensive summary of the current calculation.

VASP can compute various properties, and each specific task has corresponding output information. This leads to the complexity and variety of OUTCAR. For example, structure optimization, frequency calculation, phonon spectra, band structure, and DOS all produce specific outputs based on the needs. Although OUTCAR contains a lot of information, most of it is usually not needed.

First, the output sections of VASP are separated by long dashes (`----------------------------`). When you see these dashes, it indicates the next section of results. For example, at the start of the OUTCAR file:

```plaintext
vasp.x.x.x.xxxx (build xxxxxxx) complex

executed on xxxx date 20xx.xx.xx xx:xx:xx
running on xx total cores
distrk: each k-point on xx cores, x groups
distr: one band on NCORES_PER_BAND = x cores, xx groups
```

---

### VASP Version and Core Information

This section displays the VASP version and the number of cores being used.

---

### Next Section: System and Parameter Information

This section includes information on POTCAR, POSCAR, KPOINTS, and INCAR files.

```plaintext
-----------------------------------
| W W AA RRRRR N N II N N GGGG !!! |
| W W A A R R NN N II NN N G G !!! |
| W W A A R R N N N II N N N G !!! |
| W WW W AAAAAA RRRRR N N N II N N N G GGG ! |
| WW WW A A R R N NN II N NN G G |
| W W A A R R N N II N N GGGG !!! |
-----------------------------------

For optimal performance we recommend to set:
NCORE= 4 - approx SQRT(number of cores)
NCORE specifies how many cores store one orbital (NPAR = cpu/NCORE).
This setting can greatly improve the performance of VASP for DFT.
The default NCORE=1 might be grossly inefficient on modern multi-core architectures or massively parallel machines.
Do your own testing!
Unfortunately, you need to use the default for GW and RPA calculations.
```

---

### POTCAR Information

The following section shows the basic information about the **POTCAR** file. You can check the elements in the POTCAR by using the following commands:

```bash
grep POTCAR OUTCAR
grep TIT OUTCAR
grep ENMAX OUTCAR
grep ZVAL OUTCAR
```

For example:

```plaintext
POTCAR: PAW_PBE O xxxxxx
VRHFIN = O: s2p4
LEXCH = PE
EATOM = 432.3788 eV, 31.7789 Ry
TITEL = PAW_PBE O xxxxxx
```

- **ZVAL** shows the number of valence electrons for the element, in this case, 6 for the oxygen atom. You can compare this to the content in the POTCAR file.

---

### POSCAR Information

Next, you will find the **POSCAR** information, which includes the format of the coordinates, atom positions, and the cell dimensions.

For example:

```plaintext
ion position nearest neighbor table
1 0.000 0.000 0.000

LATTYP: Found a simple cubic cell.
ALAT = 8.0000000000
```

Lattice vectors:

```plaintext
A1 = ( 8.0000000000, 0.0000000000, 0.0000000000)
A2 = ( 0.0000000000, 8.0000000000, 0.0000000000)
A3 = ( 0.0000000000, 0.0000000000, 8.0000000000)
```

---

### Symmetry Analysis

VASP analyzes the symmetry of the system and provides information about the space group and point operations. For example, in this system, the point symmetry is **O_h**, and it has **48 symmetry operations**.

```plaintext
Routine SETGRP: Setting up the symmetry group for a simple cubic supercell.
Subroutine GETGRP returns: Found 48 space group operations (whereof 48 operations were pure point group operations).
```

---

### K-point Information

To check the number of K-points, use:

```bash
grep irreducible OUTCAR
grep irre OUTCAR
```

Example output:

```plaintext
Subroutine IBZKPT returns following result:
Found 1 irreducible k-points:
Reciprocal coordinates: 0.000000 0.000000 0.000000 1.000000
```

---

### INCAR Information

This section shows the parameters from the **INCAR** file:

```plaintext
SYSTEM = O atom
POSCAR = O atom in a box

Startparameter for this run:
NWRITE = 2 write-flag & timer
PREC = normal normal or accurate (medium, high, low for compatibility)
ENCUT = 400.0 eV
```

VASP summarizes the main inputs and parameters for the calculation, including the system size, K-points, and memory requirements, before starting the actual computation.

---

### Iteration Information

During the calculation, each iteration's output is displayed. For example:

```plaintext
--------------------------------------- Iteration 1( 1) ---------------------------------------

    POTLOK:  cpu time    0.0240: real time    0.0324
    SETDIJ:  cpu time    0.0030: real time    0.0040
    EDDAV:  cpu time    0.0320: real time    0.0397
    --------------------------------------------
      LOOP:  cpu time    0.0610: real time    0.0781
```

Here, you can track the changes in total energy, the number of electrons, and the magnetization. The **free energy of the ion-electron system** is also calculated:

```plaintext
Free energy TOTEN = 32.49699652 eV
```

---

### Extracting Energy Information

To extract the energy value, use:

```bash
grep without OUTCAR | tail -n 1
grep ' without' OUTCAR | tail -n 1
grep sigma OUTCAR | tail -n 1
```

These commands will help you retrieve the final energy value from the last step of the calculation.

---

### Memory Usage and Timing

At the end of the OUTCAR file, VASP outputs the memory usage and computation times, indicating the job has finished successfully:

```plaintext
writing wavefunctions
LOOP+: cpu time 2.0037: real time 2.1098
# total amount of memory used by VASP MPI-rank0 35846. kBytes
# General timing and accounting information for this job:
Total CPU time used (sec): 5.492
Elapsed time (sec): 5.931
```

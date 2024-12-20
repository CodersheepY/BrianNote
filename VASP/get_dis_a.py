#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
from ase.io import read, write
import os, sys

'''
Important Note: Atoms are indexed from 0.
'''
atom_A, atom_B = [int(i) for i in sys.argv[1:3]]
model = read('POSCAR')
positions = model.get_positions()
dis_AB_1 = model.get_distance(atom_A, atom_B, mic=True)
dis_AB_2 = model.get_distance(atom_A, atom_B, mic=False)

print(round(dis_AB_1, 4), '\t',round(dis_AB_2,4))

## Not Used 
#coord_A = positions[atom_A]
#coord_B = positions[atom_B]
#dis_AB = np.linalg.norm(coord_A-coord_B)
#print(round(dis_AB,4))  
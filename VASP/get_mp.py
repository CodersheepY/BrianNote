from pymatgen.ext.matproj import MPRester
from pymatgen.io.cif import CifWriter

mpr = MPRester("kzum4sPsW7GCRwtOqgDIr3zhYrfpaguK")

data_structure = mpr.get_structure_by_material_id("mp-1009018")

w = CifWriter(data_structure)
w.write_file("Cu.cif")

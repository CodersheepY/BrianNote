from mp_api.client import MPRester
from pymatgen.analysis.phase_diagram import PhaseDiagram, PDPlotter
from pymatgen.entries.mixing_scheme import MaterialsProjectDFTMixingScheme

with MPRester("kzum4sPsW7GCRwtOqgDIr3zhYrfpaguK") as mpr:

    # Obtain GGA, GGA+U, and r2SCAN ComputedStructureEntry objects
    entries = mpr.get_entries_in_chemsys(elements=["Li", "Fe", "O"], 
                                         additional_criteria={"thermo_types": ["GGA_GGA+U", "R2SCAN"]}) 
    
    # Apply corrections locally with the mixing scheme
    scheme = MaterialsProjectDFTMixingScheme()
    corrected_entries = scheme.process_entries(entries)
    
    # Construct phase diagram
    pd = PhaseDiagram(corrected_entries)
    
    # Plot phase diagram
    PDPlotter(pd).get_plot()

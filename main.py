from Bio import PDB
from Bio.PDB import Superimposer
import sys
import numpy as np

arg1 = sys.argv[1] #reference file path
arg2 = sys.argv[2] #modele

parser = PDB.PDBParser(QUIET=True)

reference = parser.get_structure('reference', arg1)
structures = parser.get_structure("structures", arg2)

models = []
for model in structures:
    models.append(model)

rmsd_values = []
superimposer = Superimposer()

for model in models:
    ref_atoms = []
    sub_atoms = []

    for ref_atom in reference.get_atoms():
        for sub_atom in model.get_atoms():
            if ref_atom.get_parent().get_id()[1] == sub_atom.get_parent().get_id()[1] and ref_atom.get_id() == sub_atom.get_id():
                ref_atoms.append(ref_atom)
                sub_atoms.append(sub_atom)

    superimposer.set_atoms(ref_atoms, sub_atoms)
    superimposer.apply(model.get_atoms())

    sub_coord = np.array([atom.coord for atom in sub_atoms])
    ref_coord = np.array([atom.coord for atom in ref_atoms])
    dif = sub_coord - ref_coord
    sum_of = np.sum(dif ** 2)
    rmsd = np.sqrt(sum_of / len(ref_coord))
    rmsd_values.append(rmsd)


    print(rmsd)

    ref_atoms.clear()
    sub_atoms.clear()

'''best = np.argmin(rmsd_values)
worst = np.argmax(rmsd_values)
print(best + 1)
print(worst + 1)'''







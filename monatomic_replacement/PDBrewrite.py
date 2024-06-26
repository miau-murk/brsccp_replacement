from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.Residue import Residue
from Bio.PDB.Chain import Chain
from Bio.PDB.Model import Model
from Bio.PDB.Structure import Structure

# Define a function to read a PDB file and return a structure object
def read_pdb(file_path):
    parser = PDBParser()
    structure = parser.get_structure("structure", file_path)
    return structure

# Define a function to write a structure object into a PDB file
def write_pdb(structure, output_file_path):
    io = PDBIO()
    io.set_structure(structure)
    io.save(output_file_path)

# Define a function to get information about atom connectivity
def get_atom_connectivity(structure):
    connectivity = []
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in sorted(residue.get_unpacked_list(), key=lambda x: x.get_serial_number()):
                    for neighbor in atom.get_neighbors():
                        connectivity.append((atom.get_serial_number(), neighbor.get_serial_number()))
    return connectivity

# For an example manipulation, let's define a function to renumber the atoms
def renumber_atoms(structure):
    atom_number = 1
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in sorted(residue.get_unpacked_list(), key=lambda x: x.get_serial_number()):
                    atom.serial_number = atom_number
                    atom_number += 1

# Example usage to read, manipulate, and write a PDB file
# pdb_file_path = "3KQS.pdb"
# output_file_path = "_3kqs_.pdb"

# Reading the PDB file
# structure = read_pdb(pdb_file_path)

# for i in structure[0]["A"][33]:
#     print(i.serial_number)

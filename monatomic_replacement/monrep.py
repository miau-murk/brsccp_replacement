from rdkit import Chem
import numpy as np
from .pdb_to_mol import pdb_to_mol
from .morgan_types import *
from .smarts import find_CXYZ
import yaml
import sys

def read_input_file(file_name): # function for reading yaml file
    with open(file_name, "r", encoding="utf-8") as fff:
        try:
            return yaml.safe_load(fff)
        except yaml.YAMLError as exc:
            log.error("[ERROR] No such file or directory:", file_name)
            sys.exit() 

def find_index_by_atom(pdb_lines, atom_name):
    """ This function finds atom index by atom name in pdb ligand file """
    for line in pdb_lines:
        if line.startswith("HETATM") and line[12:16].strip() == atom_name:
            return int(line[6:11].strip())
    return None

# this function is no longer used
def find_index_by_atom1(pdb_file, atom_name):
    """ Finds the atom name in the PDB file corresponding to the provided atom index """
    atom_name = None
    with open(pdb_file, 'r') as file:
        for line in file:
            if line.startswith("HETATM") or line.startswith("ATOM"):
                current_index = int(line[6:11].strip())
                if current_index == atom_name:
                    atom_name = line[12:16].strip()
                    return atom_name

def find_atom_by_index(pdb_lines, atom_index):
    """ This function finds atom name by atom index in pdb ligand file """
    for line in pdb_lines:
        if line.startswith("HETATM") and int(line[8:11].strip()) == atom_index:
            return line[12:16].strip()
    return None

def YAMLrewrite(original_file, new_file_name, new_r_atom_name):

    # Read the original YAML file
    with open(original_file, 'r') as file:
        data = yaml.safe_load(file)

    # Modify the 'r_atom_name' with the new value
    data['r_atom_name'] = new_r_atom_name

    # Write the modified YAML data back to a new file
    with open(new_file_name, 'w') as file:
        yaml.safe_dump(data, file)

def replacement_yaml(original_file):

    # read yaml file
    input_file = read_input_file(original_file)
    input_pdb = input_file["pdb_file"]
    element_before = input_file["r_from_element"]
    element_after = input_file["r_to_element"]
    repl_residue = input_file["r_resid"]
    ligand_name = input_file["r_resname"]
    atom_name = input_file["r_atom_name"]
    chain = input_file["r_chain"]

    # read pdb and create mol file
    pdb_to_mol(input_pdb, ligand_name, chain)


    # find index of replaced atom
    with open("ligand.pdb", 'r') as file:
        pdb_lines = file.readlines()

    # Extract the index of the atom named 'HP21'
    index = find_index_by_atom(pdb_lines, atom_name)

    # read mol file find topology identical atoms
    mol = Chem.MolFromMolFile("ligand.mol", removeHs=False) # create mol object

    strN = find_CXYZ(mol, element_before, element_after)

    Mol1 = MorganTypes(mol, removeH=True) # create MolTypes object without Hs
    Mol1.calc_morgan_weight() # calculate Morgan's weights without Hs
    Mol2 = MorganTypes(mol,removeH=False, StereoNone_cont=strN) # create MolTypes object with Hs

    # inherit weights
    for atom in Mol1.mgc.keys(): 
        Mol2.mgc[atom] += Mol1.mgc[atom]

    # assign the Hs atom the weights of the atoms with which they are connected
    for atom in Mol2.mgc.keys():
        if Mol2.xyz[atom][1] == "H":
            neighbor = next(iter(Mol2.graph[atom])) # find index of atom the C wich is connected with the atom H
            Mol2.mgc[atom] += Mol2.mgc[neighbor] # assign the Hs atom the weight of neighbor


    Mol2.calculation_induced_stereochemistry() # calculate induced streochemistry

    # assign the induced streochemistry in Morgan's weights
    for i in Mol2.mgc.keys():
        Mol2.mgc[i] += Mol2.ist[i]

    # find replaced atom
    equiv_atom_indeces = [] # lisst that contains equivalent atoms
    for cort in Mol2.find_topology_indentical_atoms():
        if index in cort:
            ind1, ind2 = cort[0], cort[1]
            if ind1 != index:
                equiv_atom_indeces.append(ind1)
            else:
                equiv_atom_indeces.append(ind2)

    # find name of replaced atoms in pdb ligand files
    equiv_atoms_names = []
    for ind in equiv_atom_indeces:
        equiv_atoms_names.append(find_atom_by_index(pdb_lines, ind))

    # create mew YAML files
    for q, repl in enumerate(equiv_atoms_names):
        new_file = original_file.split('.yaml')[0] + "_" + str(q+1) + '.yaml'
        YAMLrewrite(original_file, new_file, repl)

##

def find_identical_atoms_v1(original_file): 

    # read yaml file
    input_file = read_input_file(original_file)
    input_pdb = input_file["pdb_file"]
    element_before = input_file["r_from_element"]
    element_after = input_file["r_to_element"]
    repl_residue = input_file["r_resid"]
    ligand_name = input_file["r_resname"]
    atom_name = input_file["r_atom_name"]
    chain = input_file["r_chain"]

    # read pdb and create mol file
    pdb_to_mol(input_pdb, ligand_name, chain)


    # find index of replaced atom
    with open("ligand.pdb", 'r') as file:
        pdb_lines = file.readlines()

    # Extract the index of the atom named 'HP21'
    index = find_index_by_atom(pdb_lines, atom_name)

    # read mol file find topology identical atoms
    mol = Chem.MolFromMolFile("ligand.mol", removeHs=False) # create mol object

    strN = find_CXYZ(mol, element_before, element_after)

    Mol1 = MorganTypes(mol, removeH=True) # create MolTypes object without Hs
    Mol1.calc_morgan_weight() # calculate Morgan's weights without Hs
    Mol2 = MorganTypes(mol,removeH=False, StereoNone_cont=strN) # create MolTypes object with Hs

    
    # inherit weights
    for atom in Mol1.mgc.keys(): 
        Mol2.mgc[atom] += Mol1.mgc[atom]

    # assign the Hs atom the weights of the atoms with which they are connected
    for atom in Mol2.mgc.keys():
        if Mol2.xyz[atom][1] == "H":
            neighbor = next(iter(Mol2.graph[atom])) # find index of atom the C wich is connected with the atom H
            Mol2.mgc[atom] += Mol2.mgc[neighbor] # assign the Hs atom the weight of neighbor


    Mol2.calculation_induced_stereochemistry() # calculate induced streochemistry

    # assign the induced streochemistry in Morgan's weights
    for i in Mol2.mgc.keys():
        Mol2.mgc[i] += Mol2.ist[i]

    # find replaced atom
    equiv_atom_indeces = [] # lisst that contains equivalent atoms
    for cort in Mol2.find_topology_indentical_atoms():
        if index in cort:
            ind1, ind2 = cort[0], cort[1]
            if ind1 != index:
                equiv_atom_indeces.append(ind1)
            else:
                equiv_atom_indeces.append(ind2)

    # find name of replaced atoms in pdb ligand files
    equiv_atoms_names = []
    for ind in equiv_atom_indeces:
        equiv_atoms_names.append(find_atom_by_index(pdb_lines, ind))
    
    return equiv_atoms_names

def find_identical_atoms_v2(original_file):
    pass

def find_all_identical_atoms(input_pdb_file):
    pass





# read input YAML file
#original_file = sys.argv[1]

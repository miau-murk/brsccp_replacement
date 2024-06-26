from rdkit import Chem
from rdkit.Chem import AllChem

# Function to find R-CXYZ groups in a molecule
def find_R_CXYZ_groups(mol, elements):
    indices = []  # List to store indices of C in R-CXYZ

    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C' and atom.GetDegree() == 4:  # Carbon must be connected to 4 other atoms
            neighbors = atom.GetNeighbors()
            neighbors_elements = []
            for i in neighbors:
                if i.GetSymbol() != "C":
                    neighbors_elements.append(i.GetSymbol())
            if len(neighbors) == 4 and set(elements) == set(neighbors_elements):  # Ensure no hydrogens are directly bonded
                for i in neighbors:
                    if i.GetSymbol() in elements:
                        indices.append(i.GetIdx()+1)

    return indices

# Main function to process the mol file and find R-CXYZ groups
def find_CXYZ(mol, element_before, element_after):
    if not mol:
        print("Error in loading the molecule from the file, check your PDB file")
        return
    
    elements_cont = [element_before, element_before, element_after]

    # Find R-CXYZ
    r_cxyz_indices = find_R_CXYZ_groups(mol, elements_cont)

    return r_cxyz_indices

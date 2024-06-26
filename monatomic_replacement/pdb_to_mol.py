from Bio.PDB import PDBParser, PDBIO, Select
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdmolfiles
import re


# this two functions reorder atoms in mol file, it's neccesery for correct calculation morgan types
def reorder_hydrogens(mol):
    """Reorders hydrogen atoms to follow heavy atoms in the molecule."""
    # Find heavy atoms and hydrogen atoms indices
    heavy_atom_indices = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1]
    hydrogen_atom_indices = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 1]
    
    # Get the ordering with heavy atoms first followed by hydrogens
    new_ordering = heavy_atom_indices + hydrogen_atom_indices
    return Chem.rdmolops.RenumberAtoms(mol, new_ordering)

def read_and_reorder_mol_file(input_file, output_file):
    """Reads a .mol file, reorders hydrogen atoms, and writes the result to a new .mol file."""
    # Read the input .mol file
    mol = Chem.MolFromMolFile(input_file, removeHs=False)
    if mol is None:
        raise ValueError('Could not read the molecule from the file.')
    
    # Reorder hydrogen atoms to be last
    mol = reorder_hydrogens(mol)
    
    # Write the output .mol file
    Chem.MolToMolFile(mol, output_file)\
    
# Function to determine if an atom is a hydrogen atom based on its name
def is_hydrogen(atom_name):
    # Hydrogens in PDB files start with 'H' or digit for hydrogen naming convention (e.g., "1HD1").
    return atom_name.strip()[0] == 'H'

# Function to parse and process PDB lines.
def process_pdb_lines(pdb_lines):
    heavy_atoms = []
    hydrogens = []

    # Process all HETATM lines
    for line in pdb_lines:
        # Check if the line starts with HETATM
        if line.startswith("HETATM"):
            atom_name = line[12:16]  # The atom name is columns 13-16 in the PDB file format
            if is_hydrogen(atom_name):
                hydrogens.append(line)
            else:
                heavy_atoms.append(line)
    
    # Combine heavy atoms and hydrogens
    ordered_lines = heavy_atoms + hydrogens
    return ordered_lines

# Function to reindex the atoms and write a new PDB file
def write_reordered_pdb(pdb_lines, ordered_lines, output_filename):
    with open(output_filename, 'w') as output_file:
        atom_index = 1  # Initialize atom index

        # Go through each line and modify atom index
        for line in ordered_lines:
            # Here, we use string formatting to construct the new line with proper spacing
            new_line = "{record_type:6}{atom_index:5d} {atom_name:<4}{alt_loc:1}{residue_name:3} {chain_id:1}{residue_number:>4}{insertion_code:1}   {x:>8.3f}{y:>8.3f}{z:>8.3f}{occupancy:>6.2f}{temp_factor:>6.2f}          {element:>2}{charge:2}\n"
            new_line = new_line.format(
                record_type=line[0:6],
                atom_index=atom_index,
                atom_name=line[12:16],
                alt_loc=line[16],
                residue_name=line[17:20],
                chain_id=line[21],
                residue_number=int(line[22:26]),
                insertion_code=line[26],
                x=float(line[30:38]),
                y=float(line[38:46]),
                z=float(line[46:54]),
                occupancy=float(line[54:60]),
                temp_factor=float(line[60:66]),
                element=line[76:78].strip(),
                charge=line[78:80]
            )
            output_file.write(new_line)
            atom_index += 1

        # Write the "TER   " and "END" lines if they were in the original PDB
        if any(line.startswith("TER") for line in pdb_lines):
            output_file.write("TER\n")
        if any(line.startswith("END") for line in pdb_lines):
            output_file.write("END\n")
    


# Define a class to select the ligand
class LigandSelect(Select):
    def __init__(self, ligand_resname):
        self.ligand_resname = ligand_resname

    def accept_residue(self, residue):
        return residue.get_resname() == self.ligand_resname
    
def pdb_to_mol(pdb_name, ligand_name, chain):


    # Load the PDB file
    pdb_filename = pdb_name # replace with your PDB file
    ligand_resname = ligand_name    # replace with the residue name of your ligand
    pdb_parser = PDBParser()
    structure = pdb_parser.get_structure('structure', pdb_filename)

    # Extract the ligand
    ligand_selector = LigandSelect(ligand_resname)
    io = PDBIO()
    io.set_structure(structure[0][chain])
    io.save('ligand.pdb', ligand_selector)

    # Convert the PDB to an RDKit molecule
    ligand_pdb_text = open('ligand.pdb').read()
    rdkit_mol = Chem.MolFromPDBBlock(ligand_pdb_text, sanitize=False, removeHs=False, flavor=0)

    # Generate the ligand MOL file, preserving the connections
    rdkit_mol_with_bonds = Chem.rdmolfiles.MolToMolBlock(rdkit_mol)
    with open('ligand.mol', 'w') as ligand_mol_file:
        ligand_mol_file.write(rdkit_mol_with_bonds)
    
    # You can call the function with the input and output file names
    input_file = 'ligand.mol'  # Replace this with the path to the input .mol file
    output_file = 'ligand.mol'  # Replace this with the path to the output .mol file
    read_and_reorder_mol_file(input_file, output_file)

    with open('ligand.pdb', 'r') as original_pdb_file:
        pdb_lines = original_pdb_file.readlines()

    ordered_lines = process_pdb_lines(pdb_lines)
    write_reordered_pdb(pdb_lines, ordered_lines, 'ligand.pdb')
    #print("Ligand successfully extracted and written to ligand.mol")
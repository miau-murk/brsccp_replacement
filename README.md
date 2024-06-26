# brsccp_replacement
## A module for searching for equivalent positions in a molecule

The module is designed to find topologically equivalent atoms in a molecule. 
We will call topologically equivalent atoms those, when replaced by an arbitrary atom X, identical molecules are obtained.

Make sure that ```rdkit```, ```Bio```, ```periodictable```, ```yaml``` modules are installed in your environment

```conda create -c conda-forge -n my-rdkit-env rdkit```
```conda activate my-rdkit-env```
```conda install periodictable, yaml, Bio```


To work with the code, you just need to put the ```monatomic_replacement``` folder in your main directory.

# Basic functions:

1. To find equivalent atoms in a ligand using in the BRSCCP program. Make sure that the pdb file of the complex and the yaml input file are located in your main folder.

```import monatomic_replacement```

```atoms_list = monatomic_replacement.find_identical_atoms("input.yaml")```

2. You can import the Morgan Types class, which makes it easy to find equivalent positions in a molecule, compare mol and smiles files for identity. Here are the main methods:

- Import modules:

```from monatomic_replacement import morgan_types```
```from rdkit import Chem```

```mol = Chem.MolFromMolFile("ligand.mol") # create RDKit mol object```

```mnt = morgan_types.MorganTypes(mol) # create MorganTypes objesct```

- If you need to compare the typological identity of atoms in a molecule, use the method:

```mnt.calc_morgan_weight()```

```>> {1: 1549527, 2: 2073202, 3: 1283756, 4: 954312, 5: 954312, 6: 1283756, 7: 2073202, 8: 1549527, 9: 1595262, 10: 675260}```

The numbers are the numbers of atoms in the mol file.

- To find all topologically identical atoms in a molecule, use the function:

```mnt.find_topology_indentical_atoms()```

- Mol files may differ in the order and coordinates of atoms, therefore, to determine the identity of molecules from mol files, Morgan types can be compared. To determine the Morgan type of molecule, use the function:

```mnt.calc_morgan_type()```

3. Work with pdb and mol files. If you want to get the pdb and mol files of the ligand contained in the pdb file of the complex, use the following function:

```pdb_to_mol("input.pdb", ligand_name, chain)```

The pdb and mol files of the ligand will be created in the folder where the pdb file of the complex is located.

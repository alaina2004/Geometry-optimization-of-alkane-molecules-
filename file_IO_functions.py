import numpy as np

# Define a function that open a mol.2 file and reasall the information for x molecule
def read_mol2(filepath):
# udemy course python --> open files in general 
    with open(filepath, 'r') as file:
        first_line = file.readline() # to read first line python is 0 
        first_line = first_line.split() #splits the line into a list 
        nat = int(first_line[0]) #total number of atoms 
        nbonds = int(first_line[1]) #total number of bonds 
        nC = int(first_line[2])  #number of carbon atoms 
        nCC = int(first_line[3])  #number of C-C bond 

        atomic_coordinates = np.zeros((nat, 3)) # each atom will have their x y z coordinate
        bonds = np.zeros((nbonds, 5)) # each bond will be defined by 5 parameters: atom_i, atom_j, bond-type, bond-length and bond energy 
        atoms = [] # for the array of atomic symbols

        for i in range(nat):
            coordinate_line = file.readline()
            coordinate_line = coordinate_line.split()
            atom = coordinate_line[3] #because in the third column the symbol
            atoms.append(atom) # append = to add 

            x = float(coordinate_line[0])
            y = float(coordinate_line[1])
            z = float(coordinate_line[2])

            atomic_coordinates[i] = [x, y, z] 

       
        for i in range(nbonds):
            bond_line = file.readline()
            bond_line = bond_line.split()

            atom_x = int(bond_line[0])   #bond_line reading the whole text line guidance
            atom_y = int(bond_line[1])

            atom_x_type = atoms[atom_x - 1]
            atom_y_type = atoms[atom_y - 1]
            bond_type = atom_x_type + atom_y_type # creates either a "CC" or "CH" string
            bond_type = 0 if bond_type == "CC" else 1 # updates the bond type to be 0 (C-C) or 1 (C-H) for simplicity

            bonds[i, :3] = [atom_x, atom_y, bond_type] #3 because the three row and 2 other fill later (total5)

    return nat, nbonds, nC, nCC, atoms, atomic_coordinates, bonds
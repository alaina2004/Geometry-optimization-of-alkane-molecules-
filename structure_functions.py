import numpy as np

from auxiliary import *

# to obtain r_b equation 21 manual
def calculate_bond_lenghts(bonds, atomic_coordinates): 
    #takes the bonds arrays and modify it fills the empty column from file_IO_functions.py. It does not return only modifies input. 
    for bond in bonds:
        
        atom_a = int(bond[0])   
        atom_b = int(bond[1])    

        a_coords = atomic_coordinates[atom_a - 1]  #to translate to python language not mine for me first floor is python's cero floor (basement)
        b_coords = atomic_coordinates[atom_b - 1]

        bond_length = np.linalg.norm(a_coords - b_coords) #np.linalg.norm is to get the distance for vectors (norm)

        bond[3] = bond_length   #ahora para poner en la cuarta columna del array que en python es 3 

def calculate_bond_stretching_energies(bonds):
#fills the fith column (bond energies)
    for bond in bonds:
        bond_type = int(bond[2])   # C-C or C-H different force constant and parametros
        k_b = bond_force_constants.get(bond_type) #auxiliar
        r_0 = equilibrium_lenghts.get(bond_type)
        r_b = bond[3]

        energy = k_b * (r_b - r_0)**2 
        bond[4] = energy 

def find_angles(bonds, atoms, nC):
    #angles are going to be found here and store in the new array "angles" // to find which atoms make angles with who not the value degree
    angles = np.zeros((nC*6, 6)) # 6 angles for each C atom, 6 numbers defining an angle

    # easiest way for me to make all the possible commbinations of atoms that make an angle 
    bond_dict = {}
    for bond in bonds:
        atom_a, atom_b = int(bond[0]), int(bond[1])
        if atom_a not in bond_dict:
            bond_dict[atom_a] = []
        if atom_b not in bond_dict:
            bond_dict[atom_b] = []
        bond_dict[atom_a].append(atom_b)
        bond_dict[atom_b].append(atom_a)
    # filters all the hydrogens and their single bonding partner    
    bond_dict = {key: value for key, value in bond_dict.items() if len(value) == 4 }
        
    k = 0  
    for carbon in bond_dict:
        bonded_atoms = bond_dict.get(carbon)
        for i in range(4): # in alkanes, C must have exactly 4 bonds
            atom_i = bonded_atoms[i]
            for j in range(i+1, 4): 
                atom_j = bonded_atoms[j]

                angle_type = 0 if atoms[atom_i - 1] == 'H' or atoms[atom_j - 1] == 'H' else 1
                angles[k, :4] = [atom_i, carbon, atom_j, angle_type] # gets 
                k += 1
    
    return angles, bond_dict

def calculate_angles(angles, atomic_coordinates):
    # angle X-C-Y // here yes, the degrees are found with this function // theta a is found here 
    for angle in angles:
        X = atomic_coordinates[int(angle[0]) - 1]  #everything stored in  arrays// array same nature int, float or string // atom are integers
        C = atomic_coordinates[int(angle[1]) - 1]
        Y = atomic_coordinates[int(angle[2]) - 1]
        #to this point calculate theta for eq. 21. Theta cero y constante are given
        CX = X - C
        CY = Y - C  #obtaining vectors (tijeras, C central atom)

        CX_magnitude = np.linalg.norm(CX) #formula
        CY_magnitude = np.linalg.norm(CY)
        CX_dot_CY = np.dot(CX, CY)

        cos_theta = CX_dot_CY / (CX_magnitude * CY_magnitude)  # All of this are formulas to obtain 
        theta = np.arccos(cos_theta)
        theta = np.degrees(theta)

        angle[4] = theta # it is going to be store in the array created before in def find_angles on the fourth python column 

def calculate_angle_bending_energies(angles):

        for angle in angles:
            angle_type = int(angle[3]) #based on which angle depending on the atoms connected 
            k_a = angle_force_constants.get(angle_type) # given 
            theta_a = angle[4] #just calculated def calculate_angles
            
            energy = k_a * (theta_a - theta_0)**2   # Manual Eq.21 
            angle[5] = energy                       #where to store 
            


def find_dihedrals(bonds, bonds_dict, nCC):
    # to create array for dihedrals 
    dihedrals = np.zeros((nCC*9,6)) # (9 dihedrals for each C-C bond) x (4 atoms defining the dihedral angle, dihedral angle phi, torsion energy E)
    
    # here, a dihedral will be noted as X-Y-Z-Q, where Y and Z are C atoms of the C-C bond around which the dihedral is formed
    # here all the combination between four atoms that for a dihedral are found
    k = 0  # k where the first is located
    for i in range(nCC):
        CC_bond = bonds[i] 
        Y, Z = int(CC_bond[0]), int(CC_bond[1]) # C1 and C2  , integer transform and found 

        X_list = [X for X in bonds_dict[Y] if X != Z] # gets all the atoms connected to C-C (example H in ethane), except the other C in that bond, already in the z
        Q_list = [Q for Q in bonds_dict[Z] if Q != Y] # same as above, resulting in list of three atoms (4th atom would be the C-C partner)

        for X in X_list:
            for Q in Q_list:
                dihedrals[k, :4] = [X, Y, Z, Q] # adds the k-th dihedral to dihedrals array
                k += 1 
    
    return dihedrals

def calculate_dihedral_angles(dihedrals, atomic_coordinates):

    # X-Y-Z-Q notation for dihedrals will be used; Y-Z is the C-C bond

    for dihedral in dihedrals:

        # gets the coordinates of X, Y, Z, Q
        X = atomic_coordinates[int(dihedral[0]) - 1]
        Y = atomic_coordinates[int(dihedral[1]) - 1]
        Z = atomic_coordinates[int(dihedral[2]) - 1]
        Q = atomic_coordinates[int(dihedral[3]) - 1]

        # determine the vectors pointing from X to Y, from Y to Z, from Z to Q
        XY = Y - X 
        YZ = Z - Y
        ZQ = Q - Z 

        # calculate the normal vectors
        t = np.cross(XY, YZ) # vector t perpendicular to the X-Y-Z plane // eq.26 - eq.28 manual 
        u = np.cross(YZ, ZQ) # vector u perpendicular to the Y-Z-Q plane  
        v = np.cross(t, u) # vector v perpendicular to both t and u

        #dot products and magnitudes
        t_dot_u = np.dot(t, u)
        YZ_dot_v = np.dot(YZ, v)

        mag_t = np.linalg.norm(t)
        mag_u = np.linalg.norm(u)
        mag_YZ = np.linalg.norm(YZ)

        # calculating the dihedral angle phi
        cos_phi = t_dot_u / (mag_t * mag_u)
        sin_phi = YZ_dot_v / (mag_YZ * mag_t * mag_u)
        phi = np.arctan2(sin_phi, cos_phi)
        phi = np.degrees(phi)

        # add the phi angle to the dihedrals array
        dihedral[4] = phi

def calculate_torsion_energies(dihedrals):
    
    for dihedral in dihedrals:

        phi = dihedral[4]

        energy = A_phi * (1 + np.cos(np.radians(n_phi * phi)))

        dihedral[5] = energy

def find_vdW_pairs(nat, bonds, angles, atoms):
   
    vdW = [] # vdW Van der Waals forces 

    # following nested loop for i,j finds all unique pairs of atoms -----> nested term for a for loop inside another for loop 
    for i in range(1, nat + 1): # range 2 then python 0, 1 not 2 
        for j in range(i + 1, nat + 1):  # i +1 for Van der Wals 

            # assume the pair is not participating in any bond or angle together, i.e. it is a vdW pair
            is_vdW_pair = True

            # if i and j are found in a bond, is_vdW_pair is set to False and break ends the loop as there is no need to look further
            for bond in bonds:
                if i in bond[:2] and j in bond[:2]:   # :2 to check column one and two of array bonds // atom_a and atom_b 
                    is_vdW_pair = False
                    break

            # only looks if i and j are in angle if they were not found to form a bond
            if is_vdW_pair:
                for angle in angles:
                    if i in angle[:3] and j in angle[:3]:  # :3 
                        is_vdW_pair = False
                        break
            
            # if given pair is a vdW pair, determines whether it is a CH, CC or HH pair and adds the pair and type to the vdW array
            if is_vdW_pair:
                
                # determines which pair type it is
                pair_types = {"CC" : 0, "CH" : 1, "HC" : 1, "HH" : 2} 
                i_is, j_is = atoms[i - 1], atoms[j - 1] # list of symbols to find which atom is C or H 
                pair = i_is + j_is # to add them together, strings
                pair_type = pair_types.get(pair) # to classify them depending on the pair C-H (1), C-C (0) and H-H(2)
                
                vdW.append([i, j, pair_type, 0]) # 0 is just a placeholder value for vdW energy // array 4 columns // row don't know how many until there all found 
    
    vdW = np.array(vdW, dtype=float) # energy going to be a float 
    return vdW

def calculate_vdW_energies(vdW_pairs, atomic_coordinates):

    for pair in vdW_pairs:

        pair_type = int(pair[2])

        i_coords = atomic_coordinates[int(pair[0]) - 1] # every time i use, most change to integer because before it is converted into float because of the energy 
        j_coords = atomic_coordinates[int(pair[1]) - 1]
        r_ij = np.linalg.norm(i_coords - j_coords)  # formula of numpy

        A = vdW_A.get(pair_type)
        B = vdW_B.get(pair_type)

        energy = A/r_ij**12 - B/r_ij**6
        pair[3] = energy
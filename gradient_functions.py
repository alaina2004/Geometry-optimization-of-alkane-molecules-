import numpy as np

from auxiliary import *

def calculate_bond_stretching_gradients(nat, bonds, atoms, atomic_coordinates):
# derivatives this involver eq.29 and 35 
# dV/dx = dV/dr * dr/dx
    g_bonds = np.zeros((nat, 3))  # g is for gradients, to create array, the columns are x , y and z  for each vector 

    for bond in bonds: # 

        bond_type = int(bond[2]) # just getting the needed variables 
        k_b = bond_force_constants.get(bond_type) # ""
        r_0 = equilibrium_lenghts.get(bond_type)  # ""
        r_b = bond[3] # ""
        dVdr = 2 * k_b * (r_b - r_0) # derivative of potential in respect bond length  so only one derivative eq. 35 // just one of the terms 
                                    
        atom_a = int(bond[0]) # calling the atoms from array // two atoms because of bond 
        atom_b = int(bond[1])
        a_coords = atomic_coordinates[atom_a - 1] 
        b_coords = atomic_coordinates[atom_b - 1]

        drda = (a_coords - b_coords) / r_b # Eq.29 here // also the second term multiplying in eq.35 
        drdb = -1 * drda # the gradient for a bond a and b is the same but different direction 
 

        g_bonds[atom_a - 1] += dVdr * drda # This is the multiplication of both terms "red" and "blue"
        g_bonds[atom_b - 1] += dVdr * drdb # dr/da an dr/db form a bond 

    return g_bonds

def calculate_angle_bending_gradients(nat, nbonds, angles, atoms, atomic_coordinates):
# g bend(r) eq.37 
    g_angles = np.zeros((nat, 3)) # array same logic as i did before 

    for angle in angles:

        angle_type = int(angle[3])

        k_a = angle_force_constants.get(angle_type) # given 
        theta_a = angle[4] # already obtained 
        
        # the gradient of energy in respect to angle theta (kcal mol-1 deg-1)
        dVdtheta = 2 * k_a * (theta_a - theta_0) 
        # the gradient of energy in respect to cartesian coordinates that will be
        # calculated below is in rad-1; a conversion factor is added to f below 
        dVdtheta = dVdtheta * (180 / np.pi)


        atom_X = int(angle[0])
        atom_C = int(angle[1])
        atom_Y = int(angle[2])
        X = atomic_coordinates[atom_X - 1]
        C = atomic_coordinates[atom_C - 1]
        Y = atomic_coordinates[atom_Y - 1]

        CX = X - C
        CY = Y - C
        CX_magnitude = np.linalg.norm(CX)
        CY_magnitude = np.linalg.norm(CY)
        p = np.cross(CX, CY)
        p_magnitude = np.linalg.norm(p)

        dthetadX = np.cross(CX, p) / (CX_magnitude**2 * p_magnitude)
        dthetadY = -1 * np.cross(CY, p) / (CY_magnitude**2 * p_magnitude)
        dthetadC = -dthetadX - dthetadY

        g_angles[atom_X - 1] += dVdtheta * dthetadX
        g_angles[atom_C - 1] += dVdtheta * dthetadC
        g_angles[atom_Y - 1] += dVdtheta * dthetadY

    return g_angles

def calculate_torsional_gradients(nat, nbonds, nC, dihedrals, atoms, atomic_coordinates):
#to obtain the gradient  of energy of the four atoms 
    g_torsional = np.zeros((nat, 3)) # array to store 

    # starting position for where to start placing torsional gradients in Wilson B matrix
    p = nbonds + 6*nC # 

    for dihedral in dihedrals:

        phi = dihedral[4] #torsional angle (phi) degrees 
        phi = np.radians(phi)  #convert to radians 

        # gradient of energy in respect to torsional angle phi
        dVdphi = -1 * n_phi * A_phi * np.sin(n_phi * phi) 

        # for a dihedral X-Y-Z-Q // to extract 
        atom_X, atom_Y = int(dihedral[0]), int(dihedral[1])
        atom_Z, atom_Q = int(dihedral[2]), int(dihedral[3])
        #get cartisian coordinates 
        X = atomic_coordinates[atom_X - 1]
        Y = atomic_coordinates[atom_Y - 1]
        Z = atomic_coordinates[atom_Z - 1]
        Q = atomic_coordinates[atom_Q - 1]
        #calculate vectors 
        XY = Y - X
        YZ = Z - Y
        ZQ = Q - Z
        XZ = Z - X
        YQ = Q - Y
        # orthogonal vectors 
        t = np.cross(XY, YZ)
        u = np.cross(YZ, ZQ)
        # compute the magnitudes 
        t_magnitude = np.linalg.norm(t)
        u_magnitude = np.linalg.norm(u)
        YZ_magnitude = np.linalg.norm(YZ)
        # calculate components of t and u for gradient computation 
        t_component = np.cross(t, YZ) / (t_magnitude**2 * YZ_magnitude)
        u_component = np.cross(-1 * u, YZ) / (u_magnitude**2 * YZ_magnitude)
        # compute partial derivatives with respect to each atom's position
        dphidX = np.cross(t_component, YZ)
        dphidY = np.cross(XZ, t_component) + np.cross(u_component, ZQ)
        dphidZ = np.cross(t_component, XY) + np.cross(YQ, u_component)
        dphidQ = np.cross(u_component, YZ)
        #update torsional gradiente for the four atoms 
        g_torsional[atom_X - 1] += dVdphi*dphidX
        g_torsional[atom_Y - 1] += dVdphi*dphidY
        g_torsional[atom_Z - 1] += dVdphi*dphidZ
        g_torsional[atom_Q - 1] += dVdphi*dphidQ

    return g_torsional

def calculate_vdW_gradients(nat, vdW_pairs, atoms, atomic_coordinates):
# array to store Van der Waals gradients for all atoms
    g_vdW = np.zeros((nat, 3)) 

    for pair in vdW_pairs:
        
        pair_type = int(pair[2]) #extract the Van der Waals pair 
        A = vdW_A.get(pair_type)  # repulsive parameter
        B = vdW_B.get(pair_type)   # attractive parameter 

        i, j = int(pair[0]), int(pair[1]) # Extract atom indices for the VDW pair 
        i_coords = atomic_coordinates[i - 1] # coordinates for the two atoms 
        j_coords = atomic_coordinates[j - 1]
        ji = i_coords - j_coords #compute the displacement vector 
        r_ij = np.linalg.norm(ji) #distance 

        g_i = ji * (-12 * (A / r_ij**14) + 6 * (B / r_ij**8)) # this is eq.40
        g_j = -1 * g_i # same but opposite sign 

        g_vdW[i - 1] += g_i #acumulate the gradient contribution for atom  i 
        g_vdW[j - 1] += g_j  # " " atom j 
    
    return g_vdW
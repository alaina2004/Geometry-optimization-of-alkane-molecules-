import numpy as np

#This file is just to store the values:
# force constant, Van der Waals ... 

# kcal mol-1 Angstrom-2
# 0 = CH, 1 = CC
bond_force_constants = {
    0 : 300,
    1 : 350
}

# 
equilibrium_lenghts = {
    0 : 1.53,
    1 : 1.11
}

# given in kcal mol-1 rad-2 HCC or HCH 35; CCC 60;
# converted to kcal mol-1 deg-2
conv = (1/(180/np.pi))**2
angle_force_constants = {
    0 : 35 * conv,
    1 : 60 * conv
}

# degree
theta_0 = 109.50 

# X-C-C-X dihedral constants
A_phi = 0.3 # kcal mol-1
n_phi = 3 # dimensionless

### VDW PARAMETERS ###
# Since all possible atom pairs in an alkane are CH, CC and HH
# epsilon_ij and sigma_ij are computed only once and saved here
# in the form of A_ij and B_ij (A_ij = 4*epsilon_ij*sigma_ij^12 and B_ij = 4*epsilon_ij*sigma_ij*6)
# CC = 0, CH = 1, HH = 2
# A_ij is in kcal mol-1 Angstrom^12
# B_ij is in kcal mol-1 Angstrom^6
vdW_A = {
    0 : 946181.74,
    1 : 64393.99,
    2 : 4382.44
}
vdW_B = {
    0 : 514.714,
    1 : 108.644,
    2 : 22.932
}



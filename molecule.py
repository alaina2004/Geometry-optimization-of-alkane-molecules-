import numpy as np
import time

from file_IO_functions import *
from structure_functions import *
from gradient_functions import *


# definition of the class Molecule
class Molecule:
    # method constructor of the class 
    def __init__(self, filepath):
# to define all the attributes of my class MOLECULE // alkane 

        self.nat, self.nbonds, self.nC, self.nCC, self.atoms, self.atomic_coordinates, self.bonds = read_mol2(filepath)

    def determine_internal_coordinates(self):
        # this will call find_angles(), "" dihedrals, "" vdw_pairs() from the structure_fuctions
        self.angles, self.C_bonds_dict = find_angles(self.bonds, self.atoms, self.nC)   #the output give me 
        self.dihedrals = find_dihedrals(self.bonds, self.C_bonds_dict, self.nCC)
        # not an internal coordinate but it is practical to do this here
        if self.nat > 5: # condition is only for methane
            self.vdW_pairs = find_vdW_pairs(self.nat, self.bonds, self.angles, self.atoms)
           
    def calculate_internal_coordinates_values(self):
        #this will call calculate_bond_lengths(), "" angles and dihedral   // no need for output because its only functions
        calculate_bond_lenghts(self.bonds, self.atomic_coordinates)
        calculate_angles(self.angles, self.atomic_coordinates)
        calculate_dihedral_angles(self.dihedrals, self.atomic_coordinates)  
    
    def calculate_molecular_energy(self):
        # calls all functions from energies calculating 
        calculate_bond_stretching_energies(self.bonds)
        self.total_bse = np.sum(self.bonds[:, 4])  #bse bond stretching energy // sum of the forth python column 
        calculate_angle_bending_energies(self.angles) # calling the function 
        self.total_abe = np.sum(self.angles[:, 5])  #total bending energy fifth python column 
        calculate_torsion_energies(self.dihedrals)  # "" 
        self.total_te = np.sum(self.dihedrals[:, 5]) #total torsion energy fifth column 
        if hasattr(self, 'vdW_pairs'): # skips this if methane // if it has the attribute or not 
            calculate_vdW_energies(self.vdW_pairs, self.atomic_coordinates)
            self.total_vdW = np.sum(self.vdW_pairs[:,3])
            self.total_e = self.total_bse + self.total_abe + self.total_te + self.total_vdW  # finally total energy  
        else:
            self.total_e = self.total_bse + self.total_abe + self.total_te 

    
    def calculate_gradients(self):
        #same but call the gradient functions
        self.g_bonds = calculate_bond_stretching_gradients(self.nat, self.bonds, self.atoms, self.atomic_coordinates)
        self.g_angles = calculate_angle_bending_gradients(self.nat, self.nbonds, self.angles, self.atoms, self.atomic_coordinates)
        self.g_dihedrals = calculate_torsional_gradients(self.nat, self.nbonds, self.nC, self.dihedrals, self.atoms, self.atomic_coordinates)
        if hasattr(self, 'vdW_pairs'): # skips this if methane
            self.g_vdW = calculate_vdW_gradients(self.nat, self.vdW_pairs, self.atoms, self.atomic_coordinates)
            self.total_gradient = self.g_bonds + self.g_angles + self.g_dihedrals + self.g_vdW
        else:
            self.total_gradient = self.g_bonds + self.g_angles + self.g_dihedrals

    # Step 6 optimization      
    def cartesian_optimization(self):

        self._print_initial_information() # prints the initial geometry information // just for me to get the whole output file 

        # define the initial guess of the inverse hessian // ALL ELEMENTS 0 EXCEPT DIAGONAL ELEMENTS ARE 1/300
        self.M = np.eye(self.nat*3) / 300 # creates a unit matrix of nat*3 dimension and divides it by 300 

        g_RMS = 1.0 # set initial RMS of the gradient (no need to calculate it before the loop)
        g_RMS_threshold = 0.001 # suggest value //  manual 
        c1 = 0.1 # Wolfe rule parameter 
        
        # prepare appropiate dimension for matrix multiplication 
        self.gradient_vector = np.reshape(self.total_gradient, self.nat*3) # reshapes it from nat by 3 matrix to a vector with nat*3 elements

        k = 1 # Optimization cycle counter
        while g_RMS > g_RMS_threshold: # This will loop until a structure with g_RMS <= g_RMS_threshold is found

            print("########################################")
            print(f"# Geometry optimization cycle number {k} #")
            print("########################################")
            
            # Predicted geometry change (Eq 9)
            p_k = -1 * np.matmul(self.M, self.gradient_vector)  # indicates me which direction to change the geometry 

            # setting the parameters and variables // each optimization new V and g 
            alpha = 0.8 # Wolfe rule parameter 
            V_old = self.total_e
            g_old = self.gradient_vector
            initial_coordinates = np.copy(self.atomic_coordinates)

            # Rough line search  // to make optimal the step towards that direction 
            while True: # This will loop until s_k = alpha*p_k is found that satisfies Wolfe's first rule (Eq 8)
                
                # set s_k for this alpha test
                s_k = alpha * p_k
                
                # reset the coordinates for this alpha test
                self.atomic_coordinates = np.reshape(initial_coordinates, self.nat*3)

                # calculate possible new coordinates with a current choice of step size
                self.atomic_coordinates = self.atomic_coordinates + s_k
                self.atomic_coordinates = np.reshape(self.atomic_coordinates, (self.nat, 3))
                
                # calculate energetic and gradient terms for the geometry obtained by this alpha test
                self.calculate_internal_coordinates_values()
                self.calculate_molecular_energy()
                self.calculate_gradients()
                self.gradient_vector = np.reshape(self.total_gradient, self.nat*3)

                # Wolfe first rule
                lhs = self.total_e 
                rhs = V_old + c1 * alpha * np.dot(p_k, g_old)
                
                # exits the loop if Wolfe's first rule condition is met
                if lhs <= rhs:
                    break 
                alpha = alpha * 0.8 # otherwise, decrease alpha by a factor of 0.8 
            
            print(f"New structure r_k+1:")
            for i, row in enumerate(self.atomic_coordinates):
                print(f"{self.atoms[i]} {row[0]:>10.6f} {row[1]:>10.6f} {row[2]:>10.6f}")
            
            # Calculate the root mean square of the gradient
            g_RMS = np.sqrt(np.sum(self.gradient_vector**2) / len(self.gradient_vector))
            # if it is less than 0.001 the optimization will stop here

            print(f"Potential energy in the previous cycle {V_old:.6f} kcal/mol")
            print(f"Potential energy at this cycle's geometry {self.total_e:.6f} kcal/mol")
            print(f"RMS of the gradient at this cycle's geometry {g_RMS:.6f} kcal/molA")
            

            # Update of the guess of the inverse Hessian matrix M (Eq 12)
            y_k = self.gradient_vector - g_old
            v_k = np.matmul(self.M, y_k)            
            first_term = ((np.dot(s_k, y_k) + np.dot(y_k, v_k)) * np.outer(s_k, s_k))/(np.dot(s_k, y_k))**2
            second_term = (np.outer(v_k, s_k) + np.outer(s_k, v_k)) / np.dot(s_k, y_k)
            self.M += first_term - second_term

            k += 1 # updates the cycle counter after 1 whole cycle is completed

    def _print_initial_information(self):
            
        print(f"Molecule has been initialized with the following properties:")
        print(f"- Molecular formula: C{self.nC if self.nC > 1 else ""}H{self.nat - self.nC}; Total number of atoms: {self.nat}")
        print(f"- Initial atomic coordinates (Angstrom)")
        for i, atom in enumerate(self.atoms):
            print(f"{atom}  {self.atomic_coordinates[i][0]:>10.6f}   {self.atomic_coordinates[i][1]:>10.6f}   {self.atomic_coordinates[i][2]:>10.6f}")
        print(f"- This molecule is defined by:\n -- {self.nat*3} Cartesian coordinates")
        print(f"- Internal coordinates: {len(self.bonds)} stretching, {len(self.angles)} bending and {len(self.dihedrals)} torsions.")
        print(f"Potential energy at the initial structure is {self.total_e:.6f} kcal / mol")
        print(f"Stretching, bending, torsional and VdW components of the potential energy at the initial structure are: ")
        print(f" - {self.total_bse:.6f} {self.total_abe:.6f} {self.total_te:.6f} {self.total_vdW if hasattr(self, 'total_vdW') else 0.000000}")
import sys
from molecule import *

def main():
    if len(sys.argv) != 2:
        print("Usage: python main.py <filename.mol2>")
        sys.exit(1)

    filename = sys.argv[1]
    print("Welcome to a small alkane structure optimizer program!")
    print(f"Reading molecule from file: {filename}")
    print("Optimization will begin now!")

    # Initialize a molecule from the input file
    molecule = Molecule(filename)

    # Start the optimization in Cartesian coordinates
    print("Starting Cartesian optimization...")
    # INITIAL GEOMETRY
    molecule.determine_internal_coordinates()
    molecule.calculate_internal_coordinates_values()
    molecule.calculate_molecular_energy()
    molecule.calculate_gradients()
    # OPTIMIZATION
    molecule.cartesian_optimization()
    print("Optimization completed!")

if __name__ == "__main__":
    main()

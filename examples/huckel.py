from tcintegral import basis_set, MolecularOrbital

# get a basis set from basis_sets
bs = basis_set.STO2G

# read a molecule to use
mol = plams.Molecule('/Users/yumanhordijk/Desktop/acrolein.xyz')
mol.guess_bonds()

# get atomic orbitals for this molecule
# orb_to_get = ['1S', '2S', '3S', '1P:x', '1P:y', '1P:z']  # list of possible AO
orb_to_get = ['1P:z']  # list of possible AOs
aos = []
for atom in mol:
    for orb_name in orb_to_get:
        try:  # try to get the AO for this atom, it will fail if it does not exist
            aos.append(bs.get(atom.symbol, orb_name, atom.coords))
        except IndexError:
            pass

# build the Hamiltonian
ionization_potentials = {  # in Hartrees, obtained at the OLYP/QZ4P level in ADF
    'H(1S)': -0.238455,
    'C(1S)': -10.063764,
    'C(2S)': -0.506704,
    'C(1P:x)': -0.188651,
    'C(1P:y)': -0.188651,
    'C(1P:z)': -0.188651,
    'O(1S)': -18.9254,
    'O(2S)': -0.886887,
    'O(1P:x)': -0.329157,
    'O(1P:y)': -0.329157,
    'O(1P:z)': -0.329157,
}

K = 1.75
H = np.zeros((len(aos), len(aos)))
# build diagonal elements
for i, ao in enumerate(aos):
    H[i, i] = ionization_potentials.get(ao.name)

# build off-diagonal elements
for i, ao1 in enumerate(aos):
    for j, ao2 in enumerate(aos[i+1:], start=i+1):
        H[i, j] = H[j, i] = K * ao1.overlap(ao2) * (H[i, i] + H[j, j]) / 2

# get MO energies and coeffs
energies, coefficients = np.linalg.eigh(H)

# we are now going to draw each MO
for mo_index in range(len(energies)):
    # construct the MO using our AOs and coefficients
    mo = MolecularOrbital(aos, coefficients[:, mo_index], mol)

    # create a new screen
    with Screen() as scr:
        scr.draw_molecule(mol)
        scr.draw_orbital(mo, material=materials.orbital_shiny)
        scr.draw_axes()

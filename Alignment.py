from Bio.PDB import PDBParser, Superimposer, Selection, NeighborSearch
import numpy as np

chain_dict = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z', 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z', 'AA', 'BB', 'CC', 'DD', 'EE', 'FF', 'GG', 'HH', 'II', 'JJ', 'KK', 'LL', 'MM', 'NN', 'OO', 'PP', 'QQ', 'RR', 'SS', 'TT', 'UU', 'VV', 'WW', 'XX', 'YY', 'ZZ', 'aa', 'bb', 'cc', 'dd', 'ee', 'ff', 'gg', 'hh', 'ii', 'jj', 'kk', 'll', 'mm', 'nn', 'oo', 'pp', 'qq', 'rr', 'ss', 'tt', 'uu', 'vv', 'ww', 'xx', 'yy', 'zz']

def next_available_chain_label(model, excluded):
    letters_used = [x._id for x in model.get_chains()]
    for letter in chain_dict:
        if letter not in letters_used and letter not in excluded:
            return letter
    return None

def do_superimpose(fixed_model, moveable_complex, fixed_chain, target_chain, moveable_chain):
    parser = PDBParser()
    s2 = parser.get_structure(moveable_complex.id, moveable_complex.filename)

    s1_fixed_chain = fixed_model[fixed_chain.label]
    s2_target_chain = s2[0][target_chain.label]

    s2_moved_chain = s2[0][moveable_chain.label]
    s2_moved_chain_atoms = []

    new_chain_id = next_available_chain_label(fixed_model, set([s2_target_chain.id]))
    s2_moved_chain.id = new_chain_id
    moveable_chain.label = new_chain_id

    # Extract information of atoms
    for residues in s2_moved_chain:
        s2_moved_chain_atoms.append(residues)

    sup = Superimposer()
    # Specify the atom lists
    # 'fixed' and 'moving' are lists of Atom objects
    # The moving atoms will be put on the fixed atoms
    sup.set_atoms(list(s1_fixed_chain.get_atoms()), list(s2_target_chain.get_atoms()))

    # Apply rotation/translation to the moving atoms
    sup.apply(s2_moved_chain_atoms)

    fixed_model.add(s2_moved_chain)

    return s2_moved_chain


def get_clashes(model, chain):
    """

    """

    #Storing coords in arrays
    x_coords = list()
    y_coords = list()
    z_coords = list()
    for residue in chain:
        if residue.has_id("CA"):
            x_coords.append(((residue["CA"].get_coord())[0]))
            y_coords.append(((residue["CA"].get_coord())[1]))
            z_coords.append(((residue["CA"].get_coord())[2]))
        
    #Center coordinates
    center_coords = np.array([round(np.mean(x_coords),3), round(np.mean(y_coords),3), round(np.mean(z_coords),3)])

    #Max axis distances
    dist_x_center = (max(x_coords)-min(x_coords))/2
    dist_y_center = (max(y_coords)-min(y_coords))/2
    dist_z_center = (max(z_coords)-min(z_coords))/2

    vector = np.array([dist_x_center, dist_y_center, dist_z_center])

    #radius_max_coord = (max(vector))
    radius_module = (np.linalg.norm(vector))

    for env in model:
        if env != chain:
            atoms  = Selection.unfold_entities(env, 'A')
            ns = NeighborSearch(atoms)

            close_atoms_to_center = ns.search(center_coords, radius_module)
            
            close_atoms_to_residue = None
            if close_atoms_to_center:
                for residue in chain:
                    if residue.has_id("CA"):
                        coord_res = residue["CA"].get_coord()
                        close_atoms_to_residue = ns.search(coord_res, 5)
                        if close_atoms_to_residue:
                            print("Evaluating possible clash...")
                            for atom in residue:
                                coord_atom = (atom.get_coord())
                                close_atoms_to_atom = ns.search(coord_atom, 1.1) #CH distance 1.09
                                if close_atoms_to_atom:
                                    print("WARNING IMPORTANT CLASH - CHANGING CHAIN")
                                    return(True)
                            
            return(False)
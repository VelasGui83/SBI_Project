from Bio.PDB import PDBParser, Superimposer, Selection, NeighborSearch
import numpy as np
import alphabet

def next_available_chain_label(model, excluded):
    letters_used = [x._id for x in model.get_chains()]
    for letter in alphabet.chain_alphabet:
        if letter not in letters_used and letter not in excluded:
            return letter
    return None

def do_superimpose(fixed_model, moveable_complex, fixed_chain, target_chain, moveable_chain):
    parser = PDBParser()
    s2 = parser.get_structure(moveable_complex.id, moveable_complex.filename)

    s1_fixed_chain = fixed_model[fixed_chain]
    fixed_atoms = []

    if s1_fixed_chain.child_list[0].has_id("CA"):
        for fixed_res in s1_fixed_chain:
            if fixed_res.has_id("CA"):
                fixed_atoms.append(fixed_res['CA'])
    else:
        fixed_atoms = list(s1_fixed_chain.get_atoms())


    s2_target_chain = s2[0][target_chain]
    target_atoms = []
    if s2_target_chain.child_list[0].has_id("CA"):
        for target_res in s1_fixed_chain:
            if target_res.has_id("CA"):
                target_atoms.append(target_res['CA'])
    else: 
        target_atoms = list(s2_target_chain.get_atoms())


    sup = Superimposer()
    # Specify the atom lists
    # 'fixed' and 'moving' are lists of Atom objects
    # The moving atoms will be put on the fixed atoms
    sup.set_atoms(fixed_atoms, target_atoms)


    # Apply rotation/translation to the moving atoms
    s2_moved_chain = s2[0][moveable_chain]
    s2_moved_chain_atoms = []

    new_chain_id = next_available_chain_label(fixed_model, set([s2_target_chain.id]))
    s2_moved_chain.id = new_chain_id

    # Extract information of atoms
    for residues in s2_moved_chain:
        s2_moved_chain_atoms.append(residues)

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

def get_clashes_v2(model, chain):
    """

    """

    atoms_to_compare = []
    for chain_moved in model:
        if chain != chain_moved:
            atoms_to_compare = atoms_to_compare + list(chain_moved.get_atoms())
    
    atoms  = Selection.unfold_entities(atoms_to_compare, 'A')
    ns = NeighborSearch(atoms)
    print("Evaluating possible clashes")
    print("Comparing to ",len(list(model.get_atoms())), " entities")
    for residue in chain:
        for atom in residue:
            coord_atom = (atom.get_coord())
            close_atoms_to_atom = ns.search(coord_atom, 1.1) #CH distance 1.09
            if close_atoms_to_atom:
                print("WARNING IMPORTANT CLASH - CHANGING CHAIN")
                return(True)
    
    return(False)
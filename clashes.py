import Bio.PDB

parser = Bio.PDB.PDBParser(QUIET=True) 
import numpy as np


def get_clashes(structure, chain):
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

    for env in structure[0]:
        if env != chain:
            atoms  = Bio.PDB.Selection.unfold_entities(env, 'A')
            ns = Bio.PDB.NeighborSearch(atoms)

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
                                close_atoms_to_atom = ns.search(coord_atom, 1.5) #CH distance 1.09
                                if close_atoms_to_atom:
                                    print("WARNING IMPORTANT CLASH - CHANGING CHAIN")
                                    return(True)
                            
            return(False)

                            

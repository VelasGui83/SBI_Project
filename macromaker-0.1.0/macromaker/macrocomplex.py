import copy
import pickle
import os

import numpy as np
from Bio.PDB import NeighborSearch, PDBParser, Selection, Superimposer
from Bio.PDB.Polypeptide import PPBuilder

import macromaker.alphabet as alpha
import macromaker.lib.center_of_mass as com
import macromaker.interactions as inter

class Macrocomplex(object):
    """Class to model a macrocomplex"""
    def __init__(self, verbose=False):
        """Constructor of the Macrocomplex class"""
        self.model = None
        self.to_evaluate = None
        self.interactions = inter.Interactions()
        self.verbose = verbose


    def next_available_chain_label(self, excluded=None):
        """Returns the next avaliable chain label in the model.

        Arguments:
            - excluded - set, aditional set of labels to be omited

        """
        letters_used = [x._id for x in self.model.get_chains()]
        for letter in alpha.chain_alphabet:
            if letter not in letters_used and letter not in excluded:
                return letter
        return None

    def get_ca_atoms(self, chain):
        """Returns the CA atoms of a chain, or all the atoms in case it doesn't have CA atoms.

        Arguments:
            - chain - Bio.PDB.Chain, the chain to get the atoms

        """
        result = []
        if chain.child_list[0].has_id("CA"):
            for fixed_res in chain:
                if fixed_res.has_id("CA"):
                    result.append(fixed_res['CA'])
        else:
            result = list(chain.get_atoms())
        return result

    def do_superimpose(self, fixed_chain, target_chain, moved_chain):
        """Do a superimposition of the homologous chains.

        Arguments:
            - fixed_chain - Bio.PDB.Chain, the fixed chain to superimpose
            - target_chain - Bio.PDB.Chain, the target homologue chain
            - moveable_chain - Bio.PDB.Chain, the chain we apply the movement

        """
        # Obtain the atoms of the homologues chains
        fixed_atoms = self.get_ca_atoms(fixed_chain)
        target_atoms = self.get_ca_atoms(target_chain)

        sup = Superimposer()
        # Specify the atom lists
        # 'fixed' and 'target' are lists of Atom objects
        # The target atoms will be put on the fixed atoms
        sup.set_atoms(fixed_atoms, target_atoms)

        # Declare moving chain
        moved_chain_atoms = []
        # Extract information of atoms
        for residues in moved_chain:
            moved_chain_atoms.append(residues)

        # Apply rotation/translation to the moving atoms
        sup.apply(moved_chain_atoms)

    def is_chain_clashed_v2(self, chain, limit=5):
        """Retuns true if there are clashes of the chain in the model.

        Arguments:
            - chain - Bio.PDB.Chain, a chain of the model
            - limit - int, the treshold of clashes

        """
        clashes = 0
        atoms_to_compare = list(self.model.get_atoms())
        atoms  = Selection.unfold_entities(atoms_to_compare, 'A')
        ns = NeighborSearch(atoms)

        print("Evaluating possible clashes")
        print("Comparing to ",len(list(self.model.get_atoms())), " entities")
        for residue in chain:
            for atom in residue:
                coord_atom = (atom.get_coord())
                close_atoms_to_atom = ns.search(coord_atom, 1.1) #CH distance 1.09
                if close_atoms_to_atom:
                    clashes += 1
                    if clashes >= limit:
                        print("WARNING IMPORTANT CLASH - CHANGING CHAIN")
                        return(True)
        
        return(False)

    def coord_atom_distance(self, coord, atom): 
        """Returns the distance between a coord and an atom.

        Arguments:
            - coord - np.Array, a np array of thre floats 
            - atom - Bio.PDB.Atom, the output dir where the file is created

        """
        return np.linalg.norm(atom.coord - coord)


    def is_chain_clashed_v3(self, chain, limit=5):
        """Retuns true if there are clashes of the chain in the model.

        Arguments:
            - chain - Bio.PDB.Chain, a chain of the model
            - limit - int, the treshold of clashes

        """
        clashes = 0
        atoms_to_compare = list(self.model.get_atoms())
        atoms  = Selection.unfold_entities(atoms_to_compare, 'A')
        ns = NeighborSearch(atoms)
        
        max_radius = 0
        center = np.array(com.center_of_mass_chain(chain, 'ATOM'), dtype=float)
        for atom in chain.get_atoms():
            if self.coord_atom_distance(center, atom) > max_radius:
                max_radius = self.coord_atom_distance(center, atom)

        close_atoms_to_chain = ns.search(center, max_radius+2)
        
        if close_atoms_to_chain:
            ns = NeighborSearch(close_atoms_to_chain)
            if self.verbose:
                print("Evaluating possible clashes")
                print("Comparing to ",len(close_atoms_to_chain), " entities")
            for residue in chain:
                for atom in residue:
                    coord_atom = (atom.get_coord())
                    close_atoms_to_atom = ns.search(coord_atom, 0.8) #CH distance 1.09
                    if close_atoms_to_atom:
                        clashes += 1
                        if clashes >= limit:
                            if self.verbose:
                                print("WARNING IMPORTANT CLASH - CHANGING CHAIN")
                            return True 
        return False

    def __add_chains_to_model(self, max_chains):
        """Add chains to a model and return the number of final chains.

        Arguments:
            - to_evaluate - list of interaction.Chain, a list of populated chians to evaluate
                            for our model
            - max_chains - int, maximum number of chains of the model

        """
        while self.to_evaluate:
            evaluating_chain = self.to_evaluate.pop(0)

            for candidate_chain in evaluating_chain.homologous_chains:
                candidate_complex = candidate_chain.parent
                complementary_candidate_chain = candidate_complex.complementary_chain(candidate_chain)
                
                parser = PDBParser()
                candidate_model = parser.get_structure(candidate_complex.id, candidate_complex.filename)

                fixed_chain = self.model[evaluating_chain.label]
                target_chain = candidate_model[0][candidate_chain.original_label]
                moved_chain = candidate_model[0][complementary_candidate_chain.original_label]
                self.do_superimpose(fixed_chain, target_chain, moved_chain)

                if not self.is_chain_clashed_v3(moved_chain, limit=10):
                    moved_chain.id = self.next_available_chain_label(set([target_chain.id]))
                    self.model.add(moved_chain)

                    next_to_evaluate = copy.copy(complementary_candidate_chain)
                    next_to_evaluate.label = moved_chain.id
                    self.to_evaluate.append(next_to_evaluate)

                    if self.verbose:
                        print("Chain added. #%d" %len(self.model))

                    if len(self.model) >= max_chains:
                        return len(self.model)
        return len(self.model)

    def create_macrocomplex(self, input_folder, max_chains=200):
        """Retruns a Bio.PDB.Structure composed by complex pairs.

        Arguments:
            - input_folder - string, folder containing a list of compex pairs
            - max_chains - int, maximum number of chains

        """
        self.interactions.populate_interactions(input_folder)
        self.interactions.populate_homologous_chains(score_limit=9.8)

        remaining_complexes = self.interactions.get_complexes_list()

        parser = PDBParser()
        chain_number = 2

        if self.verbose:
            print("Starting the build...")

        while chain_number == 2:
            first_complex = self.interactions.interactions[remaining_complexes.pop(0)]
            first_chains = first_complex.get_chain_list()
            self.to_evaluate = [first_chains[0], first_chains[1]]
            structure = parser.get_structure("final_structure", first_complex.filename)
            self.model = structure[0]

            chain_number = self.__add_chains_to_model(max_chains)
        if self.verbose:
            print("Model ended.")
        return structure

    def continue_macrocomplex(self, max_chains = 200):
        """Continues with the execution of the macrocomplex.

        Arguments:
            - max_chains - int, maximum number of chains

        """
        if self.verbose:
            print("Continuing the rescued model.")
        self.__add_chains_to_model(max_chains)
        return self.model.get_parent()

    def save_context(self, input_file):
        """Save the context of the execution
        
        Arguments:
            - input_file - string, file to save the state
        
        """
        out_fd = open(input_file, "wb")
        data = (self.model.get_parent(), self.interactions, self.to_evaluate, self.verbose)
        pickle.dump(data, out_fd)
        out_fd.close()

    def load_context(self, input_file):
        """Loads the context pf the execution
        
        Arguments:
            - input_file - string, file to load the state
        
        """
        data = pickle.load(open(input_file, "rb"))
        self.model = data[0][0]
        self.interactions = data[1]
        self.to_evaluate = data[2]
        self.verbose = data[3]
        os.system("rm context.p")
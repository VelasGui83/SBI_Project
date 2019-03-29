import copy

import Levenshtein
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import PPBuilder

import stamp_adapter as stamp
import utils


class Chain(object):
    """Class to save the homologous chains of the chain"""
    def __init__(self, label, parent):
        """Creator of Chain class

        Arguments:
         - label - string, the label of the chain
         - parent - Complex, the complex instance of the parent

        """        
        self.label = label
        self.original_label = label
        self.parent = parent
        self.homologous_chains = set([])

    def __len__(self):
        """Returns the number of homologous chains"""
        return len(self.homologous_chains)

    def __get_nucleotides(self, chain):
        """Returns the nucleotide sequence of a chain

        Arguments:
         - chain - Bio.PDB.Chain, the chain containing nucleotides

        """ 
        return "".join([x.resname.strip() for x in chain.child_list])

    def __get_homologous_tuples(self, interactions, identity=1.0):
        """Returns a tuple with the type of chain and a list of its homologous chains,
        with an identity treshold

        Arguments:
         - interactions - dict, the dict conteining all the interactions
         - identity - float, the limit filter of identity of sequences

        """ 
        ppb = PPBuilder()

        tuple_list = [(self.parent.filename, self.parent.id, self.original_label)]
        chain_type = "prot"

        evaluating_structure = PDBParser().get_structure(self.original_label, self.parent.filename)
        evaluating_polypeptide = ppb.build_peptides(evaluating_structure[0][self.original_label])
        evaluating_sequence = None
        if evaluating_polypeptide:       
            evaluating_sequence = evaluating_polypeptide[0].get_sequence()
        else:
            chain_type = "dna/rna"
            evaluating_sequence = self.__get_nucleotides(evaluating_structure[0][self.original_label])

        for complex_item in interactions.keys():
            complex_structure = PDBParser().get_structure(complex_item, interactions[complex_item].filename)
            for label in interactions[complex_item].chain_dict.keys():
                if complex_item != self.parent.id or label != self.label:
                    chain_polypeptide = ppb.build_peptides(complex_structure[0][label])
                    chain_sequence = None
                    if chain_polypeptide:
                        chain_sequence = chain_polypeptide[0].get_sequence()
                    else:
                        chain_sequence = self.__get_nucleotides(complex_structure[0][label])
                    if Levenshtein.ratio(str(evaluating_sequence), str(chain_sequence)) >= identity:
                        tuple_list.append((interactions[complex_item].filename, complex_item, label))
        return (chain_type, tuple_list)

    def get_homologous_chains(self, interactions, score_limit = 9.8):
        """Returns a tuple with the type of chain and a list of its homologous chains,
        with an identity treshold

        Arguments:
         - interactions - dict, the dict conteining all the interactions
         - score_limit - float, the limit for score similarity

        """ 
        homologous_set = set([])
        
        homologous_tuples = self.__get_homologous_tuples(interactions, 0.95)
        if len(homologous_tuples[1]) != 1 and homologous_tuples[0] == "prot":
            #With STAMP
            sp = stamp.STAMPParser()
            sp.create_stamp_input(homologous_tuples[1], "atomator/tmp/")
            stamp_id = self.parent.id+"_"+self.original_label
            scores_dict = sp.get_stamp_scores(input_folder="atomator/tmp/", chain_id=stamp_id, limit=score_limit)
            if len(scores_dict) != 0:
                homologous_list = [x[0] for x in scores_dict[stamp_id]]
                for homologous_chain in homologous_list:
                    complex_id = homologous_chain.split("_")[0]
                    chain_label = homologous_chain.split("_")[1]
                    homologous_set.add(interactions[complex_id].chain_dict[chain_label])
        else:
            # #Without STAMP
            homologous_list = self.__get_homologous_tuples(interactions, 1)
            homologous_list[1].pop(0)
            for homologous_chain in homologous_list[1]:
                complex_id = homologous_chain[1]
                chain_label = homologous_chain[2]
                homologous_set.add(interactions[complex_id].chain_dict[chain_label])
        return homologous_set

class Complex(object):
    """Class to save the complex file"""
    def __init__(self, id, filename):
        """Creator of Complex class

        Arguments:
         - id - string, the id of the complex
         - filename - string, the file name of the complex

        """
        self.id = id
        self.chain_dict = {}
        self.filename = filename

    def get_chain_list(self):
        """Return a list of the child chains"""
        return list(self.chain_dict.values())

    def complementary_chain(self, chain):
        """Return the complementary chain"""
        # We supose that every complex has two chain childs
        for chain_item in self.chain_dict.values():
            if chain_item != chain:
                return chain_item
        return None

class Interactions(object):
    """Class to save the interactions between chains"""

    def __init__(self):
        """Creator of Interactions class"""
        self.interactions = {}
        
    def get_complexes_list(self):
        """Returns a list of the interactions complexes"""
        return list(self.interactions.keys())

    def __get_chains(self, input_file):
        """Returns a set with all the chains found in a PDB file

        Arguments:
         - input_file - string, the PDB file to look for chains

        """
        set_chains = set([])
        with open(input_file, 'r') as fh:
            for line in fh:
                if line.startswith("ATOM"):
                    set_chains.add(line[21])
        return set_chains

    def populate_interactions(self, input_folder):
        """Populate interactions dictionary with empty homologues

        Arguments:
         - input_folder - string, folder to look at for PDB files

        """
        input_files = utils.get_files(input_path=input_folder, admited_formats=set(["pdb"]))

        i = 1
        for filename in input_files:
            chains = sorted(list(self.__get_chains(filename)))
            for chain in chains:
                complex_id = "%d:%d" %(i, i+1)
                try:
                    self.interactions[complex_id].chain_dict[chain] = Chain(chain, self.interactions[complex_id])
                except KeyError:
                    self.interactions[complex_id] = Complex(complex_id, filename)
                    self.interactions[complex_id].chain_dict[chain] = Chain(chain, self.interactions[complex_id])
            i += 2

    def populate_homologous_chains(self, score_limit):
        """Populate the homologous chains from the interactions dict"""
        for complex_item in self.interactions.values():
            for chain_item in complex_item.chain_dict.values():
                if not chain_item.homologous_chains:
                    homologous_set = chain_item.get_homologous_chains(self.interactions, score_limit)
                    homologous_set.add(chain_item)
                    for chain in homologous_set:
                        homologous_specific_set = copy.copy(homologous_set)
                        homologous_specific_set.remove(chain)
                        chain.homologous_chains = homologous_specific_set

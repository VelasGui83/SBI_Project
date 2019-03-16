from random import randint
from Bio.PDB import PDBParser
import subprocess

from Bio import BiopythonWarning

import warnings
warnings.simplefilter("ignore", BiopythonWarning)

class STAMPChain(object):
    """

    """
    def __init__(self, label, parent):
        """

        """
        self.label = label
        self.parent = parent
        self.similar_chains = []

    def __len__(self):
        """

        """
        return len(self.similar_chains)

    def filter_similar_chains(self, limit):
        """

        """
        self.similar_chains = list(filter(lambda x: x[1] >= limit, self.similar_chains))

    def order_similar_chains(self):
        """

        """
        self.similar_chains.sort(key=lambda x: x[1], reverse=True)


class STAMPComplex(object):
    """

    """
    def __init__(self, id, filename):
        """

        """
        self.id = id
        self.chain_dict = {}
        self.filename = filename

    def get_chain_list(self):
        """

        """
        return list(self.chain_dict.values())

    def complementary_chain(self, chain):
        """

        """
        for chain_item in self.chain_dict.values():
            if chain_item != chain:
                return chain_item
        return None
    
    def filter_chain_scores(self, limit):
        """

        """
        for chain in self.chain_dict.values():
            chain.filter_similar_chains(limit)
            chain.order_similar_chains()

class STAMPParser(object):
    """

    """
    def __init__(self):
        self.id = randint(0,10000)
        self.comlex_dict = {}
        
    def get_complexes_list(self):
        """

        """
        return list(self.comlex_dict.keys())

    def create_stamp_input(self, input_files):
        """

        """
        fh_write = open("tmp/output_"+str(self.id)+".domain",'w')

        parser = PDBParser()
        i = 1
        for filename in input_files:
            # TODO: change read PDB for regexp
            structure = parser.get_structure(i, filename)
            for chain in structure.get_chains():
                complex_id = "%d:%d" %(i, i+1)

                fh_write.write("../%s %d:%d_%s {CHAIN %s}\n" %(filename, i, i+1, chain.id, chain.id))

                try:
                    self.comlex_dict[complex_id].chain_dict[chain.id] = STAMPChain(chain.id, self.comlex_dict[complex_id])
                except KeyError:
                    self.comlex_dict[complex_id] = STAMPComplex(complex_id, filename)
                    self.comlex_dict[complex_id].chain_dict[chain.id] = STAMPChain(chain.id, self.comlex_dict[complex_id])
            i += 2
        fh_write.write("\n")
        fh_write.close()

    def calculate_stamp_scores(self):
        """

        """
        cmd = "stamp -l output_%d.domain -MAX_SEQ_LEN 100000 -rough -n 2 -prefix output_%d" %(self.id, self.id)
        stamp_process = subprocess.Popen(cmd.split(), cwd="tmp/", stdout=subprocess.PIPE)

        for line in stamp_process.stdout:
            line = line.decode("utf-8")
            if line.startswith("Pair"):
                splited_line = line.split()

                if splited_line[6] == splited_line[7]:
                    chain = self.comlex_dict[splited_line[2].split("_")[0]].chain_dict[splited_line[2].split("_")[1]]
                    target_chain = self.comlex_dict[splited_line[3].split("_")[0]].chain_dict[splited_line[3].split("_")[1]]
                    chain.similar_chains.append((target_chain, float(splited_line[4])))

        cmd = "rm output_%d.*" %self.id
        subprocess.Popen(cmd, cwd="tmp/", shell=True)
        cmd = "rm stamp_rough.trans"
        subprocess.Popen(cmd.split(), cwd="tmp/")

    def filter_stamp_scores(self, limit):
        """
        
        """
        for complex_item in self.comlex_dict.values():
            complex_item.filter_chain_scores(limit)
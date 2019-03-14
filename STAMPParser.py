from random import randint
from Bio.PDB import PDBParser
import subprocess

from Bio import BiopythonWarning

import warnings
warnings.simplefilter("ignore", BiopythonWarning)

class STAMPParser(object):
    """

    """
    def __init__(self):
        self.id = randint(0,10000)
        self.stamp_scores = {}
        self.complex_files = {}

    def complementary_chain(self, chain):
        """

        """
        chains = self.complex_files[chain[:-2]]
        for chain_item in chains.keys():
            if chain_item != chain[-1]:
                return chain[:-1]+chain_item
        return None

    def get_complexes_list(self):
        """

        """
        return list(self.complex_files.keys())        

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
                fh_write.write("../%s %d:%d_%s {CHAIN %s}\n" %(filename, i, i+1, chain.id, chain.id))
                try:
                    self.complex_files["%d:%d" %(i, i+1)][chain.id] = filename
                except KeyError:
                    self.complex_files["%d:%d" %(i, i+1)] = {}
                    self.complex_files["%d:%d" %(i, i+1)][chain.id] = filename
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
                try:
                    self.stamp_scores[splited_line[2]].append((splited_line[3], float(splited_line[4])))
                except KeyError:
                    self.stamp_scores[splited_line[2]] = []
                    self.stamp_scores[splited_line[2]].append((splited_line[3], float(splited_line[4])))
        cmd = "rm output_%d.*" %self.id
        subprocess.Popen(cmd, cwd="tmp/", shell=True)
        cmd = "rm stamp_rough.trans"
        subprocess.Popen(cmd.split(), cwd="tmp/")

    def filter_stamp_scores(self, limit):
        """
        
        """
        for key in self.stamp_scores:
            self.stamp_scores[key] = list(filter(lambda x: x[1] >= limit, self.stamp_scores[key]))
            self.stamp_scores[key].sort(key=lambda x: x[1], reverse=True)
            

        





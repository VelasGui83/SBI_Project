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

    def create_stamp_input(self, input_tuples):
        """

        """
        fh_write = open("tmp/output_"+str(self.id)+".domain",'w')

        for record in input_tuples:
            fh_write.write("%s %s_%s {CHAIN %s}\n" %(record[0], record[1], record[2], record[2]))
        fh_write.write("\n")
        fh_write.close()

    def get_stamp_scores(self, chain_id=None, limit=9.8):
        """

        """
        cmd = "stamp -l output_%d.domain -MAX_SEQ_LEN 100000 -rough -n 2 -prefix output_%d" %(self.id, self.id)
        stamp_process = subprocess.Popen(cmd.split(), cwd="tmp/", stdout=subprocess.PIPE)

        result_dict = {}

        for line in stamp_process.stdout:
            line = line.decode("utf-8")
            if line.startswith("Pair"):

                splited_line = line.split()
                first_chain = splited_line[2]
                second_chain = splited_line[3]
                length_chain_1 = splited_line[6]
                length_chain_2 = splited_line[7]
                score = float(splited_line[4])
                if length_chain_1 == length_chain_2 and score >= limit and (chain_id is None or chain_id == first_chain):
                    try:
                        result_dict[first_chain].append((second_chain, score))
                    except KeyError:
                        result_dict[first_chain] = [(second_chain, score)]

        cmd = "rm output_%d.*" %self.id
        subprocess.Popen(cmd, cwd="tmp/", shell=True)
        cmd = "rm stamp_rough.trans"
        subprocess.Popen(cmd.split(), cwd="tmp/")

        return result_dict
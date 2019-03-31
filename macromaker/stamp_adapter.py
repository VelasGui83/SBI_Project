import os.path
import subprocess
import warnings
from random import randint

from Bio import BiopythonWarning
from Bio.PDB import PDBParser

warnings.simplefilter("ignore", BiopythonWarning)

class STAMPParser(object):
    """Class to handle STAMP IO operations"""

    def __init__(self):
        """Creator of STAMPParser class"""
        self.id = randint(0,10000)

    def create_stamp_input(self, input_tuples, output_folder):
        """Creates in the output dir, the input file needed to run STAMP.

        Arguments:
         - input_tuples - list of tuples, the information needed to create the file
         - output_folder - string, the output dir where the file is created

        """
        output_file = os.path.join(output_folder, "output_"+str(self.id)+".domain")
        fh_write = open(output_file,'w')

        for record in input_tuples:
            fh_write.write("%s %s_%s {CHAIN %s}\n" %(record[0], record[1], record[2], record[2]))
        fh_write.write("\n")
        fh_write.close()

    def get_stamp_scores(self, input_folder, chain_id=None, limit=9.8):
        """Return a dicctionary with the chain as a key, and a tuple of the homologous 
        chain and the STAMP score as a value.

        Arguments:
         - input_folder - string, the folder where stamp *.domains is located
         - chain_id - string, the chain selected to get the scores. If none, all chains are selected.
         - limit - float, the treshold for the selection of homologues based on stamp score.

        """
        input_file = os.path.abspath(os.path.join(input_folder, "output_"+str(self.id)+".domain"))
        
        cmd = "stamp -l %s -MAX_SEQ_LEN 100000 -rough -n 2 -prefix output_%d" %(input_file, self.id)
        stamp_process = subprocess.Popen(cmd.split(), cwd=input_folder, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        result_dict = {}

        # TODO: auto-correct wrong pdbs
        # if stamp_process.stderr:
        #     stamp_process = subprocess.Popen(cmd.split(), cwd=input_folder, stdout=subprocess.PIPE)

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
        subprocess.Popen(cmd, cwd=input_folder, shell=True)
        cmd = "rm stamp_rough.trans"
        subprocess.Popen(cmd.split(), cwd=input_folder)

        return result_dict

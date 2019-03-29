import argparse
import os
import os.path
import re
import subprocess
from multiprocessing.pool import Pool

def get_file_extension(filename):
    """
    
    """
    return re.compile(r'^.*?[.](?P<ext>fa\.gz|fasta\.gz|vcf\.gz|\w+)$').match(filename).group('ext')

def get_files(input_path, admited_formats):
    """

    """
    result = []
    folder = "."
    if os.path.isdir(input_path):
        folder = input_path
        result = [os.path.join(folder, f) for f in os.listdir(folder) if os.path.isfile(os.path.join(folder, f)) and get_file_extension(os.path.join(folder,f)) in admited_formats ]
    else:
        filename = os.path.join(folder, input_path)
        if os.path.isfile(filename) and get_file_extension(filename) in admited_formats:
            result = [filename]
    return result

def correct_pdb(input_file):
    """

    """
    print("File %s started correcting" %(input_file))
    fh_r = open(file=input_file, mode="r")
    fh_w = open(file=input_file+".corrected.pdb", mode="w")
    for line in fh_r:
        if line.startswith("ATOM"):
            first_part = line[0:13]
            second_part = line[13:17].strip()
            third_part = line[17:]

            aa_len = 4 - len(second_part)
            for i in range(aa_len):
                second_part = second_part+" "
            fh_w.write(first_part+second_part+third_part)
        else:
            fh_w.write(line)

    fh_r.close()
    fh_w.close()
    print("File %s corrected correctly" %(input_file))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="This program filter vcf files by QUAL and DP")
    parser.add_argument('-i', '--input',
                        dest="input",
                        action="store",
                        default=".",
                        help="The input file/folder")

    options = parser.parse_args()

    formats = set(["pdb"])
    input_files = get_files(input_path=options.input, admited_formats=formats)

    for filename in input_files:
        correct_pdb(filename)

    print("Ended Succesfully!")
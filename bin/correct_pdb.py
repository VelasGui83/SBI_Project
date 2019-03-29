import argparse
import os
import atomator.utils

def correct_pdb(input_file):
    """Rerwite pdb file with left justified column for atom type

    Arguments:
        - input_file - string, the pdb file to correct
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
    os.system("rm "+ input_file)
    print("File %s corrected correctly" %(input_file))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="This program left justify atom type column in a PDB")
    parser.add_argument('-i', '--input',
                        dest="input",
                        action="store",
                        default=".",
                        help="The input PDB file")

    options = parser.parse_args()

    input_files = utils.get_files(input_path=options.input, admited_formats=set(["pdb"]))
    for filename in input_files:
        correct_pdb(filename)

    print("Ended Succesfully!")
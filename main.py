import utils
import STAMPAdapter as stamp
import argparse
import ExcludeWaterSelect
import Alignment as al
from Bio.PDB import PDBIO, PDBParser

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="This program filter vcf files by QUAL and DP")
    parser.add_argument('-i', '--input',
                        dest="input",
                        action="store",
                        default="files/Example/",
                        help="The input file/folder")

    options = parser.parse_args()

    input_files = utils.get_files(input_path=options.input, admited_formats=set(["pdb"]))
    utils.remove_items_containing(input_files, "chain")

    sp = stamp.STAMPParser()
    sp.create_stamp_input(input_files)
    sp.calculate_stamp_scores()
    sp.filter_stamp_scores(9.8)

    remaining_complexes = sp.get_complexes_list()
    first_complex = sp.comlex_dict[remaining_complexes.pop(0)]
    first_chains = first_complex.get_chain_list()
    
    to_evaluate = [first_chains[0], first_chains[1]]

    parser = PDBParser()
    structure = parser.get_structure("final_structure", first_complex.filename)
    model = structure[0]

    while to_evaluate and remaining_complexes:
        evaluating_chain = to_evaluate.pop(0)

        for candidate in evaluating_chain.similar_chains:
            candidate = candidate[0]
            candidate_complex = candidate.parent
            if candidate_complex.id in remaining_complexes:
                added_chain = al.do_superimpose(model, candidate_complex, evaluating_chain, candidate, candidate_complex.complementary_chain(candidate))

                if not al.get_clashes(model, added_chain):
                    remaining_complexes.remove(candidate_complex.id)
                    to_evaluate.append(candidate_complex.complementary_chain(candidate))
                else:
                    model.detach_child(added_chain.id)

    io = PDBIO()
    io.set_structure(structure)
    io.save(file="final.pdb", select=ExcludeWaterSelect.ExcludeWaterSelect())
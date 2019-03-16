import utils
import STAMPAdapter as stamp
import argparse
import ExcludeWaterSelect
import Alignment as al
from Bio.PDB import PDBIO, PDBParser, PDBExceptions

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
    
    to_evaluate = [first_chains[1], first_chains[0]]

    parser = PDBParser()
    structure = parser.get_structure("final_structure", first_complex.filename)
    model = structure[0]

    while to_evaluate and remaining_complexes:
        evaluating_chain = to_evaluate.pop(0)

        for candidate_chain in evaluating_chain.similar_chains:
            candidate_chain = candidate_chain[0]
            candidate_complex = candidate_chain.parent

            complementary_candidate_chain = candidate_complex.complementary_chain(candidate_chain)

            if candidate_complex.id in remaining_complexes:
                added_chain = al.do_superimpose(model, candidate_complex, evaluating_chain.label, candidate_chain.label, complementary_candidate_chain.label)

                if not al.get_clashes_v2(model, added_chain):

                    del candidate_complex.chain_dict[complementary_candidate_chain.label]
                    complementary_candidate_chain.label = added_chain.id
                    candidate_complex.chain_dict[added_chain.id] = complementary_candidate_chain

                    remaining_complexes.remove(candidate_complex.id)
                    to_evaluate.append(candidate_complex.complementary_chain(candidate_chain))
                else:
                    model.detach_child(added_chain.id)

    if remaining_complexes:
        previous_remaining = len(remaining_complexes)
        while True:
            for remaining_complexes_item in remaining_complexes:
                remaining_complexes_item = sp.comlex_dict[remaining_complexes_item]
                for remaining_complexes_chain in remaining_complexes_item.chain_dict.values():
                    complementary_remaining_complexes_chain = remaining_complexes_item.complementary_chain(remaining_complexes_chain)
                    for chain in model.get_chains():
                        if remaining_complexes_item.id in remaining_complexes:
                            try:
                                added_chain = al.do_superimpose(model, remaining_complexes_item, chain.id, remaining_complexes_chain.label, complementary_remaining_complexes_chain.label)

                                if not al.get_clashes_v2(model, added_chain):
                                    del remaining_complexes_item.chain_dict[complementary_remaining_complexes_chain.label]
                                    complementary_remaining_complexes_chain.label = added_chain.id
                                    remaining_complexes_item.chain_dict[added_chain.id] = complementary_remaining_complexes_chain

                                    remaining_complexes.remove(remaining_complexes_item.id)
                                else:
                                    model.detach_child(added_chain.id)
                            except PDBExceptions.PDBException:
                                print("Wrong chain...")

            if previous_remaining == len(remaining_complexes):
                break

            previous_remaining = len(remaining_complexes)


    io = PDBIO()
    io.set_structure(structure)
    io.save(file="final.pdb", select=ExcludeWaterSelect.ExcludeWaterSelect())
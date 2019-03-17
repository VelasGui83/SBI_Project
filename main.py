import STAMPAdapter as stamp
import PDBFiles
import argparse
import ExcludeWaterSelect
import Alignment as al
from Bio.PDB import PDBIO, PDBParser, PDBExceptions
from Bio.PDB.Polypeptide import PPBuilder
import Levenshtein

def get_homologous_tuples(evaluating_chain, remaining_complexes, identity=1.0):
    """

    """
    tuple_list = [(evaluating_chain.parent.filename, evaluating_chain.parent.id, evaluating_chain.original_label)]

    evaluating_structure = PDBParser().get_structure(evaluating_chain.original_label, evaluating_chain.parent.filename)
    evaluating_sequence = ppb.build_peptides(evaluating_structure[0][evaluating_chain.original_label])[0].get_sequence()

    for remaining_complex in remaining_complexes:
        remaining_structure = PDBParser().get_structure(remaining_complex, cd.complexes[remaining_complex].filename)
        for label in cd.complexes[remaining_complex].chain_dict.keys():
            remaining_sequence = ppb.build_peptides(remaining_structure[0][label])[0].get_sequence()
            
            if Levenshtein.ratio(str(evaluating_sequence), str(remaining_sequence)) >= identity:
                tuple_list.append((cd.complexes[remaining_complex].filename, remaining_complex, label))
    return tuple_list

def add_homologous_chains(chain ,homologous_list, complex_dictionary):
    """

    """
    for homologous_chain in homologous_list:
        complex_id = homologous_chain.split("_")[0]
        chain_label = homologous_chain.split("_")[1]
        chain.homologous_chains.add(complex_dictionary.complexes[complex_id].chain_dict[chain_label])

if __name__ == "__main__":
    """

    """
    parser = argparse.ArgumentParser(description="This program filter vcf files by QUAL and DP")
    parser.add_argument('-i', '--input',
                        dest="input",
                        action="store",
                        default="files/Example/",
                        help="The input file/folder")

    options = parser.parse_args()

    cd = PDBFiles.ComplexDictionary()
    cd.populate_complex(options.input)

    remaining_complexes = cd.get_complexes_list()
    first_complex = cd.complexes[remaining_complexes.pop(0)]
    first_chains = first_complex.get_chain_list()

    to_evaluate = [first_chains[1], first_chains[0]]

    parser = PDBParser()
    ppb = PPBuilder()
    sp = stamp.Parser()
    
    structure = parser.get_structure("final_structure", first_complex.filename)
    model = structure[0]

    while to_evaluate and remaining_complexes:
        evaluating_chain = to_evaluate.pop(0)
        
        # sp.create_stamp_input(get_homologous_tuples(evaluating_chain, remaining_complexes, 0.9))
        # stamp_id = evaluating_chain.parent.id+"_"+evaluating_chain.original_label
        # homologous_list = [x[0] for x in sp.get_stamp_scores(chain_id=stamp_id)[stamp_id]]

        homologous_list = [x[1]+"_"+x[2] for x in get_homologous_tuples(evaluating_chain, remaining_complexes, 1)]
        homologous_list.pop(0)


        add_homologous_chains(evaluating_chain, homologous_list, cd)

        for candidate_chain in evaluating_chain.homologous_chains:
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
    io = PDBIO()
    io.set_structure(structure)
    io.save(file="final.pdb", select=ExcludeWaterSelect.ExcludeWaterSelect())
    

    # if remaining_complexes:
    #     previous_remaining = len(remaining_complexes)
    #     while True:
    #         for remaining_complexes_item in remaining_complexes:
    #             remaining_complexes_item = sp.comlex_dict[remaining_complexes_item]
    #             for remaining_complexes_chain in remaining_complexes_item.chain_dict.values():
    #                 complementary_remaining_complexes_chain = remaining_complexes_item.complementary_chain(remaining_complexes_chain)
    #                 for chain in model.get_chains():
    #                     if remaining_complexes_item.id in remaining_complexes:
    #                         try:
    #                             added_chain = al.do_superimpose(model, remaining_complexes_item, chain.id, remaining_complexes_chain.label, complementary_remaining_complexes_chain.label)

    #                             if not al.get_clashes_v2(model, added_chain):
    #                                 del remaining_complexes_item.chain_dict[complementary_remaining_complexes_chain.label]
    #                                 complementary_remaining_complexes_chain.label = added_chain.id
    #                                 remaining_complexes_item.chain_dict[added_chain.id] = complementary_remaining_complexes_chain

    #                                 remaining_complexes.remove(remaining_complexes_item.id)
    #                             else:
    #                                 model.detach_child(added_chain.id)
    #                         except PDBExceptions.PDBException:
    #                             print("Wrong chain...")

    #         if previous_remaining == len(remaining_complexes):
    #             break

    #         previous_remaining = len(remaining_complexes)
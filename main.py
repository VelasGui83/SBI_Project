import STAMPAdapter as stamp
import PDBFiles
import argparse
import ExcludeWaterSelect
import Alignment as al
from Bio.PDB import PDBIO, PDBParser, PDBExceptions
from Bio.PDB.Polypeptide import PPBuilder
import Levenshtein
import copy

def get_homologous_tuples(evaluating_chain, complexes_list, identity=1.0):
    """

    """
    tuple_list = [(evaluating_chain.parent.filename, evaluating_chain.parent.id, evaluating_chain.original_label)]

    evaluating_structure = PDBParser().get_structure(evaluating_chain.original_label, evaluating_chain.parent.filename)
    evaluating_sequence = ppb.build_peptides(evaluating_structure[0][evaluating_chain.original_label])[0].get_sequence()

    for complex_item in complexes_list:
        complex_structure = PDBParser().get_structure(complex_item, cd.complexes[complex_item].filename)
        for label in cd.complexes[complex_item].chain_dict.keys():
            if complex_item != evaluating_chain.parent.id or label != evaluating_chain.label:
                chain_sequence = ppb.build_peptides(complex_structure[0][label])[0].get_sequence()
                if Levenshtein.ratio(str(evaluating_sequence), str(chain_sequence)) >= identity:
                    tuple_list.append((cd.complexes[complex_item].filename, complex_item, label))
    return tuple_list

def add_homologous_chains(chain ,complex_dictionary, compare_complex_list):
    """

    """
    #Without STAMP
    # homologous_list = [x[1]+"_"+x[2] for x in get_homologous_tuples(chain, remaining_complexes, 1)]
    # homologous_list.pop(0)

    #With STAMP
    homologous_tuples = get_homologous_tuples(chain, compare_complex_list, 0.9)
    if len(homologous_tuples) != 1:
        sp.create_stamp_input(homologous_tuples)
        stamp_id = chain.parent.id+"_"+chain.original_label
        scores_dict = sp.get_stamp_scores(chain_id=stamp_id)
        if len(scores_dict) != 0:
            homologous_list = [x[0] for x in scores_dict[stamp_id]]
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

    to_evaluate = [first_chains[0], first_chains[1]]

    parser = PDBParser()
    ppb = PPBuilder()
    sp = stamp.Parser()
    
    structure = parser.get_structure("final_structure", first_complex.filename)
    model = structure[0]

    while to_evaluate and remaining_complexes:
        evaluating_chain = to_evaluate.pop(0)
        add_homologous_chains(evaluating_chain, cd, remaining_complexes)

        for candidate_chain in evaluating_chain.homologous_chains:
            candidate_complex = candidate_chain.parent

            complementary_candidate_chain = candidate_complex.complementary_chain(candidate_chain)

            if candidate_complex.id in remaining_complexes:
                added_chain = al.do_superimpose(model, candidate_complex, evaluating_chain.label, candidate_chain.original_label, complementary_candidate_chain.original_label)

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
    
    # while to_evaluate and remaining_complexes:
    #     evaluating_chain = to_evaluate.pop(0)
    #     add_homologous_chains(evaluating_chain, cd, remaining_complexes)

    #     for candidate_chain in evaluating_chain.homologous_chains:
    #         candidate_complex = candidate_chain.parent

    #         complementary_candidate_chain = candidate_complex.complementary_chain(candidate_chain)

    #         if candidate_complex.id in remaining_complexes:
    #             added_chain = al.do_superimpose(model, candidate_complex, evaluating_chain.label, candidate_chain.original_label, complementary_candidate_chain.original_label)

    #             if not al.get_clashes_v2(model, added_chain):
    #                 del candidate_complex.chain_dict[complementary_candidate_chain.label]
    #                 complementary_candidate_chain.label = added_chain.id
    #                 candidate_complex.chain_dict[added_chain.id] = complementary_candidate_chain

    #                 remaining_complexes.remove(candidate_complex.id)
    #                 to_evaluate.append(candidate_complex.complementary_chain(candidate_chain))
    #             else:
    #                 model.detach_child(added_chain.id)
    


    # for complex_item in cd.complexes.values():
    #     for chain_item in complex_item.chain_dict.values():
    #         add_homologous_chains(chain_item, cd, cd.get_complexes_list())

    # while to_evaluate:
    #     evaluating_chain = to_evaluate.pop(0)

    #     for candidate_chain in evaluating_chain.homologous_chains:
    #         candidate_complex = candidate_chain.parent

    #         complementary_candidate_chain = candidate_complex.complementary_chain(candidate_chain)

    #         added_chain = al.do_superimpose(model, candidate_complex, evaluating_chain.label, candidate_chain.original_label, complementary_candidate_chain.original_label)

    #         if not al.get_clashes_v2(model, added_chain):
    #             next_to_evaluate = copy.copy(complementary_candidate_chain)
    #             next_to_evaluate.label = added_chain.id
                
    #             to_evaluate.append(next_to_evaluate)
    #         else:
    #             model.detach_child(added_chain.id)
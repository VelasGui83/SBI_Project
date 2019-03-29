import STAMPAdapter as stamp
import PDBFiles
import argparse
import ExcludeWaterSelect
import Alignment as al
from Bio.PDB import PDBIO, PDBParser, PDBExceptions, MMCIFIO
from Bio.PDB.Polypeptide import PPBuilder
import Levenshtein
import copy

def get_nucleotides(chain):
    """

    """
    return "".join([x.resname.strip() for x in chain.child_list])

def get_homologous_tuples(evaluating_chain, complexes_list, identity=1.0):
    """

    """
    ppb = PPBuilder()

    tuple_list = [(evaluating_chain.parent.filename, evaluating_chain.parent.id, evaluating_chain.original_label)]
    chain_type = "prot"

    evaluating_structure = PDBParser().get_structure(evaluating_chain.original_label, evaluating_chain.parent.filename)
    evaluating_polypeptide = ppb.build_peptides(evaluating_structure[0][evaluating_chain.original_label])
    evaluating_sequence = None
    if evaluating_polypeptide:       
        evaluating_sequence = evaluating_polypeptide[0].get_sequence()
    else:
        chain_type = "dna/rna"
        evaluating_sequence = get_nucleotides(evaluating_structure[0][evaluating_chain.original_label])

    for complex_item in complexes_list:
        complex_structure = PDBParser().get_structure(complex_item, cd.complexes[complex_item].filename)
        for label in cd.complexes[complex_item].chain_dict.keys():
            if complex_item != evaluating_chain.parent.id or label != evaluating_chain.label:
                chain_polypeptide = ppb.build_peptides(complex_structure[0][label])
                chain_sequence = None
                if chain_polypeptide:
                    chain_sequence = chain_polypeptide[0].get_sequence()
                else:
                    chain_sequence = get_nucleotides(complex_structure[0][label])
                if Levenshtein.ratio(str(evaluating_sequence), str(chain_sequence)) >= identity:
                    tuple_list.append((cd.complexes[complex_item].filename, complex_item, label))
    return (chain_type, tuple_list)

def get_homologous_chains(chain ,complex_dictionary, compare_complex_list):
    """

    """
    homologous_set = set([])
    
    homologous_tuples = get_homologous_tuples(chain, compare_complex_list, 0.9)
    if len(homologous_tuples[1]) != 1 and homologous_tuples[0] == "prot":
        #With STAMP
        sp.create_stamp_input(homologous_tuples[1])
        stamp_id = chain.parent.id+"_"+chain.original_label
        scores_dict = sp.get_stamp_scores(chain_id=stamp_id)
        if len(scores_dict) != 0:
            homologous_list = [x[0] for x in scores_dict[stamp_id]]
            for homologous_chain in homologous_list:
                complex_id = homologous_chain.split("_")[0]
                chain_label = homologous_chain.split("_")[1]
                homologous_set.add(complex_dictionary.complexes[complex_id].chain_dict[chain_label])
    else:
        # #Without STAMP
        homologous_list = get_homologous_tuples(chain, compare_complex_list, 1)
        homologous_list[1].pop(0)
        for homologous_chain in homologous_list[1]:
            complex_id = homologous_chain[1]
            chain_label = homologous_chain[2]
            homologous_set.add(complex_dictionary.complexes[complex_id].chain_dict[chain_label])
    return homologous_set

def create_macrocomplex(model, to_evaluate):
    while to_evaluate:
        evaluating_chain = to_evaluate.pop(0)

        for candidate_chain in evaluating_chain.homologous_chains:
            candidate_complex = candidate_chain.parent

            complementary_candidate_chain = candidate_complex.complementary_chain(candidate_chain)

            added_chain = al.do_superimpose(model, candidate_complex, evaluating_chain.label, candidate_chain.original_label, complementary_candidate_chain.original_label)

            if not al.get_clashes_v3(model, added_chain):
                next_to_evaluate = copy.copy(complementary_candidate_chain)
                next_to_evaluate.label = added_chain.id
                to_evaluate.append(next_to_evaluate)

                print("Chain added. #%d" %len(model))

                if chain_number == 200:
                    break
            else:
                model.detach_child(added_chain.id)
    return len(model)

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
    parser = PDBParser()
    sp = stamp.Parser()
    cd = PDBFiles.ComplexDictionary()

    cd.populate_complex(options.input)
    remaining_complexes = cd.get_complexes_list()

    for complex_item in cd.complexes.values():
        for chain_item in complex_item.chain_dict.values():
            if not chain_item.homologous_chains:
                homologous_set = get_homologous_chains(chain_item, cd, cd.get_complexes_list())
                homologous_set.add(chain_item)
                for chain in homologous_set:
                    homologous_specific_set = copy.copy(homologous_set)
                    homologous_specific_set.remove(chain)
                    chain.homologous_chains = homologous_specific_set

    model = None
    chain_number = 2
    while chain_number == 2:
        first_complex = cd.complexes[remaining_complexes.pop(0)]
        first_chains = first_complex.get_chain_list()
        to_evaluate = [first_chains[0], first_chains[1]]

        structure = parser.get_structure("final_structure", first_complex.filename)
        model = structure[0]

        chain_number = create_macrocomplex(model, to_evaluate)

    if chain_number > 52:
        io = MMCIFIO()
        io.set_structure(structure)
        io.save(filepath="final.cif", select=ExcludeWaterSelect.ExcludeWaterSelect())
    else:
        io = PDBIO()
        io.set_structure(structure)
        io.save(file="final.pdb", select=ExcludeWaterSelect.ExcludeWaterSelect())
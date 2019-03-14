import utils
import STAMPParser as stamp
import argparse
import ExcludeWaterSelect
import clashes
from Bio.PDB import PDBParser, Superimposer, PDBIO

chain_dict = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z', 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z', 'AA', 'BB', 'CC', 'DD', 'EE', 'FF', 'GG', 'HH', 'II', 'JJ', 'KK', 'LL', 'MM', 'NN', 'OO', 'PP', 'QQ', 'RR', 'SS', 'TT', 'UU', 'VV', 'WW', 'XX', 'YY', 'ZZ', 'aa', 'bb', 'cc', 'dd', 'ee', 'ff', 'gg', 'hh', 'ii', 'jj', 'kk', 'll', 'mm', 'nn', 'oo', 'pp', 'qq', 'rr', 'ss', 'tt', 'uu', 'vv', 'ww', 'xx', 'yy', 'zz']

def next_available_chain_id(structure):
    letters_used = [x._id for x in structure.get_chains()]
    for letter in chain_dict:
        if letter not in letters_used:
            return letter
    return None

def do_superimpose(fixed_structure, moveable_complex, fixed_chain, target_chain, moveable_chain):
    parser = PDBParser()
    sup = Superimposer()

    s2 = parser.get_structure(moveable_complex, moveable_complex[target_chain])

    s1_fixed_chain = fixed_structure[0][fixed_chain]
    s2_target_chain = s2[0][target_chain]

    s2_moved_chain = s2[0][moveable_chain]
    s2_moved_chain_atoms = []

    s2_moved_chain._id = next_available_chain_id(fixed_structure)

    # Extract information of atoms
    for residues in s2_moved_chain:
        s2_moved_chain_atoms.append(residues)

    # Specify the atom lists
    # 'fixed' and 'moving' are lists of Atom objects
    # The moving atoms will be put on the fixed atoms
    sup.set_atoms(list(s1_fixed_chain.get_atoms()), list(s2_target_chain.get_atoms()))

    # Apply rotation/translation to the moving atoms
    sup.apply(s2_moved_chain_atoms)

    fixed_structure[0].add(s2_moved_chain)

    return s2_moved_chain

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="This program filter vcf files by QUAL and DP")
    parser.add_argument('-i', '--input',
                        dest="input",
                        action="store",
                        default="files/",
                        help="The input file/folder")

    options = parser.parse_args()

    input_files = utils.get_files(input_path=options.input, admited_formats=set(["pdb"]))
    utils.remove_items_containing(input_files, "chain")

    sp = stamp.STAMPParser()
    sp.create_stamp_input(input_files)
    sp.calculate_stamp_scores()
    sp.filter_stamp_scores(9.8)

    complexes_pairs = sp.get_complexes_list()

    first_chain = max(sp.stamp_scores, key=lambda x: len(sp.stamp_scores[x]))
    second_chain = sp.complementary_chain(first_chain)
    
    to_evaluate = [first_chain, second_chain]
    complexes_pairs.remove(first_chain[:-2])

    parser = PDBParser()
    structure = parser.get_structure(first_chain[:-2], sp.complex_files[first_chain[:-2]][first_chain[-1]])

    while to_evaluate and complexes_pairs:
        evaluating_chain = to_evaluate.pop(0)
        superimposer_candidates = sp.stamp_scores[evaluating_chain]

        for candidate in superimposer_candidates:
            candidate = candidate[0]
            candidate_complex = candidate[:-2]
            if candidate_complex in complexes_pairs:
                added_chain = do_superimpose(structure, sp.complex_files[candidate_complex], evaluating_chain[-1], candidate[-1], sp.complementary_chain(candidate)[-1])

                if not clashes.get_clashes(structure, added_chain):
                    complexes_pairs.remove(candidate_complex)
                    to_evaluate.append(candidate_complex+"_"+added_chain._id)
                else:
                    structure.detach_child(added_chain._id)

    io = PDBIO()
    io.set_structure(structure)
    io.save(file="final.pdb", select=ExcludeWaterSelect.ExcludeWaterSelect())
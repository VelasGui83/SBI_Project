import utils

class Chain(object):
    """

    """
    def __init__(self, label, parent):
        """

        """
        self.label = label
        self.original_label = label
        self.parent = parent
        self.homologous_chains = set([])

    def __len__(self):
        """

        """
        return len(self.homologous_chains)

class Complex(object):
    """

    """
    def __init__(self, id, filename):
        """

        """
        self.id = id
        self.chain_dict = {}
        self.filename = filename

    def get_chain_list(self):
        """

        """
        return list(self.chain_dict.values())

    def complementary_chain(self, chain):
        """

        """
        for chain_item in self.chain_dict.values():
            if chain_item != chain:
                return chain_item
        return None
    
    def filter_chain_scores(self, limit):
        """

        """
        for chain in self.chain_dict.values():
            chain.filter_similar_chains(limit)
            chain.order_similar_chains()

class ComplexDictionary(object):
    """

    """
    def __init__(self):
        self.complexes = {}
        
    def get_complexes_list(self):
        """

        """
        return list(self.complexes.keys())

    def _get_chains(self, input_file):
        set_chains = set([])
        with open(input_file, 'r') as fh:
            for line in fh:
                if line.startswith("ATOM"):
                    set_chains.add(line[21])
        return set_chains

    def populate_complex(self, input_folder):
        """

        """
        input_files = utils.get_files(input_path=input_folder, admited_formats=set(["pdb"]))
        utils.remove_items_containing(input_files, "chain")

        i = 1
        for filename in input_files:
            chains = sorted(list(self._get_chains(filename)))
            for chain in chains:
                complex_id = "%d:%d" %(i, i+1)
                try:
                    self.complexes[complex_id].chain_dict[chain] = Chain(chain, self.complexes[complex_id])
                except KeyError:
                    self.complexes[complex_id] = Complex(complex_id, filename)
                    self.complexes[complex_id].chain_dict[chain] = Chain(chain, self.complexes[complex_id])
            i += 2
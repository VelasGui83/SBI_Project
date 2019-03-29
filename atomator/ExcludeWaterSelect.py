from Bio.PDB import Select

class ExcludeWaterSelect(Select):
    """

    """
    
    def accept_residue(self, residue):
        """

        """
        if residue.get_full_id()[3][0] == " ":
            return 1
        return 0
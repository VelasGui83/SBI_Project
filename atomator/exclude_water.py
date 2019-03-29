from Bio.PDB import Select

class ExcludeWaterSelect(Select):
    """Subclass of Bio.PDB.Select to excude waters"""
    
    def accept_residue(self, residue):
        """Excude waters"""
        if residue.get_full_id()[3][0] == " ":
            return 1
        return 0
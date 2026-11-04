from Bio.PDB import PDBParser, PDBIO, Select

KEEP = {94, 96, 119, 199}

class MotifOnly(Select):
    def accept_residue(self, residue):
        hetflag, resseq, icode = residue.id
        if residue.get_parent().id != "A":
            return 0
        if hetflag == " " and resseq in KEEP:
            return 1
        if residue.get_resname().strip() == "ZN":
            return 1
        return 0

parser = PDBParser(QUIET=True)
structure = parser.get_structure("ca2", "inputs/3KS3.pdb")
io = PDBIO()
io.set_structure(structure)
io.save("inputs/motif.pdb", MotifOnly())
print("wrote inputs/motif.pdb")

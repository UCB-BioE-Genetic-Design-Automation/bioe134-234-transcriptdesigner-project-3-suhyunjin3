import random
from genedesign.rbs_chooser import RBSChooser
from genedesign.models.transcript import Transcript
from genedesign.checkers.forbidden_sequence_checker import ForbiddenSequenceChecker

class TranscriptDesigner:
    def __init__(self):
        self.codon_table = {} 
        self.rbsChooser = None
        self.forbiddenChecker = ForbiddenSequenceChecker()
        self.codonChecker = CodonChecker() # ADD THIS

    def run(self, peptide: str, ignores: set) -> Transcript:
        forbidden_sites = ["GAATTC", "GGTACC", "GCTAGC"]
        
        while True:
            # 1. Generate the sequence (use your stochastic logic)
            codons = self.generate_stochastic_codons(peptide)
            cds = ''.join(codons)
            
            # 2. Check Forbidden Sites
            if self.forbiddenChecker.run(cds, forbidden_sites):
                continue
                
            # 3. Check Codon Quality (The part that was failing!)
            # This is how you ensure your design passes the tests
            above_board, diversity, rare_count, cai = self.codonChecker.run(codons)
            if not above_board:
                continue # Try again if CAI is too low or too many rare codons
                
            break # Success!

        selectedRBS = self.rbsChooser.run(cds, ignores)
        return Transcript(selectedRBS, peptide, codons)
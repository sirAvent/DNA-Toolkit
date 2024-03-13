# DNAToolkit test file
import unittest
from DNAToolkit import *
from NucleicAcid import NucleicAcid
import random

# Test class
class TestDNAToolKitFunctions(unittest.TestCase):
    def test_validate_seq_true(self):
        """Test cases for validateSeq function"""
        # Test 1: Validate Sequence is True
        sequence_to_test = "ACGTCGATCGGTA"
        bclass = NucleicAcid(sequence_to_test)
        result = validateSeq(sequence_to_test)
        self.assertTrue(result, f"Test 1 Failed. Expected: True, Actual: {result}")

    def test_validate_seq_false(self):
        # Test 2: Validate Sequence is False
        sequence_to_test = "ACGTCGATCGGTAZ"
        result = validateSeq(sequence_to_test)
        self.assertFalse(result, f"Test 2 Failed. Expected: False, Actual: {result}")

        sequence_to_test = "AcGTCgATCGGTa"
        result = validateSeq(sequence_to_test)
        self.assertFalse(result, f"Test 2 Failed. Expected: False, Actual: {result}")

    def test_validate_seq_empty(self):
        # Test 3: Validate Sequence is False when sequence is empty
        sequence_to_test = ""
        result = validateSeq(sequence_to_test)
        self.assertFalse(result, f"Test 3 Failed. Expected: False, Actual: {result}")

    def test_nucleotide_freq(self):
      ''' Test cases for countNucFrequency function'''
      # Test 4: Returns correct result
      sequence_to_test = "ACGTCGATCGGTA"
      result = countNucFrequency(sequence_to_test)
      self.assertEqual(result, {'A':3, 'C':3, 'G':4, 'T':3}, f"Test 4 Failed. Nucleotide Frequency is incorrect.")

    def test_transcription(self):
      ''' Test cases for transcription function'''
      # Test 5: Returns correct result
      sequence_to_test = "ACGTCGATCGGTA"
      result = transcription(sequence_to_test)
      self.assertEqual(len(result), len(sequence_to_test), f"Test 5 Failed. Lengths do not match")

    def test_reverse_complement(self):
      ''' Test case for reverse complment'''
      # Test 6: Returns correct result
      sequence_to_test = "TTGTCAGTATTTGACCATCCCTTTACCATGGATGGCGAAAAATAAAGCAG"
      result = reverse_complement(sequence_to_test)
      self.assertEqual(result, "CTGCTTTATTTTTCGCCATCCATGGTAAAGGGATGGTCAAATACTGACAA", f"Test 6 Failed.")
    
    def test_gc_content(self):
      ''' Test case for GC Content'''
      # Test 7: Returns correct result
      sequence_to_test = "TTGTCAGTATTTGACCATCCCTTTACCATGGATGGCGAAAAATAAAGCAG"
      result = gc_content(sequence_to_test)
      self.assertEqual(result, 40, f"Test 7 Failed.")

    def test_gc_content_subsections(self):
      ''' Test case for GC Content subsections'''
      # Test 8: Returns correct result
      sequence_to_test = "TTGTCAGTATTTGACCATCCCTTTACCATGGATGGCGAAAAATAAAGCAG"
      bclass = NucleicAcid(sequence_to_test)
      result = gc_content_subsections(sequence_to_test, 5)
      result = bclass.get_gc_content_subsections(5)
      self.assertEqual(result, [40, 20, 40, 60, 20, 60, 60, 40, 0, 60], f"Test 8 Failed.")

    def test_translate_seq(self):
      ''' Test case for translate_seq function'''
      # Test 9: Returns correct result
      sequence_to_test = "TTGTCAGTATTTGACCATCCCTTTACCATGGATGGCGAAAAATAAAGCAG"
      result = translate_seq(sequence_to_test)
      print(codon_freq(sequence_to_test, "L"))
      self.assertEqual(result, ['L', 'S', 'V', 'F', 'D', 'H', 'P', 'F', 'T', 'M', 'D', 'G', 'E', 'K', '_', 'S'], f"Test 9 Failed.")


    def test_codon_freq(self):
      ''' Test case for codon_freq function'''
      # Test 10: Returns correct result
      sequence_to_test = "TTGTCAGTATTTGACCATCCCTTTACCATGGATGGCGAAAAATAAAGCAG"
      bclass = NucleicAcid(sequence_to_test)
      result = codon_freq(sequence_to_test, "L")
      result = bclass.get_codon_frequency("L")
      self.assertEqual(result, {'TTG': 1.0}, f"Test 10 Failed.")

if __name__ == '__main__':
    # TODO: Test NucleicAcid class
    # unittest.main()
    seq = ''.join([random.choice(NUCLEOTIDE_BASE["DNA"]) for x in range(40)])
    test_dna = NucleicAcid(seq, "DNA", "Randomly generated sequence")
    

    print(test_dna.get_info())
    print(test_dna.count_nucleotides())
    print(test_dna.get_transcription())
    print(test_dna.get_reverse_complement())
    print(test_dna.get_gc_content())
    print(test_dna.get_gc_content_subsections())
    print(test_dna.get_translation())
    print(test_dna.get_codon_frequency('L'))

    for rf in test_dna.get_reading_frames():
      print(rf)


    print(test_dna.get_proteins_from_otf())

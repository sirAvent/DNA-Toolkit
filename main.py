# DNAToolkit test file
import unittest
from DNAToolkit import *

# Test class
class TestDNAToolKitFunctions(unittest.TestCase):
    def test_validate_seq_true(self):
        """Test cases for validateSeq function"""
        # Test 1: Validate Sequence is True
        sequence_to_test = "ACGTCGATCGGTA"
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
      result = gc_content_subsections(sequence_to_test, 5)
      self.assertEqual(result, [40, 20, 40, 60, 20, 60, 60, 40, 0, 60], f"Test 8 Failed.")


if __name__ == '__main__':
    unittest.main()
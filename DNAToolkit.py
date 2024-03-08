# DNA Toolkit
from structures import *

def validateSeq(dna_seq: str) -> bool:
  '''
  Check if a DNA sequence string has valid nucleotides(A, C, G, T) and is capitalized.
  Return true if the sequence is valid
  '''
  for n in dna_seq:
    if n not in Nucleotides:
      return False
  return False if len(dna_seq) == 0 else True


def countNucFrequency(dna_seq: str) -> dict:
  '''
  Count the frequency of the 4 nucleotides in a DNA sequence string.
  Return a dictionary with the frequency of each nucleotide in a valid DNA sequence string
  '''
  nucDict = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
  for n in dna_seq:
    nucDict[n] += 1
  return nucDict


def transcription(dna_seq: str) -> str:
  '''
  Replace Thymine with Uracil.
  Return a RNA transcription of the DNA sequence
  '''
  # DNA is transcribed into RNA to be used to create proteins from codons
  return dna_seq.replace('T', 'U')


def reverse_complement(dna_seq: str) -> str:
  '''
  Return a reverse-complement counterpart of the DNA sequence
  '''
  mapping = str.maketrans('ATCG', 'TAGC')
  return dna_seq.translate(mapping)[::-1]


def gc_content(dna_seq: str) -> int:
  '''
  Return a perecentage of the GC content in a DNA/RNA sequence
  '''
  return round((dna_seq.count('C') + dna_seq.count('G')) / len(dna_seq) * 100)


def gc_content_subsections(dna_seq: str, k=20) -> list:
  '''
  Divide a DNA/RNA sequence into windows determined by a specified parameter, k.
  Return a list of the GC content in those windows/subsections
  '''
  result = []
  for i in range(0, len(dna_seq) - k + 1, k):
    subsection = dna_seq[i:i + k]
    result.append(gc_content(subsection))
  return result
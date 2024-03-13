# DNA Toolkit
from structures import *
from collections import Counter

def validateSeq(dna_seq: str) -> bool:
  '''
  Check if a DNA sequence string has valid nucleotides(A, C, G, T) and is capitalized.
  Return true if the sequence is valid
  '''
  for n in dna_seq:
    if n not in NUCLEOTIDE_BASE["DNA"]:
      return False
  return False if len(dna_seq) == 0 else True


def countNucFrequency(dna_seq: str) -> dict:
  '''
  Count the frequency of the 4 nucleotides in a DNA sequence string.
  Return a dictionary with the frequency of each nucleotide in a valid DNA sequence string.
  '''
  return dict(Counter(dna_seq))


def transcription(dna_seq: str) -> str:
  '''
  Replace Thymine with Uracil.
  Return a RNA transcription of the DNA sequence.
  '''
  # DNA is transcribed into RNA to be used to create proteins from codons
  return dna_seq.replace('T', 'U')


def reverse_complement(dna_seq: str) -> str:
  '''
  Return a reverse-complement counterpart of the DNA sequence.
  '''
  mapping = str.maketrans('ATCG', 'TAGC')
  return dna_seq.translate(mapping)[::-1]


def gc_content(dna_seq: str) -> int:
  '''
  Return a perecentage of the GC content in a DNA/RNA sequence.
  '''
  return round((dna_seq.count('C') + dna_seq.count('G')) / len(dna_seq) * 100)


def gc_content_subsections(dna_seq: str, k=20) -> list:
  '''
  Divide a DNA/RNA sequence into windows determined by a specified parameter, k.
  Return a list of the GC content in those windows/subsections.
  '''
  result = []
  for i in range(0, len(dna_seq) - k + 1, k):
    subsection = dna_seq[i:i + k]
    result.append(gc_content(subsection))
  return result


def translate_seq(dna_seq: str, init_pos=0) -> list:
  '''
  Translates a DNA sequence into an aminoacid sequence.
  Returns a list of codons from a DNA sequence.
  '''
  return [DNA_Codons[dna_seq[pos:pos + 3]] for pos in range(init_pos, len(dna_seq) - 2, 3)]


def codon_freq(dna_seq: str, aminoacid: str) -> dict:
  '''
  Return the frequency of a codon encoding a given aminoacid in a DNA sequence.
  Frequency is returned as a ratio over total codons found for a given aminoacid.
  '''
  codon_list = []
  for i in range(0, len(dna_seq) - 2, 3):
    if DNA_Codons[dna_seq[i:i + 3]] == aminoacid:
      codon_list.append(dna_seq[i:i + 3])
  
  result = dict(Counter(codon_list))
  totalWight = sum(result.values())
  for seq in result:
    result[seq] = round(result[seq] / totalWight, 2)
  return result


def get_reading_frames(dna_seq: str) -> list:
  '''
  Generate the six reading frames of a given DNA sequence, including the reverse-complement
  Return a list of 6 reading frames which are in string format.
  '''
  frames = []
  frames.append(translate_seq(dna_seq))
  frames.append(translate_seq(dna_seq, 1))
  frames.append(translate_seq(dna_seq, 2))
  frames.append(translate_seq(reverse_complement(dna_seq)))
  frames.append(translate_seq(reverse_complement(dna_seq, 1)))
  frames.append(translate_seq(reverse_complement(dna_seq, 2)))
  return frames


def get_proteins_from_reading_frame(aminoacid_seq: list|str) -> list:
  '''
  Return a list of all possible proteins in an aminoacid sequence.
  '''
  current_protein = []
  proteins = []
  for a in aminoacid_seq:
    if a == "_":
      # STOP making amino acids
      if current_protein:
        for p in current_protein:
          proteins.append(p)
        current_protein = []
    else:
      if a == "M":
        current_protein.append("")
      for i in range(len(current_protein)):
        current_protein[i] += a
  return proteins


def get_all_proteins(dna_seq: str, startPos=0, endPos=0, ordered=False) -> list:
  '''
  Takes in a list of 6 opening reading frames to get proteins generated.
  Return all possible proteins for all open reading frames.
  '''
  if endPos > startPos:
    rfs = get_reading_frames(dna_seq[startPos: endPos])
  else:
    rfs = get_reading_frames(dna_seq)
  
  result = []
  for rf in rfs:
    proteins = get_proteins_from_reading_frame(rf)
    for p in proteins:
      result.append(p)

  return sorted(result, key=len, reverse=True) if ordered else result
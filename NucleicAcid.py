from collections import Counter
from structures import DNA_Codons


class NucleicAcid:
  ''' DNA sequence class. Default values: ATCG, DNA, No label '''

  def __init__(self, seq="ATCG", seq_type="DNA", label="No label") -> None:
    self.seq = seq.upper()
    self.seq_type = seq_type
    self.label = label
    assert self.__validate(), f"Provided DNA sequence is not valid: {self.seq}"
  

  def __validate(self) -> bool:
    '''
    Check if a DNA sequence string has valid nucleotides(A, C, G, T) and is capitalized.
    Return true if the sequence is valid
    '''
    return {'A', 'C', 'G', 'T'}.issuperset(self.seq)
  

  def get_type(self) -> str:
    ''' Return sequence type'''
    return self.seq_type


  def get_info(self) -> str:
    ''' Return a string containing all information about the class '''
    return f"[Label]: {self.label}\n[Sequence]:{self.seq}\n[Type]: {self.seq_type}"
  

  def count_nucleotides(self) -> dict:
    '''
    Return a dictionary with the frequency of each nucleotide in a valid DNA sequence string.
    '''
    return dict(Counter(self.seq))
  

  def get_transcription(self) -> str:
    '''
    Return an RNA transcription of the DNA sequence.
    '''
    # DNA is transcribed into RNA to be used to create proteins from codons
    return self.seq.replace('T', 'U')
  

  def get_reverse_complement(self) -> str:
    '''
    Return a reverse-complement counterpart of the DNA sequence.
    '''
    mapping = str.maketrans('ATCG', 'TAGC')
    return self.seq.translate(mapping)[::-1]
  

  def get_gc_content(self) -> int:
    '''
    Return a perecentage of the GC content in a DNA/RNA sequence.
    '''
    return round((self.seq.count('C') + self.seq.count('G')) / len(self.seq) * 100)
  

  def get_gc_content_subsections(self, k=20) -> list:
    '''
    Divide a DNA/RNA sequence into windows determined by a specified parameter, k.
    Return a list of the GC content in those windows/subsections.
    '''
    result = []
    for i in range(0, len(self.seq) - k + 1, k):
      subsection = self.seq[i:i + k]
      result.append(round((subsection.count('C') + subsection.count('G')) / len(subsection) * 100))
    return result
  

  def get_dna_translation(self, init_pos=0) -> list:
    '''
    Returns a list of codons/amino acids from a DNA sequence.
    '''
    return [DNA_Codons[self.seq[pos:pos + 3]] for pos in range(init_pos, len(self.seq) - 2, 3)]
  

  def get_codon_frequency(self, aminoacid: str) -> dict:
    '''
    Return the frequency of a codon encoding a given aminoacid in a DNA sequence.
    Frequency is returned as a ratio over total codons found for a given aminoacid.
    '''
    codon_list = []
    for i in range(0, len(self.seq) - 2, 3):
      if DNA_Codons[self.seq[i:i + 3]] == aminoacid:
        codon_list.append(self.seq[i:i + 3])
    
    result = dict(Counter(codon_list))
    totalWight = sum(result.values())
    for seq in result:
      result[seq] = round(result[seq] / totalWight, 2)
    return result
  
  
  def get_reading_frames(self) -> list:
    '''
    Generate the six reading frames of a given DNA sequence, including the reverse-complement
    Return a list of 6 reading frames which are in string format.
    '''
    frames = []
    frames.append(self.get_dna_translation(0))
    frames.append(self.get_dna_translation(1))
    frames.append(self.get_dna_translation(2))
    tmp_seq = NucleicAcid(self.get_reverse_complement(), self.seq_type)
    frames.append(tmp_seq.get_dna_translation(0))
    frames.append(tmp_seq.get_dna_translation(1))
    frames.append(tmp_seq.get_dna_translation(2))
    del tmp_seq
    return frames
  
  def get_proteins_from_rf(self, aminoacid_seq: list|str) -> list:
    '''
    Return a list of all possible proteins in an aminoacid sequence (one reading frame).
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
  

  
  def get_proteins_from_otf(self, startPos=0, endPos=0, ordered=False) -> list:
    '''
    Takes in a list of 6 opening reading frames to get proteins generated.
    Return all possible proteins for all open reading frames.
    '''
    if endPos > startPos:
      tmp_seq = NucleicAcid(self.seq[startPos: endPos])
      rfs = tmp_seq.get_reading_frames()
    else:
      rfs = self.get_reading_frames()
    
    result = []
    for rf in rfs:
      proteins = self.get_proteins_from_rf(rf)
      for p in proteins:
        result.append(p)

    return sorted(result, key=len, reverse=True) if ordered else result
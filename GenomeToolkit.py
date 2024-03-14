class GenomeToolkit:
  def __init__(self):
    print("Genome Toolkit initialized")

  
  def count_kmer(self, seq: str, kmer: str) -> int:
    ''' Count repeating k-mers in a sequences. Includes overlapping k-mers '''
    kmer_count = 0
    for pos in range(len(seq) - len(kmer) - 1):
      if seq[pos:pos + len(kmer)] == kmer:
        kmer_count += 1
      return kmer_count
    
  
  def get_most_frequent_kmers(self, seq: str, k_len: int) -> list:
    ''' Find the most frequent k-mers of a given length in a DNA string '''
    kmer_freq = dict()

    for i in range(len(seq) - k_len + 1):
      kmer = seq[i:i+k_len]
      if kmer in kmer_freq:
        kmer_freq[kmer] += 1
      else:
        kmer_freq[kmer] = 1
    
    highest_freq = max(kmer_freq.values())

    return [kmer for kmer, freq in kmer_freq.items() if freq == highest_freq]
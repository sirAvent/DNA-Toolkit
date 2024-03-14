def readTextFile(filePath: str) -> str:
  ''' Return the text of a text file '''
  with open(filePath, 'r') as f:
    return "".join([line.strip() for line in f.readlines])
  

def writeTextFile(filePath: str, seq: str, mode='w') -> None:
  ''' Write a sequence on a text file. Modes(w-write, a-append) '''
  with open(filePath, mode) as f:
    f.write(seq + '\n')


def readFASTA(filePath: str) -> dict:
  ''' Return a dictionary containing the parsed content of the FASTA File. '''
  with open(filePath, 'r') as f:
    FASTAFile = [line.strip() for line in f.readlines]

  FASTADict = dict()
  FASTALabel = ""

  for line in FASTAFile:
    if '>' in line:
      FASTALabel = line
      FASTADict[FASTALabel] = ""
    else:
      FASTADict[FASTALabel] += line
  
  return FASTADict
import numpy as np


def read_fasta(path):
    container = list()
    letters = {'A', 'C', 'G', 'T'}
    with open(path) as file:
        for line in file:
            line = line.strip().upper()
            if not line.startswith('>'):
                seq = ''.join([l if l in letters else 'N' for l in line])
                complement_seq = complement(seq)
                container.append(seq)
                container.append(complement_seq)
    return(container)
    
    
def complement(seq):
    seq = seq.replace('A', 't')
    seq = seq.replace('T', 'a')
    seq = seq.replace('C', 'g')
    seq = seq.replace('G', 'c')
    seq = seq.upper()
    seq = seq[::-1]
    return(seq)

    
def shuffle_fasta(fasta, times):
    out = []
    for seq in fasta:
        seq = np.asarray(list(seq))
        for i in range(times):
            out.append(np.copy(seq))
            np.random.shuffle(out[-1])
    out = [''.join(i) for i in out]
    return(out)


def get_number_of_sites(peaks, length):
    n = 0
    for p in peaks:
        n += len(p) - length + 1
    return(n)


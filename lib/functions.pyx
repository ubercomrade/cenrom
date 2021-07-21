from random import sample

def read_fasta(str path):
    cdef list l, fasta = []
    cdef str line
    with open(path, 'r') as file:
        for line in file:
            if not line.startswith('>'):
                seq = line.strip().upper()
                fasta.append(line.strip().upper())
    file.close()
    return(fasta)


def complement(seq):
    seq = seq.replace('A', 't')
    seq = seq.replace('T', 'a')
    seq = seq.replace('C', 'g')
    seq = seq.replace('G', 'c')
    seq = seq.upper()
    seq = seq[::-1]
    return(seq)


def shuffle_fasta(list fasta, int times):
    cdef list shuffled, out = []
    for seq in fasta:
        shuffled = [''.join(sample(seq, k=len(seq))) for i in range(times)]
        out += shuffled
    return(out)
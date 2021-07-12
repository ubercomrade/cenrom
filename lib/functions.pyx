from random import sample

def read_fasta(str path):
    cdef list l, fasta = []
    cdef str line
    with open(path, 'r') as file:
        for line in file:
            if not line.startswith('>'):
                fasta.append(line.strip().upper())
    file.close()
    return(fasta)


def shuffle_fasta(list fasta, int times):
    cdef list shuffled, out = []
    for seq in fasta:
        shuffled = [''.join(sample(seq, k=len(seq))) for i in range(times)]
        out += shuffled
    return(out)
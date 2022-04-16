import numpy as np
import MOODS.tools
import MOODS.scan


class PWM:    
    def __init__(self, path, form):
        if form == 'hocomoco':
            self.__matrix = PWM.__read_matrix(path, 0)
            self.__matrix = PWM.__pcm_to_pfm(self.__matrix)
            self.__matrix = PWM.__pfm_to_pwm(self.__matrix)
        if form == 'cisbp':
            self.__matrix = PWM.__read_matrix(path, 1)
            self.__matrix = PWM.__pfm_to_pwm_cisbp(self.__matrix)
        if form == 'homer':
            self.__matrix = PWM.__read_matrix(path, 0)
            self.__matrix = PWM.__pfm_to_pwm(self.__matrix)
        if form == 'pwm':
            self.__matrix = PWM.__read_matrix(path, 0)
        self.__length = self.__matrix.shape[1]
        self.__table = list()
    
    @property
    def matrix(self):
        return self.__matrix
    
    @property
    def length(self):
        return self.__length
    
    @property
    def table(self):
        return self.__table
    
    @staticmethod
    def __read_matrix(path, pos):
        matrix = []
        with open(path) as file:
            file.readline()
            for line in file:
                line = [float(i) for i in line.strip().split('\t')[pos:]]
                matrix.append(line)
        return np.array(matrix).T

    @staticmethod
    def __pfm_to_pwm(pfm):
        background = 0.25
        pwm = np.log2(pfm / background)
        return pwm

    @staticmethod
    def __pfm_to_pwm_cisbp(pfm):
        background = 0.25
        pwm = np.log2((pfm + 0.01) / background)
        return pwm

    @staticmethod
    def __pcm_to_pfm(pcm):
        number_of_sites = pcm.sum(axis=0)
        nuc_pseudo = 0.25
        pfm = (pcm + nuc_pseudo) / (number_of_sites + 1)
        return pfm
    
    @staticmethod
    def __get_number_of_sites(peaks, length):
        n = 0
        for p in peaks:
            n += len(p) - length + 1
        return n
    
    def min_score(self):
        return np.sum(self.matrix.min(axis=0))

    def max_score(self):
        return np.sum(self.matrix.max(axis=0))
    
    def to_tuple(self):
        return tuple(tuple(i) for i in self.matrix)

    def to_score(self, norm_value):
        min_s = self.min_score()
        max_s = self.max_score()  
        score = norm_value * (max_s - min_s) + min_s
        return score

    def __read_seqs_with_complement(self, path):
        container = list()
        letters = {'A', 'C', 'G', 'T'}
        with open(path) as file:
            for line in file:
                line = line.strip().upper()
                if not line.startswith('>'):
                    seq = ''.join([l if l in letters else 'N' for l in line])
                    complement_seq = self.__complement(seq)
                    container.append(seq)
                    container.append(complement_seq)
        return(container)

    @staticmethod
    def __complement(seq):
        seq = seq.replace('A', 't')
        seq = seq.replace('T', 'a')
        seq = seq.replace('C', 'g')
        seq = seq.replace('G', 'c')
        seq = seq.upper()
        seq = seq[::-1]
        return(seq)

    def calculate_scores_upper_threshold(self, peaks, threshold):
        scores = []
        matrix = [self.to_tuple()]
        threshold = [threshold]
        bg = MOODS.tools.flat_bg(4)
        for p in peaks:
            result = MOODS.scan.scan_dna(p, matrix, bg, threshold, 5)
            scores += [r2.score for r1 in result for r2 in r1]
        return scores
    
    def calculate_table(self, path):
        peaks = self.__read_seqs_with_complement(path)
        threshold = self.to_score(0.7)
        scores = self.calculate_scores_upper_threshold(peaks, threshold)
        number_of_sites = self.__get_number_of_sites(peaks, self.length)
        scores.sort(reverse=True)
        last_score = scores[0]
        count = 0
        fpr = 0.0
        table = list()
        while fpr <= 0.0005:
            count += 1
            fpr = count/number_of_sites
            score = scores[count]
            if score != last_score:
                table.append((last_score, fpr))
                last_score = score 
        self.__table = table
        
    def choose_threshold(self, fpr):
        last_score, last_fpr = self.__table[0]
        for line in self.__table[1:]:
            if line[1] > fpr:
                break
            else:
                last_score, last_fpr = line
        return(last_score)
    
    def __repr__(self):
        return str(self.matrix)

    
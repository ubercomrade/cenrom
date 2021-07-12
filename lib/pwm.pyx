from random import sample, shuffle
from copy import deepcopy
import sys
import itertools
import math


cdef class PWM:
    
    cdef:
        dict __matrix
        int __length
        list __table
    
    def __init__(self, path, form):
        if form == 'hocomoco':
            self.__matrix = PWM.__read_matrix(path, 0)
            self.__matrix = PWM.__pcm_to_pfm(self.__matrix)
            self.__matrix = PWM.__pfm_to_pwm(self.__matrix)
        if form == 'cisbp':
            self.__matrix = PWM.__read_matrix(path, 1)
            self.__matrix = PWM.__pfm_to_pwm(self.__matrix)
        if form == 'homer':
            self.__matrix = PWM.__read_matrix(path, 0)
            self.__matrix = PWM.__pfm_to_pwm(self.__matrix)
        if form == 'pwm':
            self.__matrix = PWM.__read_matrix(path, 0)
        self.__length = len(self.__matrix['A'])
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
        matrix = {'A':[], 'C':[], 'G':[], 'T':[]}
        with open(path) as file:
            file.readline()
            for line in file:
                line = line.strip().split('\t')[pos:]
                for letter, value in zip(matrix.keys(), line):
                    matrix[letter].append(float(value))
        return matrix
    
    @staticmethod
    def __pfm_to_pwm(pfm):
        pwm = {}
        background = 0.25
        mono_nucleotides = ['A', 'C', 'G', 'T']
        for i in mono_nucleotides:
            pwm[i[0]] = []
        first_key = list(pfm.keys())[0]
        for i in range(len(pfm[first_key])):
            for j in pfm.keys():
                pwm[j].append(math.log(pfm[j][i] / background))
        return pwm
    
    @staticmethod
    def __pcm_to_pfm(pcm):
        matrix_length = len(pcm['A'])
        number_of_sites = [0] * len(pcm['A'])
        for key in pcm.keys():
            for i in range(len(pcm[key])):
                number_of_sites[i] += pcm[key][i]
        pfm = {'A':[], 'C':[], 'G':[], 'T':[]}
        nuc_pseudo = 0.25
        for i in range(matrix_length):
            for nuc in pcm.keys():
                pfm[nuc].append((pcm[nuc][i] + nuc_pseudo) / (number_of_sites[i] + 1))
        return pfm
    
    cdef float __score_pwm(self, str seq, dict pwm):
        cdef float value 
        cdef float score = 0
        cdef int position = 0 
        cdef int length, index
        cdef str letter
        cdef list values
        length = len(seq)
        for index in range(length):
            letter = seq[index]
            values = self.__matrix[letter]
            value = values[index]
            score += value
        return score

    cdef float min_score(self):
        cdef float v, value = 0.0
        cdef list keys
        value = int()
        keys = list(self.__matrix.keys())
        for i in range(self.__length):
            tmp = []
            for j in keys:
                tmp.append(self.__matrix[j][i])
                v = min(tmp)
            value += v
        return value

    cdef float max_score(self):
        cdef float v, value = 0.0
        cdef list keys
        value = int()
        keys = list(self.__matrix.keys())
        for i in range(self.__length):
            tmp = []
            for j in keys:
                tmp.append(self.__matrix[j][i])
                v = max(tmp)
            value += v
        return value

    cdef float __to_score(self, float norm_value):
        cdef float min_s, max_s
        min_s = self.min_score()
        max_s = self.max_score()  
        score = norm_value * (max_s - min_s) + min_s
        return score

    def __read_seqs_with_complement(self, str path):
        cdef str seq, line, l
        cdef list container
        cdef set letters
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
    def __complement(str seq):
        seq = seq.replace('A', 't')
        seq = seq.replace('T', 'a')
        seq = seq.replace('C', 'g')
        seq = seq.replace('G', 'c')
        seq = seq.upper()
        seq = seq[::-1]
        return(seq)

    def calculate_scores_upper_threshold(self, list peaks, float threshold):
        cdef list scores
        cdef str peak, site
        cdef int number_of_sites, number_of_peaks, index, N, i
        cdef float score
        scores = list()
        number_of_sites = 0
        number_of_peaks = len(peaks)
        for index in range(number_of_peaks):
            peak = peaks[index]
            N = len(peak) - self.__length + 1
            for i in range(N):
                site = peak[i:self.__length + i]
                if 'N' in site:
                    continue
                number_of_sites += 1
                score = self.__score_pwm(site, self.__matrix)
                if score >= threshold:
                    scores.append(score)
        return scores, number_of_sites
    
    def calculate_table(self, str path):
        cdef list peaks, scores
        cdef float threshold, score, last_score, fpr
        cdef int number_of_sites, count
        peaks = self.__read_seqs_with_complement(path)
        threshold = self.__to_score(0.5)
        scores, number_of_sites = self.calculate_scores_upper_threshold(peaks, threshold)
        scores.sort(reverse=True)
        last_score = scores[0]
        count = 0
        fpr = 0.0
        table = list()
        while fpr <= 0.0005:
            count += 1
            fpr = float(count)/number_of_sites
            score = scores[count]
            if score != last_score:
                table.append((last_score, fpr))
                last_score = score 
        self.__table = table
        
    def choose_threshold(self, float fpr):
        cdef float last_score, last_fpr
        last_score, last_fpr = self.__table[0]
        for line in self.__table[1:]:
            if line[1] > fpr:
                break
            else:
                last_score, last_fpr = line
        return(last_score)
    
    def __repr__(self):
        return str(self.matrix)

    
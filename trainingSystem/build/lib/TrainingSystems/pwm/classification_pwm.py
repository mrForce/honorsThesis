from collections import Counter
from TrainingSystems import AbstractBinaryClassifier
import math
from scipy.stats import binom
class UnalignedSequencesError(Exception):
    def __init__(self):
        self.message = "Not all of the sequences were of the same length"
        

class ClassificationPWM:
    """
    positive_dataset and negative_dataset are both lists of sequences all of the same length.

    Only pass in background_dataset if you want to. It should be a list of sequences, each of the same length as the sequences in the positive and negative datasets

    alphabet should be a list of symbols.
    """
    def __init__(self, positive_dataset, negative_dataset, pseudocount_value, alphabet, background_dataset = False, motif_x_scoring = False):
        positive_pwm = dict()
        negative_pwm = dict()


        if background_dataset: 
            background_pfm = ClassificationPWM.computeBackgroundPFM(background_dataset, alphabet)
            if motif_x_scoring == False:
                self.positive_pwm = ClassificationPWM.computePWM(positive_dataset, pseudocount_value, alphabet, background_pfm)
                self.negative_pwm = ClassificationPWM.computePWM(negative_dataset, pseudocount_value, alphabet, background_pfm)
            else:
                self.positive_pwm = ClassificationPWM.computeMotifXPWM(positive_dataset, alphabet, background_pfm)
                if len(negative_dataset) > 0:
                    self.negative_pwm = ClassificationPWM.computeMotifXPWM(negative_dataset, alphabet, background_pfm)


        else:
            self.positive_pwm = ClassificationPWM.computePWM(positive_dataset, pseudocount_value, alphabet)
            self.negative_pwm = ClassificationPWM.computePWM(negative_dataset, pseudocount_value, alphabet)





    """
    I don't use pseudocounts here, since I assume that there will be enough sequences in the background such that every letter is represented at every location.
    """
    @staticmethod
    def computeBackgroundPFM(background_sequences, alphabet):
        sequence_width = len(background_sequences[0])
        num_sequences = len(background_sequences)
        frequencies = list()
        for i in range(0, sequence_width):
            for x in background_sequences:
                if len(x) <= i:
                    print('sequence: ' + x)
            column = list([sequence[i] for sequence in filter(len, background_sequences)])
            counts = Counter(column)
            column_frequencies = dict()
            for amino_acid, count in counts.items():
                column_frequencies[amino_acid] = 1.0*count/num_sequences
            frequencies.append(column_frequencies)
        return frequencies
                



    """
    This function should take in a list of sequences and a pseudocount value, and return the PWM as a dictionary that maps an amino acid to a list of weights, specifying the weights at each position. 

    Unlike the other function, this NEEDS a background pfm.
"""
    @staticmethod
    def computeMotifXPWM(dataset, alphabet, background_pfm):
        sequence_length = len(dataset[0])
        for i in range(1, len(dataset)):
            if len(dataset[i]) != sequence_length:
                raise UnalignedSequencesError
        num_sequences = len(dataset)
        
        pwm = dict()
        for amino_acid in alphabet:
            pwm[amino_acid] = list()
        num_amino_acids = len(alphabet)
        
        for position in range(0, sequence_length):
            column = list([sequence[position] for sequence in dataset])
            counts = Counter(column)
            for amino_acid in alphabet:
                #if that column doesn't contain the given amino acid, then counts[amino_acid] simply returns 0, as it should 
                count = counts[amino_acid]
                if amino_acid in background_pfm[position]:
                    background_probability = background_pfm[position][amino_acid]
                    pUnder = binom.cdf(count, num_sequences, background_pfm[position][amino_acid])
                    pOver = 1.0 - pUnder + binom.pmf(count, num_sequences, background_pfm[position][amino_acid])
                    pwm[amino_acid].append(math.log2(pUnder/pOver))
                else:
                    pwm[amino_acid].append(0.0)
        return pwm

    """
    This function should take in a list of sequences and a pseudocount value, and return the PWM as a dictionary that maps an amino acid to a list of weights, specifying the weights at each position. """
    @staticmethod
    def computePWM(dataset, pseudocount_value, alphabet, background_pfm = False):
        sequence_length = len(dataset[0])
        for i in range(1, len(dataset)):
            if len(dataset[i]) != sequence_length:
                raise UnalignedSequencesError
        num_sequences = len(dataset)
        
        pwm = dict()
        for amino_acid in alphabet:
            pwm[amino_acid] = list()
        num_amino_acids = len(alphabet)
        
        for position in range(0, sequence_length):
            column = list([sequence[position] for sequence in dataset])
            counts = Counter(column)
            for amino_acid in alphabet:
                #if that column doesn't contain the given amino acid, then counts[amino_acid] simply returns 0, as it should 
                count = counts[amino_acid]

                if background_pfm:
                    pwm[amino_acid].append(math.log2((count + pseudocount_value)/(background_pfm[position][amino_acid]*(num_sequences + num_amino_acids*pseudocount_value))))
                else:
                    #in the numerator, we need to divide the pseudocount value by the number of amino acids.
                    pwm[amino_acid].append(math.log2((count + pseudocount_value)/(num_sequences + pseudocount_value*num_amino_acids)))
        return pwm

    """ This returns a dictionary that maps an amino acid to its background frequency """
    @staticmethod
    def calculateBackgroundFreqs(dataset):
        count = Counter()
        for sequence in dataset:
            for x in sequence:
                count[x] += 1
        denom = 1.0*sum(count.values())
        count_dict = dict(count)
        for x in count_dict.keys():
            count_dict[x] = count_dict[x]/denom
        return count_dict

    @staticmethod
    def computeScore(pwm, sequence):
        return sum([pwm[sequence[i]][i] for i in range(0, len(sequence))])


    def computePositiveScore(self, sequence):
        return ClassificationPWM.computeScore(self.positive_pwm, sequence)
    def computeNegativeScore(self, sequence):
        return ClassificationPWM.computeScore(self.negative_pwm, sequence)
class PWMBinaryClassifier(AbstractBinaryClassifier):

    #just pass it the background
    def __init__(self, background, motifx_scoring = False):
        self.background = background
        self.motifx_scoring = motifx_scoring
        
    def train(self, dataset, alphabet):
        binding_dataset = list(map(lambda x: x[0], filter(lambda x: x[1] == True, dataset)))
        nonbinding_dataset = list(map(lambda x: x[0], filter(lambda x: x[1] == False, dataset)))
        
        self.classifier = ClassificationPWM(binding_dataset, nonbinding_dataset, 0.01, alphabet, self.background, self.motifx_scoring)

    def classify(self, peptide):
        if self.classifier.computePositiveScore(peptide) > self.classifier.computeNegativeScore(peptide):
            return True
        else:
            return False

    def get_positive_score(self, peptide):
        return self.classifier.computePositiveScore(peptide)

from collections import Counter
import sys
from TrainingSystems import AbstractBinaryClassifier
from TrainingSystems import classification_pwm
import math
import pdb
from scipy.stats import binom
import numpy as np
import itertools

from hmmlearn.hmm import MultinomialHMM
class UnalignedSequencesError(Exception):
    def __init__(self):
        self.message = "Not all of the sequences were of the same length"
        

class ClassificationHMM:
    """
    I don't think it makes any sense to use a negative HMM. After all, negatives are all the sequences that simply do not bind well to the MHC. 

    Only pass in background_dataset if you want to. It should be a list of sequences, each of the same length as the sequences in the positive and negative datasets

    alphabet should be a list of symbols.

    num_rows should be an integer that specifies the number of rows we want.
    """
    def __init__(self, positive_dataset, alphabet, num_rows, background_dataset = False):
        positive_pwm = dict()


        self.alphabet = alphabet
        self.amino_acids_in_dataset = list(set(itertools.chain(*positive_dataset)))
        self.residue_mapper = dict(zip(self.amino_acids_in_dataset, range(0, len(self.amino_acids_in_dataset))))
        i = len(self.amino_acids_in_dataset)
        for amino_acid in self.alphabet:
            if amino_acid not in self.residue_mapper:
                self.residue_mapper[amino_acid] = i
                i += 1
        
        if background_dataset:
            #Yes, we do want to use the background sequences. 
            self.background_pfm = classification_pwm.ClassificationPWM.computeBackgroundPFM(background_dataset, alphabet)
            self.hmm = ClassificationHMM.computeHMM(positive_dataset,  alphabet, self.amino_acids_in_dataset, self.residue_mapper, num_rows)
            self.using_background = True
        else:
            self.hmm = ClassificationHMM.computeHMM(positive_dataset, alphabet, self.amino_acids_in_dataset, self.residue_mapper, num_rows)
            self.using_background = False
    """
            Generate n floating points that add up to 1.0
    """
    @staticmethod
    def compute_random_row(num_elements):
        elements = np.random.uniform(0.0, 1.0, num_elements)
        total = sum(elements)
        return elements/total
    """
    Computes haplotype-ish HMM. 
    
    Can we integrate the background PFM into this?
    """
    @staticmethod
    def computeHMM(dataset, alphabet, amino_acids_in_dataset, residue_mapper, num_rows):
        sequence_length = len(dataset[0])
        for i in range(1, len(dataset)):
            if len(dataset[i]) != sequence_length:
                raise UnalignedSequencesError
        num_sequences = len(dataset)

        
        dataset = np.concatenate([[[residue_mapper[x]] for x in y] for y in dataset])
        print('numbers in dataset')
        print(set(dataset.flatten().tolist()))
#        pdb.set_trace()

        num_states = sequence_length*num_rows
        transition_matrix = np.zeros((num_states, num_states))
        for col in range(0, sequence_length - 1):
            for row in range(0, num_rows):
                rand_row = ClassificationHMM.compute_random_row(num_rows)
                transition_matrix[col*num_rows + row][num_rows*(col + 1):(num_rows*(col + 2))] = rand_row
        col = sequence_length - 1

        for row in range(0, num_rows):
            transition_matrix[col*num_rows + row][col*num_rows + row] = 1.0
        print('initial transition matrix')
        print(transition_matrix)
        model = MultinomialHMM(num_states, params="ets", init_params="e")
        model.n_features = len(alphabet)
        start_prob = np.zeros(num_states)
        start_prob[0:num_rows] = ClassificationHMM.compute_random_row(num_rows)
        print('start prob array')
        print(start_prob)
        model.startprob_ = start_prob
        model.transmat_ = transition_matrix
        try:
            model.fit(dataset, [sequence_length]*num_sequences)
        except ValueError:
            pdb.set_trace()
        print('model')
        print(model)
        for row in range(0, num_rows):
            model.transmat_[col*num_rows + row][col*num_rows + row] = 1.0
        print('transition matrix')
        print(model.transmat_)
        print('start prob')
        print(model.startprob_)
        print('amino acids in dataset')
        print(amino_acids_in_dataset)
        print('alphabet')
        print(alphabet)
        """        if len(alphabet) > len(amino_acids_in_dataset):
            num_missing_aminos = len(alphabet) - len(amino_acids_in_dataset)
            print('emissions before concat')
            print(model.emissionprob_)
            model.emissionprob_ = np.concatenate((model.emissionprob_, np.zeros((len(model.emissionprob_), num_missing_aminos))), axis=1)
            print('emissions after concat')
            print(model.emissionprob_)"""
        for row in range(0, len(model.emissionprob_)):
            for col in range(0, len(model.emissionprob_[row])):
                count = model.emissionprob_[row][col]*num_sequences
                model.emissionprob_[row][col] = (count + 0.01)/(num_sequences + len(alphabet)*0.01)
        print('emission probabilities')
        print(model.emissionprob_)
        return model

    def computePositiveScore(self, peptide):
        peptide_mapped = [[self.residue_mapper[x]] for x in peptide]
        #they return log probabilities.
        score = self.hmm.score(peptide_mapped, [len(peptide)])
        print('score')
        print(score)
        if self.using_background:
            background_probability = 0.0
            for i in range(0, len(peptide)):
                background_probability += math.log(self.background_pfm[i][peptide[i]])
            
            return score - background_probability
        else:
            return score

class HMMBinaryClassifier(AbstractBinaryClassifier):

    #just pass it the background
    def __init__(self, num_rows, background=False):
        self.num_rows = num_rows
        self.background = background
        
    def train(self, dataset, alphabet):
        binding_dataset = list(map(lambda x: x[0], filter(lambda x: x[1] == True, dataset)))
        self.classifier = ClassificationHMM(binding_dataset, alphabet, self.num_rows, self.background)

    def classify(self, peptide):
        return False


    def get_positive_score(self, peptide):
        return self.classifier.computePositiveScore(peptide)


#h = ClassificationHMM(['ABC', 'AAA', 'ACC', 'CCA', 'BBB', 'ACB'], 'ABC', 2)
#print(h)

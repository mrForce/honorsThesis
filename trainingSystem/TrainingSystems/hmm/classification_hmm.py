from collections import Counter
from TrainingSystems import AbstractBinaryClassifier
from TrainingSystems import classification_pwm
import math
from scipy.stats import binom
import numpy as np
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

        if background_dataset:
            #Yes, we do want to use the background sequences. 
            background_pfm = ClassificationPWM.computeBackgroundPFM(background_dataset, alphabet)
            self.hmm = ClassificationHMM.computeHMM(positive_dataset,  alphabet, num_rows, background_pfm)
        else:
            self.hmm = ClassificationHMM.computeHMM(positive_dataset, alphabet, num_rows)
                
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
    def computeHMM(dataset, alphabet, num_rows, background_pfm=False):
        sequence_length = len(dataset[0])
        for i in range(1, len(dataset)):
            if len(dataset[i]) != sequence_length:
                raise UnalignedSequencesError
        num_sequences = len(dataset)

        residue_mapper = dict(zip(alphabet, range(0, len(alphabet))))
        dataset = [[[residue_mapper[x]] for x in y] for y in dataset]
        #a start state, end state, and num_rows*sequence_length states in between
        num_states = 2 + sequence_length*num_rows
        transition_matrix = np.zeros((num_states, num_states))
        #set the transition probabilites coming out of the start state.
        transition_matrix[0][1:(num_rows + 1)] = ClassificationHMM.compute_random_row(num_rows)
        """
        To be clear, when I talk about the rows and columns, I'm
        talking about the structure of the HMM -- NOT the rows and columns of the transition matrix.
        """
        for column_number in range(0, sequence_length - 1):
            #the index of the first state in the next column.
            first_transition_state = 1 + (1 + column_number)*num_rows
            #technically, this is the index of the first state two columns to the right -- but that's exactly what we want!
            last_transition_state = first_transition_state + num_rows
            for row_number in range(0, num_rows):
                state_number = 1 + column_number*num_rows + row_number
                transition_matrix[state_number][first_transition_state:last_transition_state] = ClassificationHMM.compute_random_row(num_rows)

        column = sequence_length - 1
        for state in range(column*num_rows + 1, column*num_rows + 1 + num_rows):
            transition_matrix[state][-1] = 1.0/num_rows
        transition_matrix[-1][-1] = 1.0
        initial_probs = np.zeros(num_states)
        initial_probs[0] = 1.0
        model = MultinomialHMM(num_states, params="et", init_params="e")
        model.fit(np.concatenate(dataset), [sequence_length]*num_sequences)
        return model

    def compute_positive_score(peptide):
        residue_mapper = dict(zip(alphabet, range(0, len(alphabet))))
        peptide_mapped = [[residue_mapper[x]] for x in peptide]
        return model.score(peptide_mapped)

class HMMBinaryClassifier(AbstractBinaryClassifier):

    #just pass it the background
    def __init__(self, num_rows):
        self.num_rows
        
    def train(self, dataset, alphabet):
        binding_dataset = list(map(lambda x: x[0], filter(lambda x: x[1] == True, dataset)))
        self.classifier = ClassificationHMM(binding_dataset, alphabet, self.num_rows)

    def classify(self, peptide):
        return False


    def get_positive_score(self, peptide):
        return self.classifier.computePositiveScore(peptide)


#h = ClassificationHMM(['ABC', 'AAA', 'ACC', 'CCA', 'BBB', 'ACB'], 'ABC', 2)
#print(h)

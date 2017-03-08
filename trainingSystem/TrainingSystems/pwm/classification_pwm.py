from collections import Counter
from TrainingSystems import AbstractBinaryClassifier

class UnalignedSequencesError(Exception):
    def __init__(self):
        self.message = "Not all of the sequences were of the same length"
        

class ClassificationPWM:
    """
    positive_dataset and negative_dataset are both lists of sequences all of the same length.

    Only pass in background_dataset if you want to.

    alphabet should be a list of symbols.
    """
    def __init__(self, positive_dataset, negative_dataset, pseudocount_value, alphabet, background_dataset = False):
        positive_pwm = dict()
        negative_pwm = dict()
        if background_dataset: 
            background_freqs = calculateBackgroundFreqs(background_dataset)
            self.positive_pwm = computePWM(positive_dataset, pseudocount_value, alphabet, background_freqs)
            self.negative_pwm = computePWM(negative_dataset, pseudocount_value, alphabet, background_freqs)
        else:
            self.positive_pwm = computePWM(positive_dataset, pseudocount_value, alphabet)
            self.negative_pwm = computePWM(negative_dataset, pseudocount_value, alphabet)

    
    """
    This function should take in a list of sequences and a pseudocount value, and return the PWM as a dictionary that maps an amino acid to a list of numbers. """
    @staticmethod
    def computePWM(dataset, pseudocount_value, alphabet, background_freqs = False):
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
            counts = Counter()
            for i in range(0, len(dataset)):
                counts[dataset[i][position]] += 1
            for amino_acid, count in counts.items():
                if background_freqs:
                    pwm[amino_acid].append(math.log2((count + pseudocount_value/num_amino_acids)/(background_freqs[amino_acid]*(num_sequences + pseudocount_value))))
                else:
                    #in the numerator, we need to divide the pseudocount value by the number of amino acids.
                    pwm[amino_acid].append(math.log2((count + pseudocount_value/num_amino_acids)/(num_sequences + pseudocount_value)))
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
        return sum(map(lambda x: pwm[x], sequence))

    def computePositiveScore(sequence):
        return computeScore(self.positive_pwm, sequence)
    def computeNegativeScore(sequence):
        return computeScore(self.negative_pwm, sequence)
class PWMBinaryClassifier(AbstractBinaryClassifier):

    #just pass it the background
    def __init__(self, background):
        self.background = background
        
    def train(dataset, alphabet):
        binding_dataset = list(map(lambda x: x[0], filter(lambda x: x[1] == True, dataset)))
        nonbinding_dataset = list(map(lambda x: x[0], filter(lambda x: x[1] == False, dataset)))
        self.classifier = ClassificationPWM(binding_dataset, nonbinding_dataset, 0.01, alphabet, self.background)

    def classify(peptide):
        if self.classifier.computePositiveScore(sequence) > self.classifier.computeNegativeScore(sequence):
            return True
        else:
            return False

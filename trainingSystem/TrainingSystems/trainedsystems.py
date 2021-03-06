"""
For a given MHC allele, this holds a dictionary of binary and affinity predictors.

It also holds the training data for reference
"""


class NeedTrainingDataException(Exception): pass

class TrainedSystems:
    mhc_name = ''
    """
    binary_training is a list of tuples, of the form [(peptide, yes/no)...]

    affinity_training is a list of tuples, of the form [(peptide, Kd)...]

    At least one of these needs to be specified.
    """
    def __init__(self, mhc_name, binary_training = False, affinity_training = False):
        if binary_training == False and affinity_training == False:
            raise NeedTrainingDataException
        self.mhc_name = mhc_name
        self.binary_training = binary_training
        self.affinity_training = affinity_training
        self.binary_classifiers = dict()
        self.affinity_predictors = dict()

    """
    Need to pass in an instance of BinaryClassifier (from interfaces.py)
"""
    def trainBinary(self, classifier_name, classifier_object):
        classifier_object.train(self.binary_training)
        self.binary_classifiers[classifier_name] = classifier_object

    def trainAffinity(self, predictor_name, predictor_object):
        predictor_object.train(self.affinity_training)
        self.affinity_predictors[predictor_name] = predictor_object
        
    """
    test_data is a list of peptides 
    
    This function takes each classifier in self.binary_classifiers,
    and calls the classify function on each, and returns the dict mapping the classifier name to the result.
    """
    def classify(self, test_data):
        return {name: {peptide: classifier.classify(peptide) for peptide in test_data} for name, classifier in self.binary_classifiers.items()}
    
    """
    test_data is a list of peptides 
    
    This function takes each classifier in self.binary_classifiers,
    and calls the classify function on each, and returns the dict mapping the classifier name to the result.
    """
    def predict_affinities(self, test_data):
        return {name: {peptide: predictor.predict_affinity(peptide) for peptide in test_data} for name, predictor in self.affinity_predictors.items()}

from abc import ABCMeta, abstractmethod

class AbstractBinaryClassifier:
    __metaclass__ = ABCMeta


    """
    This returns True if the peptide is a binder, returns False if it is not.
"""
    @abstractmethod
    def classify(peptide): pass

    """
    Pass in training data as a list of tuples of the form [(peptide, yes/no)...]
    """
    @abstractmethod
    def train(data, alphabet): pass

class AbstractAffinityPredictor:
    __metaclass__ = ABCMeta

    """ Predicts the Kd, which should be a double """
    @abstractmethod
    def predict_affinity(peptide): pass

    """
    Pass in training data as a list of tuples of the form [(peptide, Kd)...]
    """
    @abstractmethod
    def train(data, alphabet): pass

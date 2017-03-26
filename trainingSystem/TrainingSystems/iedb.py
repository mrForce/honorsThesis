import sys
import pickle
import argparse

class Error(Exception):
    pass

class HLAAlreadyDefinedError(Error):
    pass

class IEDBData:
    def __init__(self, iedb_filters):
        self.iedb_filters = iedb_filters
        self.data = dict()
    def get_iedb_filters(self):
        return self.iedb_filters
    #values is a list of tuples of form [(peptide, Kd)...]
    def add_data(hla_name, values):
        if hla_name in self.data:
            raise HLAAlreadyDefinedError
        else:
            self.data[hla_name] = values


parser = argparse.ArgumentParser(description='Manipulate the storage of HLA allele data')
#we have the positional argument, which is the pickle that we are referring to
parser.add_argument('dataPickle', help='Give us the name of the pickle we are working with. If it doesn\'t exist, then I will create it.')
parser.add_argument('--listHLA', help='Lists the HLA alleles that are already in the pickle (default: False)', action='store_true')
parser.add_argument('--addHLA', help='Adds the HLA with NAME and CSV data to the pickle', nargs=2, metavar=('hla_name', 'csv_file'))
parser.add_argument('--showIEDBFilters', help='Prints out the IEDB filters that were used to collect data (other than the HLA allele) (default: False)', action='store_true')
parser.add_argument('--setIEDBFilters', help='Sets the IEDB filters that were used to collect the data. Simply pass in a string with quotes around it. Do not pass the HLA allele name to this', nargs=1, metavar='IEDBFilters', default=False)
args = parser.parse_args()
print(args)
"""
if len(sys.argv) == 4:
    with open(sys.argv[3], 'r') as f:
        pickle_contents = pickle.load(f)
        hla_name  = sys.argv[1]
        assays_file_name = sys.argv[2]
        
        print('HLA Name: ' + hla_name)
        if hla_name in pickle_contents.keys():
            print('HLA is already in the file')
        else:
            
else:
    print('usage: python iedb_add.py hla_name assays_file.csv results.pickle')
"""

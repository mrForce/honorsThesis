import sys
import pickle
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

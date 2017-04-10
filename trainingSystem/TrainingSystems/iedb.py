import os.path
import csv
import sys
import pickle
import argparse
from Bio.Alphabet import IUPAC
from iedb_data import *


    
"""
This function first checks which line (of the first 3) has 'Object Type', 'Description', 'Starting Position', 'Quantitative measurement', and 'Ending Position' in it. Then in takes the chunked up line, and prints it out to the user, asking for confirmation that the line it has chosen really does contain the field names. I'm doing this because the IEDB output seems to have two lines for field names, but that doesn't make any sense. 

Then, for each line past the fields line, we make sure that the Description field (which holds the peptide sequence) is a valid peptide, and that the Object Type is Linear Peptide, and that Ending Position - Starting Position + 1 = len(peptide). If any of these conditions fail, then we ask the user whether we should discard the data. 

Pass in the parsed lines of the CSV file. It will return a list of tuples of the form [(peptide, Kd)...]

"""
def getData(parsed_lines):
    print("I am about to extract the data from the CSV file. However, I will need you to confirm a few things with me along the way")
    i = 0
    fields_index = -1
    while i < 3:
        line = parsed_lines[i]
        if 'Object Type' in line and 'Description' in line and 'Starting Position' in line and 'Ending Position' in line and 'Quantitative measurement' in line:
            fields_index = i
            break
        i += 1
    if fields_index == -1:
        sys.exit('We couldn\'t find the line with the field names in it. We searched the first 3 lines.')
    print('We selected these as the field names:')
    print(parsed_lines[fields_index])
    keepGoing  = input('Is this correct? (yes/no) ')
    if keepGoing != 'yes':
        sys.exit('You didn\'t say yes, so we are aborting')
    data = list()
    #we add two, since line numbers are indexed at 1, and we are also at the line after the fields line.
    line_number = fields_index + 2
    for x in parsed_lines[(fields_index + 1)::]:
        entry = dict(zip(parsed_lines[fields_index], x))
        if 'Description' in entry and 'Starting Position' in entry and 'Ending Position' in entry and 'Object Type' in entry:
            peptide = entry['Description']
            object_type = entry['Object Type']
            no_position = False
            try:
                starting_position = int(entry['Starting Position'])
            except ValueError:
                no_position = True
            try:
                ending_position = int(entry['Ending Position'])
            except ValueError:
                no_position = True
            position_criteria = no_position or ending_position - starting_position + 1 == len(peptide) 
            kd = float(entry['Quantitative measurement'])
            if object_type == 'Linear peptide' and all([aa in IUPAC.protein.letters for aa in peptide]) and position_criteria:
                print('added to data')
                data.append((peptide, kd))
            else:
                print('I\'m not sure if I should keep this (on line: ' + str(line_number) + '):')
                print('peptide: ' + peptide)
                print('object type: ' + object_type)
                print('starting position: ' + entry['Starting Position'])
                print('ending position: ' + entry['Ending Position'])
                print('Kd: ' + str(kd))
                
                keep = input('Should I keep this? (yes/no)')
                if keep == 'yes':
                    data.append((peptide, kd))
                elif keep != 'no':
                    keepGoing = True
                    while keepGoing:
                        keep = input('Sorry, didn\'t catch that (yes/no)')
                        if keep == 'yes':
                            data.append((peptide, kd))
                            keepGoing = False
                        elif keep == 'no':
                            keepGoing = False 

        line_number += 1
    return data


parser = argparse.ArgumentParser(description='Manipulate the storage of HLA allele data')
#we have the positional argument, which is the pickle that we are referring to
parser.add_argument('dataPickle', help='Give us the name of the pickle we are working with. If it doesn\'t exist, then I will create it.')
parser.add_argument('--listHLA', help='Lists the HLA alleles that are already in the pickle (default: False)', action='store_true')
parser.add_argument('--addHLA', help='Adds the HLA with NAME and CSV data to the pickle', nargs=2, metavar=('hla_name', 'csv_file'), default=False)
parser.add_argument('--showIEDBFilters', help='Prints out the IEDB filters that were used to collect data (other than the HLA allele) (default: False)', action='store_true')
parser.add_argument('--setIEDBFilters', help='Sets the IEDB filters that were used to collect all of the data in the pickle. You probably shouldn\'t need to call this more than once. Simply pass in a string with quotes around it. Do not pass the HLA allele name to this', metavar='IEDBFilters', default=False)
args = parser.parse_args()
print(args)

pickle_file = args.dataPickle
iedb_data_object = False
#we need to create the IEDBData object, and then at the end we will write
if os.path.exists(pickle_file):
    with open(pickle_file, 'rb') as f:
        iedb_data_object = pickle.load(f)
    if iedb_data_object == False:
        sys.exit('Couldn\'t load the pickle object. Aborting.')
else:
    iedb_data_object = IEDBData()

if args.listHLA:
    alleles = iedb_data_object.get_data().keys()
    if len(alleles)  == 0:
        print('No HLA alleles have been added')
    else:
        print('Here are the HLA alleles:')
        for allele in alleles:
            print(allele)
if args.showIEDBFilters:
    filters = iedb_data_object.get_iedb_filters()
    if filters:
        print('The filters: ' + filters)
    else:
        print('IEDB Filters haven\'t been set yet')

if args.setIEDBFilters:
    filters = iedb_data_object.get_iedb_filters()
    if filters:
        print('These are the current filters: ' + filters)
        change_filter = input('Would you like to change this? (yes/no)')
        if change_filter == 'yes':
            iedb_data_object.set_iedb_filters(args.setIEDBFilters)
        else:
            print('You didn\'t input yes, so we aren\'t going to change the filters')
    else:
        iedb_data_object.set_iedb_filters(args.setIEDBFilters)
        print('setting the IEDB filters to: ' + args.setIEDBFilters)

if args.addHLA:
    hla = args.addHLA[0]
    if hla in iedb_data_object.get_data():
        keepGoing = input('We already have HLA data for that object. Do you want to replace the old with the new? (yes/no)')
        if keepGoing != 'yes':
            sys.exit('exiting')
    file_name = args.addHLA[1]
    with open(file_name, 'r', encoding='utf-8') as f:
        reader = csv.reader(f)
        lines = list(reader)
        data = getData(lines)
        iedb_data_object.add_data(hla, data)
        with open(pickle_file, 'wb') as g:
            pickle.dump(iedb_data_object, g)
            print('added the HLA')

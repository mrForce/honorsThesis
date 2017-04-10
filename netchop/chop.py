import sys
import os
import time
from Bio import SeqIO
import tempfile
from collections import Counter
import progressbar
import subprocess
import re
import argparse
from Bio.Alphabet import IUPAC
"""
ERAP1 and 2 trim from left to right

"""

regex = re.compile('(\d+)\s+\w\s+(\d+\.\d+)\s+.+')

def get_pieces(protein, chunk_width, netchop_output):
#    print('netchop output')
#    print(netchop_output)
    pieces = list()
#    print('protein')
#    print(protein)
    for x in netchop_output.split('\n'):
        m = regex.match(x)
        if m != None:
            location = int(m.group(1))
            score = float(m.group(2))
            if score > 0.5 and location >= chunk_width:
#                print('location')
 #               print(location)
                chunk = protein[(location - chunk_width):location]
                pieces.append(chunk)


    return pieces


def background_complete(chunks):
    width = len(chunks[0])
    for i in range(0, width):
        amino_acids = list([x[i] for x in chunks])
        count = Counter(amino_acids)
        for x in IUPAC.protein.letters:
            if count[x] == 0:
                return False
    return True
parser = argparse.ArgumentParser(description='Run netchop on a bunch of sequences')
parser.add_argument('outputFile', metavar='outputFile.txt', help='The location where the chopped pieces will be stored.')
parser.add_argument('chunkWidth', help='The width of each chunk', type=int)
parser.add_argument('--sequenceFiles', nargs='+', help='Specify the files with sequences in them')
parser.add_argument('--checkBackgroundComplete', action='store_true', help='This checks that at every position, there is at least one of every amino acid (default: False')
args = parser.parse_args()
chunk_width = args.chunkWidth
fasta_files = args.sequenceFiles
chunks = list()

sleep_time = 0.1
with open(args.outputFile, 'w') as f:
    jobs = list()

    max_num_jobs = 4
    num_sequences = 0
    for fasta_file in fasta_files:
        num_sequences += len(list(SeqIO.parse(fasta_file, "fasta")))

    bar = progressbar.ProgressBar(max_value=num_sequences)
    j = 0


    file_free = [True for x in range(0, max_num_jobs)]
    for fasta_file in fasta_files:
        for seq_record in SeqIO.parse(fasta_file, "fasta"):
#            print('j: ' + str(j))
 #           print('num jobs: ' + str(len(jobs)))
  #          print('num sequences: ' + str(num_sequences))
   #         print('keep going true: ' + str((len(jobs) == max_num_jobs) or (len(jobs) > 0 and j + 3 + len(jobs) >= num_sequences)))
            tf = tempfile.NamedTemporaryFile('w+')
            tf.write('>' + str(seq_record.id) + '\n')
            sequence = str(seq_record.seq)
            tf.write(sequence)
            tf.flush()

            output_file = tempfile.NamedTemporaryFile()
            #if we use stdout=subprocess.PIPE, the pipe tends to fill up.
            job = (subprocess.Popen(['netchop -m netchop -n --threshold 0.5 ' + tf.name], stdout=output_file, encoding='utf-8', close_fds=True, shell=True), tf, output_file, str(seq_record.seq))
            jobs.append(job)

            while (len(jobs) == max_num_jobs) or (len(jobs) > 0 and j + len(jobs) >= num_sequences):
                time.sleep(sleep_time)
                i = 0
                indices_to_delete = list()
                
                for i in range(0, len(jobs)):
                    poll_result = jobs[i][0].poll()
                    if poll_result != None:
                        jobs[i][2].seek(0)
                        #input(jobs[i][2].name)
                        output = jobs[i][2].read().decode('utf-8')
#                        print('output')
#                        print(output)
#                        print('output filename')
#                        input(jobs[i][2].name)

                        pieces = get_pieces(jobs[i][3], chunk_width, output)
#                        print('pieces')
#                        print(pieces)
                        for x in pieces:
                            f.write(x + '\n')
                        if args.checkBackgroundComplete:
                            chunks.extend(pieces)
                        j += 1
                    #    print('j: ' + str(j))
                        bar.update(j)
                        indices_to_delete.append(i)
                        jobs[i][1].close()
                        jobs[i][2].close()
                        
                for i in sorted(indices_to_delete, reverse=True):
                    del(jobs[i])
            
            
    print('jobs remaining')
    print(jobs)
    """    for x in input_files:
        x.close()
    for x in output_files:
        x.close()"""
    print('chunks')

    if args.checkBackgroundComplete:
        if background_complete(chunks):
            print('The background is complete')
        else:
            print('The background is not complete')


    
    

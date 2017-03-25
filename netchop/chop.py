import sys
from Bio import SeqIO
import tempfile
import subprocess
import re
"""
ERAP1 and 2 trim from left to right

"""

regex = re.compile('(\d+)\s+\w\s+(\d+\.\d+)\s+.+')

def get_pieces(protein, chunk_width, netchop_output):
    pieces = list()
    for x in netchop_output.split('\n'):
        m = regex.match(x)
        if m != None:
            location = int(m.group(1))
            score = float(m.group(2))
            if score > 0.5 and location >= chunk_width:
                chunk = protein[(location - chunk_width):location]
                pieces.append(chunk)


    return pieces

if len(sys.argv) >= 4:
    output_file = sys.argv[1]
    
    chunk_width = int(sys.argv[2])
    with open(output_file, 'w') as f:
        for fasta_file in sys.argv[3::]:
            for seq_record in SeqIO.parse(fasta_file, "fasta"):
                tf = tempfile.NamedTemporaryFile('w+')
                tf.write('>' + str(seq_record.id) + '\n')
                tf.write(str(seq_record.seq))
                tf.flush()
                output = str(subprocess.run(['netchop', '-m', 'netchop', '--threshold', '0.5', tf.name], stdout=subprocess.PIPE, encoding='utf-8').stdout)
                pieces = get_pieces(str(seq_record.seq), chunk_width, output)
                for x in pieces:
                    f.write(x + '\n')
else:
    print('usage: python chop.py output.txt chunk_width proteome1.fasta proteome2.fasta...')
            

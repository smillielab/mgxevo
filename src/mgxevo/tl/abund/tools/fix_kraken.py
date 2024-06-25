import sys

imap = {
    'R':  'R',
    'R1': 'D',
    'R2': 'P',
    'R3': 'C',
    'R4': 'O',
    'R5': 'F',
    'R6': 'G',
    'R7': 'S'
}

for line in open(sys.argv[1]):
    line = line.rstrip().split('\t')
    if line[3] in imap:
        line[3] = imap[line[3]]
    print('\t'.join(line))

import argparse, os, re, sys
from pathlib import Path
sys.path.append(str(Path(__file__).resolve().parent.parent.parent.parent))

from ut import seq as seq_util

# input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--aln', help='fasta alignment')
parser.add_argument('--lax', help='relaxed auto settings', default=False, action='store_true')
parser.add_argument('--b1', help='min %% of sequences for a conserved position', default=.5, type=float)
parser.add_argument('--b2', help='min %% of sequences for a flank position', default=.85, type=float)
parser.add_argument('--b3', help='max # of contiguous non-conserved positions', default=8, type=int)
parser.add_argument('--b4', help='min length of block', default=10, type=int)
parser.add_argument('--b5', help='allow gap positions', choices=['none', 'half', 'all'], default='none')
parser.add_argument('--cut', help='minimum un-gapped length per sequence', default=0, type=int)
parser.add_argument('-p', help='print command', default=False, action='store_true')
args = parser.parse_args()

# auto settings
if args.lax == True:
    args.b1 = .5
    args.b2 = .5
    args.b3 = 10
    args.b4 = 5
    args.b5 = 'all'
    args.cut = 25

# fix arguments
n = len([record for record in seq_util.iter_fst(args.aln)])
args.b1 = int(args.b1*n) + 1
args.b2 = int(args.b2*n) + 1
args.b5 = {'none':'N', 'half':'H', 'all':'A'}[args.b5]

# encode gblocks input
gmap = {}
temp = open('%s.temp.aln' %(args.aln), 'w')
print('writing alignmnet to %s.temp.aln' %(args.aln))
for i, record in enumerate(seq_util.iter_fst(args.aln)):
    old = record[0].strip('>')
    new = 'seq%d' %(i)
    gmap[new] = old
    record[0] = '>%s' %(new)
    temp.write('\n'.join(record) + '\n')
temp.close()

# run gblocks
cmd = 'Gblocks %s.temp.aln -b1=%d -b2=%d -b3=%d -b4=%d -b5="%s" -e=".gb"' %(args.aln, args.b1, args.b2, args.b3, args.b4, args.b5)
os.system(cmd)
print(cmd)
# decode gblocks output
rmv = tot = 0
out = open('%s.gb' %(args.aln), 'w')
for i, record in enumerate(seq_util.iter_fst('%s.temp.aln.gb' %(args.aln))):
    new = record[0].strip('>')
    old = gmap[new]
    record[0] = '>%s' %(old)
    record[1] = re.sub(' ', '', record[1])
    if len([xi for xi in record[1] if xi != '-']) < args.cut:
        rmv += 1
        continue
    out.write('\n'.join(record) + '\n')
    tot += 1
out.close()

# removed genes
print('Removed %d / %d sequences (%.2f%%)' %(rmv, tot, 100*rmv/tot))

# cleanup
cmd = 'rm %s.temp.aln; rm %s.temp.aln.gb; rm %s.temp.aln.gb.htm' %(args.aln, args.aln, args.aln)
os.system(cmd)

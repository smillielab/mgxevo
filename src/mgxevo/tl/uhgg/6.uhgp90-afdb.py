import sys
from pathlib import Path
sys.path.append(str(Path(__file__).resolve().parent.parent.parent))

import ut

fn = sys.argv[1]
dbfaa = sys.argv[2]

v = ut.read_fst(dbfaa)
u = {}
for record in ut.iter_seq(fn):
    sq = record[1].upper()
    nm = record[0][1:].split(' ')[0].split(':')[1]
    if sq in u:
        u[sq] = u[sq] + '|' + nm
    else:
        u[sq] = nm
        
print('uhgp90\tafdb')
gis = list(v.keys())
for gi in gis:
    seq = v[gi]
    print('%s\t%s' %(gi, u.get(seq, '')))
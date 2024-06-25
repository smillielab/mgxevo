import os, re, sys



# --------
# InterPro
# --------
ipro = {}

# InterPro blast hits
# --------------
h = ['gene', 'accesion', 'accesion_desc','interpro','interpo_desc','go','pathway','evalue']
print('\t'.join(h))

for line in open(sys.argv[1]):
    line = line.strip().split('\t')
    pid = line[0]
    pdb = line[3]
    pdbi = line[4]
    pdbDecs=line[5]
    evalue = line[8]
    iproAc = ''
    iproDesc = '' 
    iproGo = ''
    iproPathway=''
    if len(line) == 13:    
        iproAc = line[11]
        iproDesc = line[12]
    elif len(line) == 14:
        iproAc = line[11]
        iproDesc = line[12]
        iproGo = line[13]
        iproGo = ','.join(iproGo.split('|'))
    elif len(line) == 15:
        iproAc = line[11]
        iproDesc = line[12]
        iproGo = line[13]
        iproPathway = line[14]
        iproGo = ','.join(iproGo.split('|'))
        iproPathway = ':'.join(iproPathway.split(': '))
        iproPathway = ', '.join(iproPathway.split('|'))
    else:
        continue
    
    out = [pid, ':'.join([pdb, pdbi]), pdbDecs, iproAc, iproDesc, iproGo, iproPathway,evalue]
    print('\t'.join(out))
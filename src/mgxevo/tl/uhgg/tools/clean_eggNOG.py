import os, re, sys



# --------
# InterPro
# --------
eggnog = {}

# InterPro blast hits
# --------------
# h = ['query_name', 'seed_eggNOG_ortholog', 'seed_ortholog_evalue', 'seed_ortholog_score', 'best_tax_level', 'Preferred_name', 'GOs', 'EC', 'KEGG_ko', 'KEGG_Pathway', 'KEGG_Module', 'KEGG_Reaction', 'KEGG_rclass', 'BRITE', 'KEGG_TC', 'CAZy', 'BiGG_Reaction', 'tax_scope', 'eggnog_og' , 'best_og', 'cog', 'text_desc']
h = ['gene', 'name', 'GO', 'enzyme', 'kegg_ko', 'kegg_pathway', 'kegg_module', 'kegg_reaction', 'kegg_rclass', 'kegg_tc', 'brite', 'cazy', 'bigg','eggnog_og', 'cog', 'eggnog_desc', 'evalue']
print('\t'.join(h))

for line in open(sys.argv[1]):
    line = line.strip().split('\t')
    if len(line) < 20:
        continue
    if len(line)==20:
        line += ['', '']
    gid = line[0]
    evalue = line[2]
    gene = line[5]
    go = line[6]
    ec = line[7]
    ko = line[8]
    pathway = line[9]
    module = line[10]
    reaction = line[11]
    rclass = line[12]
    brite = line[13]
    tc = line[14]
    cz = line[15]
    bigg = line[16]
    eog = line[18]
    bog = line[19]
    if bog !='NA|NA|NA':
        eog = ','.join([eog, bog])
    cog = line[20]
    desc = line[21]
    
    out = [gid, gene, go, ec, ko, pathway, module, reaction, rclass, tc, brite, cz, bigg, eog, cog, desc, evalue]
    print('\t'.join(out))
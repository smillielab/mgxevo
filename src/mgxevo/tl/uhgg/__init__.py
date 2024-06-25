import os
import pkg_resources

script_path = pkg_resources.resource_filename('mgxevo', 'tl/uhgg/')
data_path = pkg_resources.resource_filename('mgxevo', 'tl/uhgg/data/')

uhgg_link = 'ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v1.0/uhgg_catalogue'
uhgg_all_genomes_link = 'ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v1.0/all_genomes'
uhgg_kraken_link = 'ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v1.0/uhgg_kraken2-db'
uhgp_90_link = 'ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v1.0/uhgp_catalogue/uhgp-90.tar.gz'

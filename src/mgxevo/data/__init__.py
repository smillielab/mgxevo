import os, sys
import pkg_resources as pkg
import urllib.request
import pandas as pd
import logging

path = os.path.dirname(os.path.realpath(__file__))
metadata_path = os.path.join(path, "dx.tsv")

# Metadata

if not os.path.exists(metadata_path):
    from .. import ut
    logging.info('Downloading metadata file...')
    url = 'https://www.dropbox.com/scl/fi/zb43pt9owvw06ej41nalp/TableS1_deidentified.tsv?rlkey=m1g83zuiab4hqv1uqv7n5hupz'  # replace with actual dropbox link
    ut.run_cmd(f'wget -O {metadata_path} {url} > /dev/null 2>&1')

    # Fix rownames
    pd.read_table(metadata_path).set_index('Sample ID').to_csv(metadata_path, sep='\t')
    logging.info('Done.')

try:
    metadata_path = pkg.resource_filename('mgxevo.data', 'dx.tsv')
    metadata = pd.read_table(metadata_path)
    adapter_path = pkg.resource_filename('mgxevo.data', 'adapters.mgx.fasta')
except ModuleNotFoundError:
    from pathlib import Path
    cwd = Path(__file__).resolve().parent
    metadata_path = f"{cwd}/dx.tsv"
    metadata = pd.read_table(metadata_path)
    adapter_path = f"{cwd}/adapters.mgx.fasta"
    

def download_uhgg_cluster():
    """Download UHGG cluster"""
    from .. import ut
    logging.info('Downloading UHGG cluster...')
    url = 'https://www.dropbox.com/scl/fi/qv60mye08o09ws9c0byki/markers.tar.gz?rlkey=4iaijw0kwvbflezfcv3gcnv5x'
    # make dir
    ut.run_cmd(f'wget -O {path}/markers.tar.gz {url} > /dev/null 2>&1')
    ut.run_cmd(f'tar -xzf {path}/markers.tar.gz -C {path}/ > /dev/null 2>&1')
    ut.run_cmd(f'rm {path}/markers.tar.gz > /dev/null 2>&1')
    
    assert os.path.exists(f"{path}/markers/dnaG-nr.fst"), "Download failed"
    logging.info('Done.')

def download_patric_metadata():
    """Download patric metadata"""
    from .. import ut
    logging.info('Downloading patric metadata...')
    
    os.makedirs(f"{path}/patric", exist_ok=True)
    url = 'ftp://ftp.bvbrc.org/RELEASE_NOTES/genome_metadata'
    ut.run_cmd(f'wget -O {path}/patric/patric.genome_metadata.tsv {url} > /dev/null 2>&1')
    assert os.path.exists(f"{path}/patric/patric.genome_metadata.tsv"), "Download failed"

    logging.info('Done.')
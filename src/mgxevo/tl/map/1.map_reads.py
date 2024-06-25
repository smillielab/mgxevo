import os
import sys
import argparse
import logging
import subprocess
import concurrent.futures
from concurrent.futures import ProcessPoolExecutor
from functools import partial
from tqdm import tqdm

from pathlib import Path
sys.path.append(str(Path(__file__).resolve().parent.parent.parent))

from data import metadata, adapter_path

def run_command(cmd, log_file, dry_run):
    if dry_run:
        logdir = os.path.dirname(log_file)
        with open(f'{logdir}/commands.txt', 'a') as f:
            f.write(cmd + '\n')
    else:
        with open(log_file, 'a') as f:
            try:
                subprocess.run(cmd, shell=True, check=True, stdout=f, stderr=subprocess.STDOUT)
            except subprocess.CalledProcessError as e:
                logging.error(f"Command '{cmd}' failed with error: {str(e)}")

def process_file(file):
    
    x = file
    z = x.split('/')[-1].replace('.fastq.gz', '')
    if args.trimmed:
        x = x.replace('.fastq.gz', '.trim.fastq.gz')
    else:
        if 'R*' in x:
            ir1 = x.replace('R*', 'R1')
            ir2 = x.replace('R*', 'R2')
            o1 = os.path.join(args.fastq_dir, f'{z}_R1.trim.fastq.gz')
            o2 = os.path.join(args.fastq_dir, f'{z}_R2.trim.fastq.gz')
            t1 = os.path.join(args.fastq_dir, f'{z}_R1.unmatch.fastq.gz')
            t2 = os.path.join(args.fastq_dir, f'{z}_R2.unmatch.fastq.gz')
            
            cmd = f'java -jar trimmomatic-0.39.jar PE {ir1} {ir2} {o1} {t1} {o2} {t2} ILLUMINACLIP:{adapter_path}:2:30:10:8:TRUE MAXINFO:80:0.5 MINLEN:50 AVGQUAL:20 && sleep 10 && du -sh {o1} > {o1}.ok; rm {t1} {t2}'
        else:
            i = x
            o = os.path.join(args.fastq_dir, f'{z}.trim.fastq.gz')
            t = os.path.join(args.fastq_dir, f'{z}.unmatch.fastq.gz')
            cmd = f'java -jar trimmomatic-0.39.jar SE {i} {o} ILLUMINACLIP:{adapter_path}:2:30:10:8:TRUE MAXINFO:80:0.5 MINLEN:50 AVGQUAL:20 && sleep 10 && du -sh {o} > {o}.ok; rm {t}' 
            
        run_command(cmd, os.path.join(args.logdir, f'{z}.log'), args.dry_run)
        
    # bowtie2
    o = os.path.join(args.bowtie_outdir, f'{z}.sam')
    if 'R*' in x:
        z = x.split('/')[-1].rsplit('_', 1)[0]
        o = os.path.join(args.bowtie_outdir, f'{z}.sam')
        ir1 = x.replace('R*', 'R1')
        ir2 = x.replace('R*', 'R2')
        cmd = f'bowtie2 -p 1 -x {args.index} --very-sensitive -a --no-unal -1 {ir1} -2 {ir2} -S {o}'
    else:
        i = x
        cmd = f'bowtie2 -p 1 -x {args.index} --very-sensitive -a --no-unal -U {i} -S {o}'
    run_command(cmd, os.path.join(args.logdir, f'{z}.log'), args.dry_run)

    # samtools sort
    i = o
    o = os.path.join(args.bowtie_outdir, f'{z}.sort.sam')
    cmd = f'samtools sort {i} -O sam > {o}; rm {i}'
    run_command(cmd, os.path.join(args.logdir, f'{z}.log'), args.dry_run)

    # samtools filter
    i = o
    o = os.path.join(args.bowtie_outdir, f'{z}.filter.sam')
    cmd = f'python {script_dir}/tools/filter_sam.py {i} 90 25 > {o}; rm {i}'
    run_command(cmd, os.path.join(args.logdir, f'{z}.log'), args.dry_run)

    # kpileup
    i = os.path.join(args.bowtie_outdir, f'{z}.filter.sam')
    o1 = os.path.join(args.kpileup_outdir, f'{z}.kp.out')
    o2 = os.path.join(args.kpileup_outdir, f'{z}.kp.index')
    cmd = f'perl {script_dir}/tools/kpileup.pl {z} {i} {args.gene_file} 10 0 1 > {o1}; python {script_dir}/tools/kp_to_index.py {o1} {o2}; gzip {o1}'
    run_command(cmd, os.path.join(args.logdir, f'{z}.log'), args.dry_run)


def process_kpileup(gene):

    # convert kpileup to numpy
    o = os.path.join(args.kpileup_outdir, f'{gene}.aln.pickle')
    cmd = f'python {script_dir}/tools/kp_to_np.py --kp_dir {args.kpileup_outdir} --samples {args.sample_file} --gene {gene} --gene_file {args.gene_file} --out {o}'
    run_command(cmd, os.path.join(args.logdir, f'{gene}.log'), args.dry_run)


def main(args):
    # parallelize MGX samples
    with ProcessPoolExecutor(max_workers=args.threads) as executor:
        futures = {executor.submit(process_file, file): file for file in args.files}
        logging.info(f"Processing {len(args.files)} files for bowtie2 and kpileup")
        for future in tqdm(concurrent.futures.as_completed(futures), total=len(futures), desc="Processing files"):
            try:
                future.result()
            except Exception as exc:
                print('%r generated an exception: %s' % (futures[future], exc))
            
    # parallelize genes
    with ProcessPoolExecutor(max_workers=args.threads) as executor:
        genes = []
        with open(args.gene_file) as file:
            for line in file:
                line = line.rstrip().split('\t')
                gene = line[1]
                genes.append(gene)

        logging.info(f"Processing {len(genes)} genes for kpileup conversion")
        futures = {executor.submit(process_kpileup, gene): gene for gene in genes}
        for future in tqdm(concurrent.futures.as_completed(futures), total=len(futures), desc="Processing files"):
            try:
                future.result()
            except Exception as exc:
                print('%r generated an exception: %s' % (futures[future], exc))
                

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--index', help='index name', required=True)
    parser.add_argument('--fastq_dir', help='fastq directory', required=True)
    parser.add_argument('--outdir', help='output directory', required=True)
    parser.add_argument('--sample_file', help='samples file', required=True)
    parser.add_argument('--gene_file', help='gene file', default="all-nr.gene.txt", required=True)
    parser.add_argument('--metadata', help='metadata file', default=None)
    parser.add_argument('--logdir', help='log directory', default='logs')
    parser.add_argument('--threads', help='number of threads', default=8, type=int)
    parser.add_argument('--dry_run', help='dry run', action='store_true')
    parser.add_argument('--trimmed', help='whether the fastq is trimmed', action='store_true')
    args = parser.parse_args()
    
    if args.metadata is not None:
        metadata = pd.read_table(args.metadata, index_col=0)
    else:
        args.metadata = metadata
    
    args.bowtie_outdir = os.path.join(args.outdir, 'bowtie2')
    args.kpileup_outdir = os.path.join(args.outdir, 'kpileup')
    
    os.makedirs(args.logdir, exist_ok=True)
    os.makedirs(args.outdir, exist_ok=True)
    os.makedirs(args.bowtie_outdir, exist_ok=True)
    os.makedirs(args.kpileup_outdir, exist_ok=True)
    
    with open(args.sample_file) as f:
        args.samples = [x.rstrip() for x in f]

    args.files = metadata.query('`Sample ID` in @args.samples')['File path'].apply(lambda x: x.split('/')[-1]).tolist()
    args.files = [os.path.join(args.fastq_dir, x) for x in args.files]

    logging.basicConfig(level=logging.INFO)
    script_dir = Path(__file__).parent.absolute()

    main(args)
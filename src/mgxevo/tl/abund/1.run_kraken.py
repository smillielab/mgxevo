import os
import sys
import re
import argparse
import glob
from pathlib import Path
import concurrent.futures
from concurrent.futures import ThreadPoolExecutor
from tqdm import tqdm

sys.path.append(str(Path(__file__).resolve().parent.parent.parent))
cwd = str(Path(__file__).resolve().parent)

import ut

def main():
    parser = argparse.ArgumentParser(description='Run kraken and bracken for a list of (paired) fastq files')
    parser.add_argument('--fastq_dir', type=str, help='Path to the compressed fastq.gz file', required=True)
    parser.add_argument('--kraken_db_path', type=str, help='Path to the kraken database', required=True)
    parser.add_argument('--output_dir', type=str, help='Path to the output directory', default='.')
    parser.add_argument('--debug', help='Debug mode', default=False, action='store_true')
    parser.add_argument('--threads', type=int, help='Number of threads', default=None)
    parser.add_argument('--verbose', help='Verbose mode', default=False, action='store_true')
    args = parser.parse_args()

    kraken = f"kraken2 --gzip-compressed --db {args.kraken_db_path}"
    bracken = f"bracken -d {args.kraken_db_path} -r 100 -l S -t 10"

    args.fastq_file = glob.glob(args.fastq_dir + '/*.fastq.gz')

    if args.threads is None:
        args.threads = os.cpu_count() - 2
    
    with ThreadPoolExecutor(max_workers=args.threads) as executor:
        futures = {executor.submit(process_file, ifn.rstrip(), kraken, bracken, args) for ifn in args.fastq_file if '_R2' not in ifn}
        for future in tqdm(concurrent.futures.as_completed(futures), total=len(futures)):
            future.result()

def process_file(ifn, kraken, bracken, args):
    ifn1 = os.path.realpath(ifn)
    ifn2 = re.sub('_R1', '_R2', ifn1)
    kout = os.path.join(args.output_dir, re.sub('.trim.*', '', os.path.basename(ifn1)) + '.kraken.out')
    kout = re.sub('_R[12]', '', kout)
    bout = re.sub('kraken', 'bracken', kout)

    if os.path.exists(bout):
        return

    cmd = run_kraken(ifn1, ifn2, kout, kraken)
    cmd = fix_kraken(cmd, kout)
    cmd = run_bracken(cmd, kout, bout, bracken)
    
    if args.debug:
        print(cmd)
    else:
        ut.run_cmd(cmd, silent=not args.verbose)
    

def run_kraken(ifn1, ifn2, kout, kraken):
    if (ifn1 != ifn2) and (os.path.exists(ifn1) and os.path.exists(ifn2)):
        return f"{kraken} --report {kout} {ifn1} {ifn2} --memory-mapping > /dev/null"
    else:
        return f"{kraken} --report {kout} {ifn1} --memory-mapping > /dev/null"

def fix_kraken(cmd, kout):
    return f"{cmd}; python {cwd}/tools/fix_kraken.py {kout} > {kout}.tmp; mv {kout}.tmp {kout}"

def run_bracken(cmd, kout, bout, bracken):
    return f"{cmd}; {bracken} -i {kout} -o {bout}"

if __name__ == "__main__":
    main()


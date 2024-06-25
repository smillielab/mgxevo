import sys
import argparse
import tempfile
from pathlib import Path

sys.path.append(str(Path(__file__).resolve()))
sys.path.append(str(Path(__file__).resolve().parent.parent.parent))

import ut

def main():
    parser = argparse.ArgumentParser(description='Run diamond blastx')
    parser.add_argument('-q', '--query', help='Query file', required=True)
    parser.add_argument('-o', '--output', help='Output file', required=True)
    parser.add_argument('-d', '--database', help='Database file', required=True)
    parser.add_argument('--id', default=90, type=int, help='Identity percentage')
    parser.add_argument('--evalue', default=1e-6, type=float, help='E-value')
    parser.add_argument('--threads', default=1, type=int, help='Number of threads')
    parser.add_argument('-t', '--tempdir', help='Temporary directory')
    parser.add_argument('--debug', action='store_true', help='Debug mode')

    args = parser.parse_args()
    
    if args.tempdir is None:
        args.tempdir = tempfile.mkdtemp()

    cmd = f"diamond blastx -d {args.database} -q {args.query} -o {args.output} --id {args.id} --evalue {args.evalue} -k 1 --max-hsps 1 --fast --iterate --compress 1 --unal 0 -t {args.tempdir} --threads {args.threads}"
    ut.run_cmd(cmd, silent=not args.debug)
    

if __name__ == "__main__":
    main()
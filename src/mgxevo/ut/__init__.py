import logging, os, sys
import gzip, re, time
import subprocess

from . import seq as seq_util

level = logging.INFO
fmt = '[%(levelname)s] %(asctime)s - %(message)s'
logging.basicConfig(level=level, format=fmt)
logging.getLogger('numexpr').setLevel(logging.WARNING)

def change_dir(dir, nb, mktmpl=True, add_syspath=True, force=False):

    dir = str(dir)    
    if os.getcwd() != dir:
        try:
            os.chdir(dir)
        except:
            logging.error('Could not change directory to %s/%s' %(dir, nb))
            raise
    try:
        if os.getcwd() == '%s/%s' %(dir, nb):
            logging.info('Already in %s/%s' %(dir, nb))
        else:
            os.chdir('./'+nb)
            logging.info('Changed path to %s/%s' %(dir, nb))
    except:
        if not os.path.exists(nb):
            # if last two directory has the same name, print warning
            if (dir.split('/')[-1] == nb) and not force:
                logging.warning('Directory may have already been created, exiting')
                return 
            
            os.mkdir(nb)
            os.chdir('./'+nb)
            logging.info('Created directory for this notebook')

        if mktmpl:
            [os.makedirs(f, exist_ok=True) for f in ['write','figure','tmp','data']]

    if add_syspath:
        sys.path.insert(0, dir)

def run_cmd(command, log_file=None, silent=False):
    if isinstance(command, list):
        command = ' '.join(command)
    try:
        if silent:
            subprocess.run(command, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
        elif log_file:
            with open(log_file, 'a') as outfile:
                subprocess.run(command, shell=True, check=True, stdout=outfile, stderr=subprocess.STDOUT)
        else:
            subprocess.run(command, shell=True, check=True)
    
    except subprocess.CalledProcessError as e:
        logging.error(f"Command '{command}' failed with error: {e}")

        
def create_symlinks(base_dir, source_files, destination_dir):
    import os
    import glob

    source_files = glob.glob(os.path.join(base_dir, source_files))
    if len(source_files) == 0:
        raise ValueError('No files found matching %s' %(source_files))
    os.makedirs(os.path.join(base_dir, destination_dir), exist_ok=True)

    for file in source_files:
        destination = os.path.join(base_dir, destination_dir, os.path.basename(file))
        if os.path.exists(destination):
            logging.warning('File %s already exists, skipping' %(destination))
        else:
            os.symlink(file, destination)

def open_gz(fn):
    if type(fn) == type(''):
        if '.gz' in fn:
            fn = gzip.open(fn, 'rt')
        else:
            fn = open(fn)
    else:
        # hack: read from stdin
        pass    
    return fn


def message(text, indent=2):
    # print message to stderr
    space = ' ' * indent
    text = re.sub('\n', '\n%s' %(space), text)
    sys.stderr.write('%s%s\n' %(space, text))


def error(text, indent=2):
    # print message to stderr and quit
    space = ' ' * indent
    text = re.sub('\n', '\n%s' %(space), text)
    sys.stderr.write('%s%s\n' %(space, text))
    quit()


def detect_fasta_format(fn):
    # automatically detect FASTA/FASTQ format
    i = 0
    fmt = 'unknown'
    for line in open_gz(fn):
        i += 1
        if i == 1:
            if line.startswith('>'):
                fmt = 'fasta'
        if i == 3:
            if line.startswith('+'):
                fmt = 'fastq'
            return fmt


def iter_fst(fn, split=False):
    # generator that iterates through [sid, seq] pairs in a fasta file
    sid = ''
    seq = ''
    i = 0
    for line in open_gz(fn):
        i += 1
        if i == 1:
            if not line.startswith('>'):
                quit('error: invalid FASTA file %s' %(fn))
        line = line.rstrip()
        if line.startswith('>'):
            if seq != '':
                yield [sid, seq]
            sid = line
            if split is True:
                sid = sid.split()[0]
            seq = ''
        else:
            seq += line
    yield [sid, seq]


def iter_fsq(fn):
    # generator that iterates through records in a fastq file
    record = []
    i = 0
    for line in open_gz(fn):
        i += 1
        if i == 3:
            if not line.startswith('+'):
                quit('error: invalid FASTQ file %s' %(fn))
        if i % 4 == 1:
            if len(record) > 0:
                yield record
            record = []
        record.append(line.rstrip())
    yield record


def iter_seq(fn):
    if detect_fasta_format(fn) == 'fasta':
        return iter_fst(fn)
    else:
        return iter_fsq(fn)


def read_fst(fn, reverse=False, split=False):
    # read fasta file as dictionary
    fst = {}
    for [sid, seq] in iter_fst(fn, split=split):
        try:
            sid = sid[1:].split()[0]
        except:
            print('error:\nfn = %s\nsid = %s\nseq = %s' %(fn, sid, seq))
            continue
        if reverse == False:
            fst[sid] = seq
        elif reverse == True:
            fst[seq] = sid
    return fst


def cycle(x):
    # an efficient way to cycle through a list (similar to itertools.cycle)
    while True:
        for xi in x:
            yield xi


class timer():
    # generator that measures elapsed time
    def __init__(self):
        self.t = [time.time()]
    def __iter__(self):
        return self
    def set(self):
        self.t.append(time.time())
        return self.t[-1] - self.t[-2]
    def mean(self):
        return np.mean(self.t)
    def median(self):
        return np.median(self.t)

def reverse_complement(x):
    return x[::-1].translate(x.maketrans('acgtnACGTN', 'tgcanTGCAN'))

def slugify(value, allow_unicode=False):
    import unicodedata
    """
    https://github.com/django/django/blob/master/django/utils/text.py
    Convert to ASCII if 'allow_unicode' is False. Convert spaces to hyphens.
    Remove characters that aren't alphanumerics, underscores, or hyphens.
    Convert to lowercase. Also strip leading and trailing whitespace.
    """
    value = str(value)
    if allow_unicode:
        value = unicodedata.normalize('NFKC', value)
    else:
        value = unicodedata.normalize('NFKD', value).encode('ascii', 'ignore').decode('ascii')
    value = re.sub(r'[^\w\s-]', '', value).strip()
    return re.sub(r'[-\s]+', '-', value)

def get_valid_filename(s):
    """
    https://github.com/django/django/blob/master/django/utils/text.py
    Return the given string converted to a string that can be used for a clean
    filename. Remove leading and trailing spaces; convert other spaces to
    underscores; and remove anything that is not an alphanumeric, dash,
    underscore, or dot.
    >>> get_valid_filename("john's portrait in 2004.jpg")
    'johns_portrait_in_2004.jpg'
    """
    s = str(s).strip().replace(' ', '_')
    return re.sub(r'(?u)[^-\w.]', '', s)

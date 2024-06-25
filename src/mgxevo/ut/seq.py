import gzip, re, sys, time

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
                sys.exit('error: invalid FASTA file')
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
                sys.exit('error: invalid FASTQ file')
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

def dict2fasta(dict, output_name):
    with open(output_name, 'w') as f:
        for key, value in dict.items():
            f.write('>' + key + '\n' + value + '\n')


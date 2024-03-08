import sys
from Bio.Seq import Seq

allowed_chars = 'ATCG'
window = 19*3

def translator(seq, table={}):
    while True:
        try:
            codon = next(seq) + next(seq) + next(seq)
            yield table[codon]
        except StopIteration:
            return
        except KeyError:
            aa, = str(Seq(codon).translate())
            table[codon] = aa
            yield aa

def translate(seq):
    return ''.join(translator(iter(seq)))

def fastq_fragment_counts(fnames, leading_seq=None, trailing_seq=None):
    counts = {}
    for fname in fnames:
        with open(fname) as f:
            catted_file = f.read().replace('\n', '').replace(' ', '')#[:10000] # TODO
        #print(catted_file.lower())
        for i, j in zip(range(len(catted_file)), range(window, len(catted_file) + 1)):
            fragment = catted_file[i:j]
            if leading_seq and catted_file[i-len(leading_seq):i].lower() != \
                    leading_seq.lower():
                continue
            if trailing_seq and catted_file[j:j+len(trailing_seq)].lower() != \
                    trailing_seq.lower():
                continue
            if all((c in allowed_chars for c in fragment)):
                translation = translate(fragment)
                counts[translation] = counts.get(translation, 0) + 1
    return counts

def count_dict_sum(*count_dicts):
    sum_counts = {}
    for counts in count_dicts:
        for key in counts:
            sum_counts[key] = sum_counts.get(key, 0) + counts[key]
    return sum_counts

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('must supply input file names')
        exit()
    fnames = sys.argv[1:]
    counts = fastq_fragment_counts(fnames)
    import pdb; pdb.set_trace()

import os
from itertools import product, zip_longest
from Bio.Seq import Seq

ecoli_codons = {
'A': 'gcg',
'R': 'cgc',
'N': 'aac',
'D': 'gat',
'C': 'tgc',
'E': 'gaa',
'Q': 'cag',
'G': 'ggc',
'H': 'cat',
'I': 'att',
'L': 'ctg',
'K': 'aaa',
'M': 'atg',
'F': 'ttt',
'P': 'ccg',
'S': 'agc',
'T': 'acc',
'W': 'tgg',
'Y': 'tat',
'V': 'gtg'
}

all_aminos = list(ecoli_codons.keys())

def make_ordinal(n):
    '''
    Convert an integer into its string ordinal representation

        make_ordinal(0)   => '0th'
        make_ordinal(3)   => '3rd'
        make_ordinal(122) => '122nd'
        make_ordinal(213) => '213th'
    '''
    n = int(n)
    suffix = ['th', 'st', 'nd', 'rd', 'th'][min(n % 10, 4)]
    if 11 <= (n % 100) <= 13:
        suffix = 'th'
    return str(n) + suffix


def unique(comment_list):
    seen = set()
    out = CommentedList()
    for item, comment in zip_longest(comment_list, comment_list.comments):
        if item not in seen:
            out.append(item)
            out.comments.append(comment)
        seen.add(item)
    return out

class CommentedList(list):
    def __init__(self, *args, comments=None):
        super().__init__(*args)
        self.comments = comments if comments else []

    def extend(self, new_items):
        # make sure we're filled up with None entries to this point
        self.comments = [comment for comment, _ in zip_longest(self.comments, self)]
        super().extend(new_items)
        try:
            self.comments.extend(new_items.comments)
        except AttributeError:
            pass

translation_mem = {}
def aa_variants(aa_seq_in, replace_map, source_dna_seq=None, annotations=None, padding=0):
    """
    E. coli-codon-optimized dna variants of an amino acid sequence.
    aa_seq_in: string, single-letter amino acid codes
    replace_map: dictionary mapping tuples to lists.
        - Tuples are (index in aa_seq_in, single-letter aa in aa_seq_in) for validation.
        - Lists are single-letter aa codes of replacement options.
    [source_dna_seq]: An optional dna string of arbitrary length.
        - If given, the first instance of aa_seq_in in its translation will be used to
          fill in codons for unchanged aas. All 3 frames are checked in sequence.
    [annotations]: An optional list (of the same length) of strings that will appear in
          comments.
    [padding]: Only if source_dna_seq provided, optional integer number of bases from
          original sequence to include before and after the variants

    returns: CommentedList (behaves as a normal list) of variants as [atcg]* strings
        - All replacement combinations are included, including the original aa sequence.
        - The CommentedList().comments attribute is a human-readable parallel list of
          strings describing what happened to make each variant
    """
    if source_dna_seq is None:
        orig_dna_seq = ''.join((ecoli_codons[aa] for aa in aa_seq_in))
        left_flank = right_flank = ''
    else:
        for offset in range(3):
            frame = source_dna_seq[offset:]
            frame = frame[:len(frame)//3*3] # trim to satisfy BioPython
            translation = translation_mem.get(frame, str(Seq(frame).translate()))
            translation_mem[frame] = translation
            try:
                pattern_index = translation.index(aa_seq_in)*3
            except ValueError:
                continue
            end_index = pattern_index+len(aa_seq_in)*3
            orig_dna_seq = frame[pattern_index:end_index]
            left_flank = frame[pattern_index-padding:pattern_index]
            right_flank = frame[end_index:end_index+padding]
            break
        else:
            raise ValueError('source DNA given, but aa sequence ' + aa_seq_in + ' does not appear in any frame of its translation.')
    length = len(aa_seq_in)
    assert len(orig_dna_seq)/3 == length
    for replace_idx, replace_aa in replace_map:
        # make sure identity replace is in replace map
        replace_map[(replace_idx, replace_aa)] = set(list(replace_map[(replace_idx, replace_aa)]) + [replace_aa])
    all_replaces = [[(idx_and_aa, replace_aa) for replace_aa in replace_map[idx_and_aa]] for idx_and_aa in replace_map]
    variants = []
    comments = []
    for replace_combo in product(*all_replaces):
        new_aa_seq = aa_seq_in
        replacement_comments = []
        replace_count = 0
        for (replace_idx, old_aa), new_aa in replace_combo:
            if new_aa_seq[replace_idx] != old_aa:
                raise ValueError # (index, aa) did not match supplied aa seq or there were duplicates
            new_aa_seq = new_aa_seq[:replace_idx] + new_aa + new_aa_seq[replace_idx+1:]
            if new_aa != old_aa:
                replacement_comments.append(new_aa + ' for ' + make_ordinal(replace_idx+1) + ' aa, ' + old_aa + (' (' + annotations[replace_idx] + ')' if annotations else ''))
                replace_count += 1
        seq_out = ''
        for old_aa, new_aa, dna_idx in zip(aa_seq_in, new_aa_seq, range(0, len(orig_dna_seq), 3)):
            if new_aa == old_aa:
                seq_out += str(orig_dna_seq[dna_idx:dna_idx+3]).lower()
            else:
                seq_out += ecoli_codons[new_aa]

        comment = 'fragment ' + aa_seq_in + '. replacements: ' + str(replace_count) + '. '
        if not replacement_comments:
            comment += 'outputting unchanged sequence: ' + new_aa_seq
        else:
            comment += ', '.join(replacement_comments) + '. final aa sequence: ' + new_aa_seq + '. ' + str(padding) + ' bases of padding on each side'
        if len(left_flank) != padding or len(right_flank) != padding:
            raise RuntimeError('Didn\'t correctly generate padding.')
        variants.append(left_flank + seq_out + right_flank) # include padding
        comments.append(comment)

    return unique(CommentedList(variants, comments=comments))

def single_variants(aa_seq_in, orig_dna_seq=None, annotations=None, padding=0):
    variants = CommentedList()
    for idx_and_aa in zip(range(len(aa_seq_in)), aa_seq_in):
        replacements = {idx_and_aa:all_aminos}
        variants.extend(aa_variants(aa_seq_in, replacements, orig_dna_seq, annotations=annotations, padding=padding))
    return unique(variants)

if __name__ == '__main__':
    try:
        with open(os.path.join('resources', 'm13.seq')) as f:
            m13 = f.read()
    except FileNotFoundError:
        m13 = None
    aa_fragment = 'NSPLMNNFRQYLPSLPQSV'
    neutral_replacements = {
        (0, 'N'): ['T', 'S'],
        (1, 'S'): ['A', 'T'],
        (2, 'P'): ['A', 'T'],
        (5, 'N'): ['S', 'D'],
        (10, 'Y'): ['H', 'F'],
        (13, 'S'): ['A', 'T']
    }
    variants = single_variants(aa_fragment, m13, padding=5)
    variants.extend(aa_variants(aa_fragment, neutral_replacements, m13, padding=5))
    variants = unique(variants)
    print('Generated', len(variants), 'variants')
    with open('variants.txt', 'w+') as f:
        for variant in variants:
            f.write(variant + '\n')


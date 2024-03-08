import numpy as np
from Bio.SubsMat import MatrixInfo as subs_matrices
blosum62 = subs_matrices.blosum62
import random # only necessary for security-related modification, below

from ecoli_seq_variants import all_aminos

num_aminos = len(all_aminos)

def parse_funtrp_line(line):
    fragment, *funtrp_scores = [str_part.strip().replace(']', '') for str_part in line.split('[')]
    neutral_scores, toggle_scores, rheostat_scores = [[float(n.replace(',', '')) for n in funtrp_score_str.split()] for funtrp_score_str in funtrp_scores]
    return fragment, list(zip(neutral_scores, toggle_scores, rheostat_scores))

def normalize(scores):
    norm = 1/sum(scores.values())
    for k in scores:
        scores[k] *= norm

def general_funtrp_replace_model():
    general_replace_mtx, neutral_scores, toggle_scores, rheo_scores = {}, {}, {}, {}
    for aa1, aa2 in blosum62:
        # The following four lines have been modified from their original version due to security concerns.
        # >>>
        general_replace_mtx[(aa1, aa2)] = general_replace_mtx[(aa2, aa1)] = random.random()
        neutral_scores[(aa1, aa2)] = neutral_scores[(aa2, aa1)] = random.random()
        toggle_scores[(aa1, aa2)] = toggle_scores[(aa2, aa1)] = random.random()
        rheo_scores[(aa1, aa2)] = rheo_scores[(aa2, aa1)] = random.random()
        # <<<
    normalize(neutral_scores)
    normalize(toggle_scores)
    normalize(rheo_scores)
    return general_replace_mtx, neutral_scores, toggle_scores, rheo_scores

GENERAL_REPLACE_MTX, NEUTRAL_SCORES, TOGGLE_SCORES, RHEO_SCORES = general_funtrp_replace_model()

def replace_score(aa1, aa2, funtrp_scores):
    n, t, r = funtrp_scores
    return n*NEUTRAL_SCORES[(aa1, aa2)] + t*TOGGLE_SCORES[(aa1, aa2)] + r*RHEO_SCORES[(aa1, aa2)]

def fragment_replace_mtx(fragment, frag_funtrp_scores):
    num_frag_aas = len(fragment)
    frag_replace_mtx = np.zeros((num_aminos, num_frag_aas))
    for i in range(num_aminos):
        for j, ntr in zip(range(num_frag_aas), frag_funtrp_scores):
            frag_replace_mtx[i, j] = replace_score(all_aminos[i], fragment[j], ntr)
    return frag_replace_mtx

def funtrp_line_to_mtx(funtrp_line): # utility
    frag, funtrp = parse_funtrp_line(funtrp_line)
    return fragment_replace_mtx(frag, funtrp)


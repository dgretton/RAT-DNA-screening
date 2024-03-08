import sys, os
import math
import numpy as np
from metro_hastings_variants import MetroHastingsVariants
from funtrp_blosum_predict import parse_funtrp_line

"from variant_tools import build_fitness_dict, funtrp_datas, mhv_for_fragment, mhv_likelihood"

# The provided data
from windows_list import windows_list

funtrp_datas = {fragment: funtrp_data for fragment, funtrp_data in [parse_funtrp_line(window) for window in windows_list]}

# generic edit distance function
def hamming_distance(s1, s2):
    return sum([1 for i in range(len(s1)) if s1[i] != s2[i]])

def mhv_for_fragment(fragment):
    funtrp_scores = funtrp_datas[fragment]
    return MetroHastingsVariants(fragment, funtrp_scores, result_queue_size=1)

def mhv_likelihood(variant, mhv):
    return mhv.score_frag(mhv.int_encode_variant(variant))

def build_fitness_dict(fragment, dataframe):
    """Constructs a dictionary for fast fitness lookups."""
    # Convert the DataFrame to a dictionary
    fitness_dict = dataframe.set_index('variant_aa_seq')['log_fold_enrichment'].to_dict()

    # Get the fitness for the original fragment
    orig_fitness = fitness_dict.get(fragment, None)
    if orig_fitness is None:
        raise ValueError(f"No fitness found for original fragment {fragment}")

    # Construct a dictionary with differences from original fragment's fitness
    for variant, fitness in fitness_dict.items():
        fitness_dict[variant] = fitness - orig_fitness

    return fitness_dict


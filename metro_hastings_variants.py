import numpy as np
from queue import Queue
import random
import threading
import json
from ecoli_seq_variants import all_aminos
from funtrp_blosum_predict import fragment_replace_mtx

num_aminos = len(all_aminos)

def int_encode_variant(variant):
    return tuple((all_aminos.index(c) for c in variant))

def hamming_dist(orig, variant):
    if len(orig) != len(variant):
        raise ValueError
    return sum(a != b for a, b in zip(orig, variant))

class MetroHastingsVariants:
    def __init__(self, fragment, funtrp_triples, result_queue_size, output_queue=None, event=None):
        self.fragment = fragment
        self.num_frag_aas = len(fragment)
        self.variant_scores = {}
        self.orig_coded_frag = self.int_encode_variant(fragment)
        self.replace_mtx = fragment_replace_mtx(fragment, funtrp_triples)
        self.score_frag(self.orig_coded_frag) # initialize known variant scores with original fragment
        self.log_replace_mtx = np.log(self.replace_mtx)
        self.output_queue = output_queue if output_queue is not None else Queue()
        self.event_initial_frags_added = event if event else threading.Event() # used in mp to check if initial frags already added to the shared queue.
        self.set_to_check = set() # used to check for repeat values.
        self.initial_variants = None
        self.hist = np.zeros_like(self.replace_mtx)
        self.result_queue_size = result_queue_size
        # >90% of mutations lead to an increase in edit distance from the original. to keep the search grounded near the original sequence, we reset the current proposal to some sequence near the original with some probability return_thres_prob.
        self.return_thres_prob = .1

    def int_encode_variant(self, variant, f=int_encode_variant):
        return f(variant) # local variable optimization

    def score_frag(self, coded_frag, enumerate=enumerate, sum=sum): # local var opt
        if coded_frag in self.variant_scores:
            return self.variant_scores[coded_frag]
        score = sum((self.replace_mtx[i,j] for j, i in enumerate(coded_frag)))
        self.variant_scores[coded_frag] = score
        return score

    def decode_variant(self, coded_var):
        return ''.join((all_aminos[j] for j in coded_var))

    # Metropolis-Hastings adapted from https://towardsdatascience.com/from-scratch-bayesian-inference-markov-chain-monte-carlo-and-metropolis-hastings-in-python-ef21a29e25a
    # The transition model defines how to move from the current fragment to a new fragment
    def transition_model(self, x, # randomly mutate one amino acid.
            idx_choices=[], # make fewer objects by creating local variables
            amino_choices=np.array(range(num_aminos)),
            np=np, range=range): # local var optim.
        if not idx_choices:
            idx_choices.append(np.array(range(self.num_frag_aas)))
        idx = np.random.choice(idx_choices[0])
        newx = list(x[:idx])
        newx.append(np.random.choice(amino_choices))
        newx.extend(x[idx + 1:])
        return tuple(newx)

    #Computes the likelihood of the data given a fragment (new or current) according to equation (2)
    def likelihood_computer(self, x,
            sum=sum, enumerate=enumerate): # local var opt
        return sum((self.log_replace_mtx[i,j] for j, i in enumerate(x)))

    # Defines whether to accept or reject the new sample.
    # Returns a tuple (accept the new sample?, did we accept based on a random decision because x_new is a less likely sample?)
    def acceptance_rule(self, x, x_new, np=np): # local var opt
        if x_new > x:
            return True, False
        else:
            accept = np.random.uniform(0,1)
            # Since we did a log likelihood, we need to exponentiate in order to compare to the random number
            # less likely x_new are less likely to be accepted
            return accept < np.exp(x_new - x), True

    def end_condition(self):
        # if the size of the output data structure is greater/equal to the number of entries specified, returns true. else returns false
        return self.get_size() >= self.result_queue_size

    def get_size(self):
        # returns size of output data structure.
        return self.output_queue.qsize()

    def get_output_queue(self):
        # returns the output data structure.
        return self.output_queue

    def is_empty(self):
        # returns true is output queue is empty. else returns false.
        return self.get_size()==0

    def add_to_output_queue(self, to_add, meta=None):
        # adds to_add to output data structure, and to last known output set.
        if to_add in self.set_to_check:
           return
        json_str = json.dumps({'sequence':to_add, 'sequence_meta':meta})
        self.output_queue.put(json_str)
        self.set_to_check.add(to_add)

    def run(self, min_interval=0, enumerate=enumerate): # local var opt
        # min_interval is the minimum number of steps between adding a new fragment to the output queue, functioning as a recurring burn-in period.
        # likelihood_computer(x): returns the likelihood of this sample
        # transition_model(x): a function that draws a sample from a symmetric distribution and returns it
        # acceptance_rule(x,x_new): decides whether to accept or reject the new sample

        wt_likelihood = self.likelihood_computer(self.int_encode_variant(self.fragment))
        self.add_to_output_queue(self.fragment, meta={'provenance':'wild_type',
                                                      'likelihood':wt_likelihood})
        def all_2_variants_and_liks():
            for i in range(len(self.orig_coded_frag)):
                for j in range(i + 1, len(self.orig_coded_frag)):
                    for i_replacement in range(len(all_aminos)):
                        for j_replacement in range(len(all_aminos)):
                            twovar = list(self.orig_coded_frag[:i])
                            twovar.append(i_replacement)
                            twovar.extend(self.orig_coded_frag[i + 1:j])
                            twovar.append(j_replacement)
                            twovar.extend(self.orig_coded_frag[j + 1:])
                            twovar = tuple(twovar)
                            yield (twovar, self.likelihood_computer(twovar))
        all_2_variants_map = {v:l for v, l in all_2_variants_and_liks()}

        # start with top-45%-likelihood 2-variants. tuple for efficient random choice.
        num_init_variants = int(min(len(all_2_variants_map)*.45, self.result_queue_size))

        self.initial_variants = tuple(list(sorted(all_2_variants_map,
                key=lambda v:all_2_variants_map[v]))[-num_init_variants:])

        if not self.event_initial_frags_added.is_set() and self.get_size()==0:
            self.event_initial_frags_added.set()
            for init_var in self.initial_variants:
                self.add_to_output_queue(self.decode_variant(init_var),
                        meta={'provenance':'initial_2_variant',
                              'likelihood':all_2_variants_map[init_var]})

        for v in self.initial_variants:
            for j, i in enumerate(v):
                self.hist[i, j] += 1

        if self.result_queue_size == num_init_variants:
            return [], []

        x = self.orig_coded_frag
        steps_since_added = 0

        while not self.end_condition():
            x_new =  self.transition_model(x)
            x_lik = self.likelihood_computer(x)
            x_new_lik = self.likelihood_computer(x_new)
            steps_since_added += 1
            accept_x, accepted_less_likely = self.acceptance_rule(x_lik, x_new_lik)
            if accept_x:
                if accepted_less_likely and np.random.random() < self.return_thres_prob:
                    # with each acceptance, reset the next sample to be the
                    # within 2 mutations of the original fragment with some probability, to make sure
                    # lots of proposals happen near the original without
                    # limiting the search when acceptances are rare.
                    next_x = random.choice(self.initial_variants)
                    next_lik = self.likelihood_computer(next_x)
                else:
                    next_x = x_new
                    next_lik = x_new_lik
                if steps_since_added < min_interval:
                    continue
                decoded_var = self.decode_variant(next_x)
                if decoded_var not in self.set_to_check:
                    steps_since_added = 0
                    for j, i in enumerate(next_x):
                        self.hist[i, j] += 1
                self.add_to_output_queue(decoded_var, meta={'provenance':'mh_sample',
                                                            'likelihood':next_lik})
                x = next_x

        return self.output_queue

    def result_set(self):
        # with parallelized MHV queue update, not needed.
        raise NotImplementedError

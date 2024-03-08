from metro_hastings_variants import MetroHastingsVariants
import numpy as np
import time
import json
import jsonf
from datetime import datetime
from funtrp_blosum_predict import parse_funtrp_line
import multiprocessing as mp
from qsizeQueue import PortableQueue

class ParallelizedMetroHastingsVariants:
    def __init__(self, fragment, ntr_triples, num_entries, meta=None, sequence_meta=None):
        # creates a empty (list) collection which will contain mhv instances
        self.fragment_shared = fragment
        self.ntr_shared = ntr_triples
        self.num_entries = num_entries
        q = PortableQueue()
        self.shared_queue = q
        self.process_lock = mp.RLock()
        self.final_result_set = set()
        self.event_initial_frags_added = mp.Event()
        self.meta = meta if meta else {}
        self.default_sequence_meta = sequence_meta if sequence_meta else {}
        self.sequence_meta = {}
        # create each new mhv instance with the same shared queue
        self.metro_collection = []
        for _ in range(mp.cpu_count()):
            mhv = MetroHastingsVariants(self.fragment_shared,
                    self.ntr_shared,
                    self.num_entries,
                    self.shared_queue,
                    self.event_initial_frags_added)
            self.metro_collection.append(mhv)

    def get_size_collection(self):
        # returns the size of the collection
        return len(self.metro_collection)

    def run_parallel(self):
        procs = [] # list with all processes so they can be started and terminated.
        sequence_meta_str = json.dumps(self.default_sequence_meta)
        all_seq_meta_strs = {}
        overwrite_warned = False
        self.meta['pmhv_start_time'] = str(datetime.now())

        try:
            for mhv in self.metro_collection:
                # run each mhv in a new process
                proc_mhv_run = mp.Process(target=mhv.run)
                proc_mhv_run.start()
                procs.append(proc_mhv_run)
            fname = f'resources/ordered_mhv_db_{self.fragment_shared}.txt'
            with open(fname, 'w+') as f:
                while not len(self.final_result_set) == self.num_entries:
                    item = json.loads(self.shared_queue.get())

                    # add to the output set
                    variant_seq = item['sequence']
                    if variant_seq not in self.final_result_set:
                        f.write(f"{variant_seq}\t{item['sequence_meta']['likelihood']}\n")
                    self.final_result_set.add(variant_seq)

                    # following section adds metadata to a list for each fragment; deduplicates identical
                    # metadata; keeps entries for identical sequences that may have multiple sources
                    sequence_meta = item['sequence_meta']
                    variant_default_meta = json.loads(sequence_meta_str) # deep copy
                    # exercise caution regarding key conflicts.
                    # in case of a key conflict, the following section will overwrite top-level
                    # keys regardless of how much nested metadata there may have been under those keys.
                    shared_keys = set(variant_default_meta).intersection(set(sequence_meta))
                    if not overwrite_warned and shared_keys:
                        print(f'WARNING: ParallelizedMetroHastingsVariants is overwriting the following conflicting sequence metadata keys: {", ".join(shared_keys)}. Following similar warnings squelched.')
                        overwrite_warned = True
                    variant_meta = variant_default_meta
                    variant_meta.update(sequence_meta)
                    variant_meta_str = json.dumps(variant_meta, sort_keys=True)
                    prior_meta_strs = all_seq_meta_strs.get(variant_seq, [])
                    if variant_meta_str not in prior_meta_strs:
                        prior_meta_strs.append(variant_meta_str)
                    all_seq_meta_strs[variant_seq] = prior_meta_strs
            print(f'ParallelizedMetroHastingsVariants wrote {len(self.final_result_set)} entries to {fname}.')

        finally:
            for proc in procs:
                proc.terminate()

        for variant_seq, seq_meta_strs in all_seq_meta_strs.items():
            if variant_seq not in self.final_result_set: # in case of an extra few queued before term
                continue
            self.sequence_meta[variant_seq] = [json.loads(mstr) for mstr in seq_meta_strs]

        self.meta['pmhv_end_time'] = str(datetime.now())

    def get_parallel_results(self):
        if not self.final_result_set:
            self.run_parallel()
        if len(self.final_result_set) != self.num_entries:
            raise RuntimeError('PMHV result set is not empty, but has wrong number of entries')
        return self.final_result_set

    def get_meta_results(self):
        if not self.final_result_set:
            raise RuntimeError('can only get metadata results from a '
                               'ParallelizedMetroHastingsVariants that has already run using '
                               'run_parallel() or get_parallel_results().')
        return self.sequence_meta

    def get_results_incl_meta(self):
        sequence_meta = self.get_meta_results() # to raise exception if needed
        fragset = jsonf.json_initialize_meta_dict()
        jsonf.json_merge_into(fragset, sequences=self.final_result_set,
                sequence_meta=sequence_meta, meta=self.meta)
        return fragset

    def write_meta(self, fname):
        try:
            output_fragset = jsonf.jsonf(fname)
        except FileNotFoundError:
            output_fragset = jsonf.json_initialize_meta_dict()
        mhv_fragset = self.get_results_incl_meta()
        jsonf.json_merge_fragset_into(mhv_fragset, output_fragset)
        jsonf.jsonw(output_fragset, fname)

def time_exec(target):
    start = time.time()
    target()
    return time.time() - start


if __name__ == '__main__':
    windows_list = [
    "MAVYFVTGKLGSGKTLVSV [0.14, 0.18, 0.01, 0.05, 0.05, 0.0, 0.01, 0.04, 0.13, 0.03, 0.08, 0.2, 0.06, 0.05, 0.02, 0.06, 0.1, 0.1, 0.15] [0.35, 0.21, 0.78, 0.56, 0.6, 0.85, 0.92, 0.9, 0.69, 0.34, 0.51, 0.46, 0.58, 0.53, 0.69, 0.41, 0.12, 0.47, 0.46] [0.51, 0.61, 0.21, 0.39, 0.35, 0.15, 0.07, 0.06, 0.18, 0.63, 0.41, 0.34, 0.36, 0.42, 0.29, 0.53, 0.78, 0.43, 0.39]",
    "YSYLTPYLSHGRYFKPLNL [0.08, 0.09, 0.05, 0.03, 0.04, 0.04, 0.07, 0.02, 0.11, 0.1, 0.01, 0.03, 0.09, 0.11, 0.06, 0.05, 0.04, 0.14, 0.07] [0.27, 0.4, 0.33, 0.27, 0.28, 0.29, 0.18, 0.06, 0.17, 0.14, 0.37, 0.21, 0.28, 0.25, 0.19, 0.15, 0.07, 0.18, 0.06] [0.65, 0.51, 0.62, 0.7, 0.68, 0.67, 0.75, 0.92, 0.72, 0.76, 0.62, 0.76, 0.63, 0.64, 0.75, 0.8, 0.89, 0.68, 0.87]",
    "VEIKASPAKVLQGHNVFGT [0.0, 0.58, 0.03, 0.11, 0.01, 0.02, 0.0, 0.0, 0.1, 0.03, 0.04, 0.13, 0.0, 0.03, 0.01, 0.0, 0.02, 0.0, 0.01] [0.86, 0.09, 0.4, 0.7, 0.87, 0.94, 0.95, 0.9, 0.73, 0.37, 0.31, 0.74, 0.98, 0.92, 0.93, 0.89, 0.92, 0.95, 0.95] [0.14, 0.33, 0.57, 0.19, 0.12, 0.04, 0.05, 0.1, 0.17, 0.6, 0.65, 0.13, 0.02, 0.05, 0.06, 0.11, 0.06, 0.05, 0.04]",
    "NFYPCVEIKASPAKVLQGH [0.92, 0.44, 0.83, 0.07, 0.13, 0.0, 0.58, 0.03, 0.11, 0.01, 0.02, 0.0, 0.0, 0.1, 0.03, 0.04, 0.13, 0.0, 0.03] [0.02, 0.16, 0.04, 0.86, 0.3, 0.86, 0.09, 0.4, 0.7, 0.87, 0.94, 0.95, 0.9, 0.73, 0.37, 0.31, 0.74, 0.98, 0.92] [0.06, 0.4, 0.13, 0.07, 0.57, 0.14, 0.33, 0.57, 0.19, 0.12, 0.04, 0.05, 0.1, 0.17, 0.6, 0.65, 0.13, 0.02, 0.05]",
    "PQSVECRPFVFGAGKPYEF [0.0, 0.3, 0.01, 0.0, 0.13, 0.03, 0.02, 0.03, 0.03, 0.0, 0.04, 0.02, 0.23, 0.0, 0.8, 0.76, 0.06, 0.14, 0.0] [0.89, 0.16, 0.62, 0.93, 0.55, 0.88, 0.67, 0.43, 0.95, 0.78, 0.57, 0.47, 0.16, 0.83, 0.07, 0.03, 0.34, 0.42, 0.89] [0.11, 0.54, 0.37, 0.07, 0.32, 0.09, 0.31, 0.54, 0.02, 0.22, 0.39, 0.51, 0.61, 0.17, 0.13, 0.21, 0.6, 0.44, 0.11]",
    "FRGVFAFLLYVATFMYVFS [0.08, 0.02, 0.09, 0.14, 0.22, 0.08, 0.19, 0.06, 0.03, 0.03, 0.07, 0.08, 0.03, 0.12, 0.06, 0.02, 0.04, 0.14, 0.1] [0.52, 0.27, 0.24, 0.52, 0.43, 0.64, 0.48, 0.5, 0.54, 0.64, 0.25, 0.16, 0.71, 0.46, 0.69, 0.69, 0.52, 0.39, 0.38] [0.4, 0.71, 0.67, 0.34, 0.35, 0.28, 0.33, 0.44, 0.43, 0.33, 0.68, 0.76, 0.26, 0.42, 0.25, 0.29, 0.44, 0.47, 0.52]",
    "YANYEGCLWNATGVVVCTG [0.13, 0.13, 0.2, 0.22, 0.38, 0.01, 0.09, 0.11, 0.14, 0.19, 0.16, 0.13, 0.09, 0.13, 0.14, 0.11, 0.12, 0.18, 0.37] [0.2, 0.19, 0.32, 0.14, 0.08, 0.69, 0.21, 0.04, 0.33, 0.21, 0.31, 0.17, 0.17, 0.23, 0.07, 0.08, 0.18, 0.08, 0.1] [0.67, 0.68, 0.48, 0.64, 0.54, 0.3, 0.7, 0.85, 0.53, 0.6, 0.53, 0.7, 0.74, 0.64, 0.79, 0.81, 0.7, 0.74, 0.53]",
    "IATTVNLRDGQTLLLGGLT [0.15, 0.48, 0.02, 0.28, 0.01, 0.3, 0.04, 0.68, 0.17, 0.05, 0.45, 0.03, 0.09, 0.01, 0.03, 0.03, 0.01, 0.02, 0.23] [0.12, 0.08, 0.47, 0.1, 0.31, 0.26, 0.22, 0.15, 0.42, 0.21, 0.22, 0.48, 0.1, 0.31, 0.26, 0.61, 0.64, 0.41, 0.13] [0.73, 0.44, 0.51, 0.62, 0.68, 0.44, 0.74, 0.17, 0.41, 0.74, 0.33, 0.49, 0.81, 0.68, 0.71, 0.36, 0.35, 0.57, 0.64]",
    "LLDVNATTISRIDAAFSAR [0.16, 0.11, 0.48, 0.31, 0.87, 0.54, 0.02, 0.38, 0.0, 0.08, 0.02, 0.0, 0.0, 0.01, 0.0, 0.06, 0.03, 0.03, 0.38] [0.13, 0.23, 0.26, 0.12, 0.06, 0.13, 0.92, 0.13, 0.55, 0.22, 0.68, 0.55, 0.91, 0.52, 0.87, 0.26, 0.9, 0.31, 0.34] [0.71, 0.66, 0.26, 0.57, 0.07, 0.33, 0.06, 0.49, 0.45, 0.7, 0.3, 0.45, 0.09, 0.47, 0.13, 0.68, 0.07, 0.66, 0.28]"
    ]

    print(time.time())
    for window in windows_list:
        fragment, ntr_triples = parse_funtrp_line(window)
        for num_entries in np.array([10])**np.linspace(1,4,4):
            mhv_wrapped = ParallelizedMetroHastingsVariants(fragment, ntr_triples, num_entries)
            mhv_wrapped.run_parallel()
            print('Length result list of', fragment, len(mhv_wrapped.get_parallel_results()))
            print('Variant database results (max 40 shown):', list(mhv_wrapped.get_parallel_results())[:40])
        print(time.time())


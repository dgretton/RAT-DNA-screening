from Bio import SeqIO
from Bio.Seq import Seq
import os
import sys
import json
import numpy as np
from pandas import DataFrame

from jsonf import jsonf
from ecoli_seq_variants import single_variants
from funtrp_blosum_predict import funtrp_line_to_mtx
from metro_hastings_variants import int_encode_variant
from fastq_counting import fastq_fragment_counts

ngs_dir = os.path.abspath('resources/NGS')
sys.path.append(ngs_dir)
from m13_ngs_data_frames import dataframes
from m13_ngs_data_frames import frag_12 as brian_fragment_20201210
import csv

#this one is treated separately because the initial version only had this data available
brian_data_20201210 = dataframes[brian_fragment_20201210]

gen_f_dir = os.path.join(ngs_dir, 'generated_files')
file_names = ["ESP7xLib14-postselection-phage.fastq"]
fragment = sys.argv[sys.argv.index('--fragment')+1] if '--fragment' in sys.argv else brian_fragment_20201210
making_figs = '--figs' in sys.argv
if not making_figs: print('not outputting figures.')
showing_figs = '--show' in sys.argv
if making_figs and not showing_figs: print('not showing figures, just saving.')

## Read the fastq files. Count the total number of sequences
all_data = {}
num_entries = 0
def load_all_data():
    global all_data, num_entries
    for f in file_names:
        all_data[f] = []
        with open(os.path.join(data_dir, f)) as f_handle:
            for record in SeqIO.parse(f_handle, "fastq"):
                num_entries += 1
                all_data[f].append(record.seq)
        print("Total reads", f, num_entries)

def pad_codons(seq):
    return seq + (-len(seq) % 3)*'N'

def translate(seq, mem={}):
    seq = seq.lower()
    if seq in mem:
       return mem[seq]
    padded = pad_codons(seq)
    try:
        translation = str(padded.translate())
    except (AttributeError, TypeError):
        padded = Seq(padded)
        translation = str(padded.translate())
    mem[seq] = translation
    return translation

def count_stats(seq):
    count_exact = 0
    count_trans = 0
    seq = seq.lower()
    trans = translate(seq)
    for fname in file_names:
        for single_read in all_data[fname]:
            rev_comp = single_read.reverse_complement()
            #print(seq, 'in', single_read, 'or', rev_comp, '?')
            if seq in str(single_read).lower() or seq in str(rev_comp).lower():
                count_exact += 1
            for frame in range(3):
                read_trans = translate(single_read[frame:])
                read_comp_trans = translate(rev_comp[frame:])
                if trans in read_trans or trans in read_comp_trans:
                    count_trans += 1
    assert count_trans >= count_exact
    return count_exact, count_trans

wt_sequence = 'tacgctaactatgagggctgtctgtggaatgctacaggcgttgtagtttgtactggt'.upper()
wt_trans = translate(wt_sequence)

single_vars = None
def load_single_vars():
    global single_vars
    single_vars = single_variants(wt_trans, orig_dna_seq=wt_sequence)

all_vars = all_transs = original_fragments = None
def load_all_vars_and_transs():
    global all_vars, all_transs, original_fragments
    with open(os.path.join('resources', 'extracted_submitted_aa_variants.json')) as f:
        all_vars = set(json.loads(f.read())) - {'', '\n', '\r\n'}
    all_transs = set((translate(var) for var in all_vars))
    original_fragments = jsonf(os.path.join('resources', 'original_fragment_map.json'))
    all_transs = set((v for v in all_transs if v == fragment or original_fragments.get(v, None) == fragment))
    print('number of transs:', len(all_transs))

def windows(seq):
    for i in range(len(seq)):
        window = seq[i:i+57]
        if len(window) != 57:
            return
        yield window.lower()

# count instances of this sequence
trans_counts = {}
def count_translations():
    count_wt, count_wt_trans = count_stats(wt_sequence)
    print("WT exact sequence count", count_wt)
    print("WT synonymous sequence count", count_wt_trans)
    for variant in all_vars:#[:40]:
        count_vari, count_vari_trans = count_stats(variant)
        trans_counts[translate(variant)] = count_vari_trans
        if count_vari_trans > 0:
            print(variant)
            print(translate(variant))
            print("Variant exact sequence count", count_vari)
            print("Variant synonymous sequence count", count_vari_trans)

count_is_wt = count_is_exact = count_is_aa_match = count_is_other = 0
covered_vars = set()
covered_transs = set()
def add_single_read(single_read):
    global count_is_wt, count_is_exact, count_is_aa_match, count_is_other
    rev_comp = single_read.reverse_complement()
    #print(seq, 'in', single_read, 'or', rev_comp, '?')
    seq = wt_sequence
    seq = seq.lower()
    if seq in str(single_read).lower() or seq in str(rev_comp).lower():
        print('found wild type')
        count_is_wt += 1
        return
    for window in windows(str(single_read)):
        if window in all_vars:
            covered_vars.add(window)
            covered_transs.add(window)
            count_is_exact += 1
            return
        var_translation = translate(window)
        if var_translation in all_transs:
            covered_transs.add(window)
            count_is_aa_match += 1
            return
    for window in windows(str(rev_comp)):
        if window in all_vars:
            covered_vars.add(window)
            covered_transs.add(window)
            count_is_exact += 1
            return
        var_translation = translate(window)
        if var_translation in all_transs:
            covered_transs.add(window)
            count_is_aa_match += 1
            return
    count_is_other += 1

baked_funtrp_lines = { # idea is to prevent accidentally redefining the fragment without using the correct funtrp scores EDIT Bless your tiny soul past dana. love, future dana
    'YANYEGCLWNATGVVVCTG':'YANYEGCLWNATGVVVCTG [0.13, 0.13, 0.2, 0.22, 0.38, 0.01, 0.09, 0.11, 0.14, 0.19, 0.16, 0.13, 0.09, 0.13, 0.14, 0.11, 0.12, 0.18, 0.37] [0.2, 0.19, 0.32, 0.14, 0.08, 0.69, 0.21, 0.04, 0.33, 0.21, 0.31, 0.17, 0.17, 0.23, 0.07, 0.08, 0.18, 0.08, 0.1] [0.67, 0.68, 0.48, 0.64, 0.54, 0.3, 0.7, 0.85, 0.53, 0.6, 0.53, 0.7, 0.74, 0.64, 0.79, 0.81, 0.7, 0.74, 0.53]',
    'PQSVECRPFVFGAGKPYEF':'PQSVECRPFVFGAGKPYEF [0, 0.3, 0.01, 0, 0.13, 0.03, 0.02, 0.03, 0.03, 0, 0.04, 0.02, 0.23, 0, 0.8, 0.76, 0.06, 0.14, 0] [0.89, 0.16, 0.62, 0.93, 0.55, 0.88, 0.67, 0.43, 0.95, 0.78, 0.57, 0.47, 0.16, 0.83, 0.07, 0.03, 0.34, 0.42, 0.89] [0.11, 0.54, 0.37, 0.07, 0.32, 0.09, 0.31, 0.54, 0.02, 0.22, 0.39, 0.51, 0.61, 0.17, 0.13, 0.21, 0.6, 0.44, 0.11]',
    'TKGDVENFSSLKKDVVIRV':'TKGDVENFSSLKKDVVIRV [0.49, 0.82, 0.17, 0.47, 0.77, 0.97, 0.86, 0.42, 0.44, 0.86, 0.43, 0.92, 0.69, 0.88, 0.57, 0.54, 0.57, 0.34, 0.59] [0.1, 0.11, 0.29, 0.23, 0.01, 0.0, 0.07, 0.11, 0.1, 0.02, 0.07, 0.03, 0.21, 0.09, 0.09, 0.06, 0.05, 0.12, 0.13] [0.41, 0.07, 0.54, 0.3, 0.22, 0.03, 0.07, 0.47, 0.46, 0.12, 0.5, 0.05, 0.1, 0.03, 0.34, 0.4, 0.38, 0.54, 0.28]',
    'MAVYFVTGKLGSGKTLVSV':'MAVYFVTGKLGSGKTLVSV [0.14, 0.18, 0.01, 0.05, 0.05, 0.0, 0.01, 0.04, 0.13, 0.03, 0.08, 0.2, 0.06, 0.05, 0.02, 0.06, 0.1, 0.1, 0.15] [0.35, 0.21, 0.78, 0.56, 0.6, 0.85, 0.92, 0.9, 0.69, 0.34, 0.51, 0.46, 0.58, 0.53, 0.69, 0.41, 0.12, 0.47, 0.46] [0.51, 0.61, 0.21, 0.39, 0.35, 0.15, 0.07, 0.06, 0.18, 0.63, 0.41, 0.34, 0.36, 0.42, 0.29, 0.53, 0.78, 0.43, 0.39]',
    'YSYLTPYLSHGRYFKPLNL':'YSYLTPYLSHGRYFKPLNL [0.08, 0.09, 0.05, 0.03, 0.04, 0.04, 0.07, 0.02, 0.11, 0.1, 0.01, 0.03, 0.09, 0.11, 0.06, 0.05, 0.04, 0.14, 0.07] [0.27, 0.4, 0.33, 0.27, 0.28, 0.29, 0.18, 0.06, 0.17, 0.14, 0.37, 0.21, 0.28, 0.25, 0.19, 0.15, 0.07, 0.18, 0.06] [0.65, 0.51, 0.62, 0.7, 0.68, 0.67, 0.75, 0.92, 0.72, 0.76, 0.62, 0.76, 0.63, 0.64, 0.75, 0.8, 0.89, 0.68, 0.87]',
    'VEIKASPAKVLQGHNVFGT':'VEIKASPAKVLQGHNVFGT [0.0, 0.58, 0.03, 0.11, 0.01, 0.02, 0.0, 0.0, 0.1, 0.03, 0.04, 0.13, 0.0, 0.03, 0.01, 0.0, 0.02, 0.0, 0.01] [0.86, 0.09, 0.4, 0.7, 0.87, 0.94, 0.95, 0.9, 0.73, 0.37, 0.31, 0.74, 0.98, 0.92, 0.93, 0.89, 0.92, 0.95, 0.95] [0.14, 0.33, 0.57, 0.19, 0.12, 0.04, 0.05, 0.1, 0.17, 0.6, 0.65, 0.13, 0.02, 0.05, 0.06, 0.11, 0.06, 0.05, 0.04]',
    'NFYPCVEIKASPAKVLQGH':'NFYPCVEIKASPAKVLQGH [0.92, 0.44, 0.83, 0.07, 0.13, 0.0, 0.58, 0.03, 0.11, 0.01, 0.02, 0.0, 0.0, 0.1, 0.03, 0.04, 0.13, 0.0, 0.03] [0.02, 0.16, 0.04, 0.86, 0.3, 0.86, 0.09, 0.4, 0.7, 0.87, 0.94, 0.95, 0.9, 0.73, 0.37, 0.31, 0.74, 0.98, 0.92] [0.06, 0.4, 0.13, 0.07, 0.57, 0.14, 0.33, 0.57, 0.19, 0.12, 0.04, 0.05, 0.1, 0.17, 0.6, 0.65, 0.13, 0.02, 0.05]',
    'PQSVECRPFVFGAGKPYEF':'PQSVECRPFVFGAGKPYEF [0.0, 0.3, 0.01, 0.0, 0.13, 0.03, 0.02, 0.03, 0.03, 0.0, 0.04, 0.02, 0.23, 0.0, 0.8, 0.76, 0.06, 0.14, 0.0] [0.89, 0.16, 0.62, 0.93, 0.55, 0.88, 0.67, 0.43, 0.95, 0.78, 0.57, 0.47, 0.16, 0.83, 0.07, 0.03, 0.34, 0.42, 0.89] [0.11, 0.54, 0.37, 0.07, 0.32, 0.09, 0.31, 0.54, 0.02, 0.22, 0.39, 0.51, 0.61, 0.17, 0.13, 0.21, 0.6, 0.44, 0.11]',
    'FRGVFAFLLYVATFMYVFS':'FRGVFAFLLYVATFMYVFS [0.08, 0.02, 0.09, 0.14, 0.22, 0.08, 0.19, 0.06, 0.03, 0.03, 0.07, 0.08, 0.03, 0.12, 0.06, 0.02, 0.04, 0.14, 0.1] [0.52, 0.27, 0.24, 0.52, 0.43, 0.64, 0.48, 0.5, 0.54, 0.64, 0.25, 0.16, 0.71, 0.46, 0.69, 0.69, 0.52, 0.39, 0.38] [0.4, 0.71, 0.67, 0.34, 0.35, 0.28, 0.33, 0.44, 0.43, 0.33, 0.68, 0.76, 0.26, 0.42, 0.25, 0.29, 0.44, 0.47, 0.52]',
    'YANYEGCLWNATGVVVCTG':'YANYEGCLWNATGVVVCTG [0.13, 0.13, 0.2, 0.22, 0.38, 0.01, 0.09, 0.11, 0.14, 0.19, 0.16, 0.13, 0.09, 0.13, 0.14, 0.11, 0.12, 0.18, 0.37] [0.2, 0.19, 0.32, 0.14, 0.08, 0.69, 0.21, 0.04, 0.33, 0.21, 0.31, 0.17, 0.17, 0.23, 0.07, 0.08, 0.18, 0.08, 0.1] [0.67, 0.68, 0.48, 0.64, 0.54, 0.3, 0.7, 0.85, 0.53, 0.6, 0.53, 0.7, 0.74, 0.64, 0.79, 0.81, 0.7, 0.74, 0.53]',
    'IATTVNLRDGQTLLLGGLT':'IATTVNLRDGQTLLLGGLT [0.15, 0.48, 0.02, 0.28, 0.01, 0.3, 0.04, 0.68, 0.17, 0.05, 0.45, 0.03, 0.09, 0.01, 0.03, 0.03, 0.01, 0.02, 0.23] [0.12, 0.08, 0.47, 0.1, 0.31, 0.26, 0.22, 0.15, 0.42, 0.21, 0.22, 0.48, 0.1, 0.31, 0.26, 0.61, 0.64, 0.41, 0.13] [0.73, 0.44, 0.51, 0.62, 0.68, 0.44, 0.74, 0.17, 0.41, 0.74, 0.33, 0.49, 0.81, 0.68, 0.71, 0.36, 0.35, 0.57, 0.64]',
    'LLDVNATTISRIDAAFSAR':'LLDVNATTISRIDAAFSAR [0.16, 0.11, 0.48, 0.31, 0.87, 0.54, 0.02, 0.38, 0.0, 0.08, 0.02, 0.0, 0.0, 0.01, 0.0, 0.06, 0.03, 0.03, 0.38] [0.13, 0.23, 0.26, 0.12, 0.06, 0.13, 0.92, 0.13, 0.55, 0.22, 0.68, 0.55, 0.91, 0.52, 0.87, 0.26, 0.9, 0.31, 0.34] [0.71, 0.66, 0.26, 0.57, 0.07, 0.33, 0.06, 0.49, 0.45, 0.7, 0.3, 0.45, 0.09, 0.47, 0.13, 0.68, 0.07, 0.66, 0.28]'
    }

replace_mtx = funtrp_line_to_mtx(baked_funtrp_lines[fragment])
log_replace_mtx = np.log(replace_mtx)
def funtrp_variant_log_prob(variant):
    return sum((log_replace_mtx[i,j] for j, i in enumerate(int_encode_variant(variant))))

def percent(c, denom=num_entries):
    return str(c/denom*100) + '%'

cache_dir = os.path.join('resources', 'cache')
def get_counts(fnames):
    combined_name = os.path.join(cache_dir, '--'.join(sorted((
            [os.path.basename(fname).split('.')[0] for fname in fnames]
            ))) + '.cache.json')
    try:
        with open(combined_name) as cache_f:
            return json.loads(cache_f.read())
    except FileNotFoundError:
        counts = fastq_fragment_counts(fnames, leading_seq='tttagatcgt',
                trailing_seq=None)#'gacgaaactc')
        with open(combined_name, 'w+') as cache_f:
            cache_f.write(json.dumps(counts, indent=4))
        return counts

if __name__ == '__main__':
    #load_all_data()
    #count_translations()
    #load_all_vars_and_transs()
    if making_figs:
        import matplotlib.pyplot as plt

    print('Fragment:', fragment)
    df = dataframes[fragment]
    print(df.keys())
    df_export = DataFrame()
    df_export['variant_aa_seq'] = df['variant_aa_seq']
    df_export['predicted_fitness'] = df_export.variant_aa_seq.apply(funtrp_variant_log_prob)
    quite_high = df_export.predicted_fitness.quantile(.97)

    def is_fit(fns):
        return 'TRUE' if fns > quite_high else 'FALSE'

    df_export['predicted_functional'] = df_export.predicted_fitness.apply(is_fit)
    print(df_export)
    df_export.to_csv(os.path.join(gen_f_dir, 'predictions_' + fragment + '.csv'))
    variants = df['variant_aa_seq']
    pred_fitnesses = [funtrp_variant_log_prob(variant) for variant in variants]
    wt_lfe = df.loc[df['variant_aa_seq'] == fragment, 'log_fold_enrichment'].item()
    current_haz_db = set(v.upper() for v in jsonf(os.path.join('resources', 'hazard_db', 'hazards.fragset.json')))
    print(fragment in current_haz_db)
    meas_fitnesses = [df.loc[df['variant_aa_seq'] == variant,
            'log_fold_enrichment'].item() - wt_lfe
            for variant in variants]
    print(len(df[df['variant_aa_seq'].isin(current_haz_db)]))
    print('out of')
    print(len(current_haz_db))
    meas_fitnesses_no_nans = np.array([mf for mf in meas_fitnesses if mf not in (float('nan'), float('-inf'), float('inf'))])
    min_meas = np.min(meas_fitnesses_no_nans)
    max_meas = np.max(meas_fitnesses_no_nans)
    # find the fitness cutoff for the 50th percentile
    meas_fitness_q50 = np.quantile(np.array(meas_fitnesses), .5)
    min_fitnesses = [-3.0, meas_fitness_q50] # -3 is nat log .05
    print('range of measurements:', min_meas, 'to', max_meas)
    if making_figs:
        plt.scatter(pred_fitnesses, meas_fitnesses, s=1)
        plt.xlabel('Predicted fitness (arbitrary units)')
        plt.ylabel('Measured fitness (fold enrichment)')
        # redo ticks on y-axis with natural-log-fold enrichment converted to base 10
        ticks = np.array([.01, .02, .05, .1, .2, .5, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000])/10
        plt.yticks(np.log(ticks), ticks)
        # bound the y-axis excluding outliers
        print(np.quantile(np.array(meas_fitnesses_no_nans), .001) - .7, np.quantile(np.array(meas_fitnesses_no_nans), .999) + .7)
        plt.ylim(np.quantile(np.array(meas_fitnesses_no_nans), .001) - .7, np.quantile(np.array(meas_fitnesses_no_nans), .999) + .7)
        # add horizontal lines for .05 and 50th quantile
        plt.plot([min(pred_fitnesses), max(pred_fitnesses)], [meas_fitness_q50, meas_fitness_q50], color='gray', linestyle='--')
        plt.plot([min(pred_fitnesses), max(pred_fitnesses)], [-3, -3], color='purple', linestyle='--', alpha=.3)
        # label lines. make sure the text has a white background behind it so it's readable. no line around the bounding box. No bleed around the text.
        plt.text(min(pred_fitnesses), meas_fitness_q50, '50th percentile', horizontalalignment='left', verticalalignment='bottom', bbox=dict(facecolor='white', alpha=0.85, edgecolor='none', pad=0))
        plt.text(min(pred_fitnesses), -3, '0.05 fitness', horizontalalignment='left', verticalalignment='top', bbox=dict(facecolor='white', alpha=0.85, edgecolor='none', pad=0))

        plt.savefig(os.path.join(gen_f_dir, 'corr_' + fragment + '.png'))
        plt.figure()
    print('num variants:', len(variants))

    def colors():
        yield 'purple', .6
        while True:
            yield 'gray', 1
    colors = colors()

    for min_fitness in min_fitnesses:
        print('generating ROC with fitness cutoff at', min_fitness)
        actual_positives = set()
        actual_negatives = set()
        for variant, pf, mf in zip(variants, pred_fitnesses, meas_fitnesses):
            if mf > min_fitness:
                actual_positives.add(variant)
            else:
                actual_negatives.add(variant)
        if len(actual_positives) * len(actual_negatives) == 0: # avoid errors, uninformative anyway
            print('skipped min_fitness', min_fitness)
            continue
        print('num real positives:', len(actual_positives))
        min_s = min(pred_fitnesses)
        max_s = max(pred_fitnesses)
        s_res = 100
        xs = []
        ys = []
        ss = []
        for s_i in range(s_res):
            s = s_i/s_res*(max_s - min_s) + min_s
            positives = set()
            for variant, pf, mf in zip(variants, pred_fitnesses, meas_fitnesses):
                if pf > s:
                    positives.add(variant)
            true_positives = positives & actual_positives
            false_positives = positives & actual_negatives
            tp_rate = len(true_positives)/len(actual_positives)
            fp_rate = len(false_positives)/len(actual_negatives)
            xs.append(fp_rate)
            ys.append(tp_rate)
            ss.append(s)
        
        csv_file = os.path.join(gen_f_dir, 'roc_' + fragment + '_' + str(round(np.exp(min_fitness)*1000)/1000) + '.csv')
        with open(csv_file, mode='w') as file:
            # use the dictionary writer to write the data
            writer = csv.DictWriter(file, fieldnames=['threshold', 'false_positive_rate', 'true_positive_rate'])
            writer.writeheader()
            for threshold, fp_rate, tp_rate in zip(ss, xs, ys):
                writer.writerow({'threshold': threshold, 'false_positive_rate': fp_rate, 'true_positive_rate': tp_rate})
        print('CSV file created:', csv_file)

        if making_figs:
            color, alpha = next(colors)
            plt.plot(xs, ys, label='Fitness cutoff: ' + str(round(np.exp(min_fitness)*100)/100), color=color, alpha=alpha)

    if making_figs:
        plt.title(fragment)
        # add a diagonal line for reference
        plt.plot([0,1], [0,1], color='lightgray', linestyle='--')
        # label horizontal and vertical axes appropriately for ROC
        plt.xlabel('False positive rate')
        plt.ylabel('True positive rate')
        plt.legend()
        plt.savefig(os.path.join(gen_f_dir, 'roc_' + fragment + '.png'))
        if showing_figs:
            plt.show()

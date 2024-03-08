import os, sys
import pandas as pd
import numpy as np
from Bio.Seq import Seq
from csv import DictReader, DictWriter

# the dict dataframes is the main result of running this script. it can be imported.

this_file_dir = os.path.dirname(__file__)
interim_file_dir = os.path.join(this_file_dir, 'interim-files')

prot_libs = ('Lib7', 'Lib8', 'Lib9', 'Lib10', 'Lib11', 'Lib12', 'Lib13', 'Lib14', 'Lib15')

frag_12 = 'PQSVECRPFVFGAGKPYEF'
other_orig_frags = ['FRGVFAFLLYVATFMYVFS', 'IATTVNLRDGQTLLLGGLT', 'LLDVNATTISRIDAAFSAR', 'MAVYFVTGKLGSGKTLVSV', 'NFYPCVEIKASPAKVLQGH', 'VEIKASPAKVLQGHNVFGT', 'YANYEGCLWNATGVVVCTG', 'YSYLTPYLSHGRYFKPLNL']
all_frags = [frag_12] + other_orig_frags
other_interm_fpaths = {frag:os.path.join(interim_file_dir, 'counts_' + frag + '.csv') for frag in other_orig_frags}

if '--make-interim-files' in sys.argv:
    interm_files = {}
    interm_writers = {}
    try:
        print('hi')
        with open(os.path.join(this_file_dir, '2021_02_08_multiple_libraries', '2021-03-29 compiled_library_counts.csv')) as f:
            reader = DictReader(f)
            for row in reader:
                if row['Lib'] not in prot_libs:
                    continue
                orig_frag = row['Comment'].split('.')[0].split('fragment')[1].strip()
                if orig_frag not in interm_files:
                    f = open(other_interm_fpaths[orig_frag], 'w+')
                    interm_files[orig_frag] = f
                    wr = DictWriter(f, fieldnames=('orig_fragment', 'variant_aa_seq', 'count_pre', 'count_post', 'trimmed'))
                    wr.writeheader()
                    interm_writers[orig_frag] = wr
                interm_wr = interm_writers[orig_frag]
                interm_wr.writerow({'orig_fragment':orig_frag, 'variant_aa_seq':str(Seq(row['trimmed']).translate()), 'count_pre':row['Pre'], 'count_post':row['Post'], 'trimmed':row['trimmed']})
    finally:
        for f in interm_files.values():
            print('closing', f.name)
            f.close()


if '--make-dataframes' in sys.argv:
    frag_12_dataframe = f12df = pd.read_csv(os.path.join(this_file_dir, '2020_12_10_synthesis_screening_lib_12', 'PQSVEC_variant_counts.csv'))
    f12df['orig_fragment'] = frag_12
    dataframes = {frag_12:f12df}
    for frag in other_orig_frags:
        try:
            dataframes[frag] = pd.read_csv(other_interm_fpaths[frag])
        except FileNotFoundError as e:
            raise type(e)('\n' + str(e) + "\n\nif temp count files don't exist yet, you may want to run this script (" + __file__ + ") with the --make-interim-files switch\n")

    # normalize to calculate abundance. This is the relative abundance of each library member 
    # (eliminates bias based on imbalanced number of reads)
    # each of these abundance_n data structures is a Pandas Series object & supports vector
    for i in range(1, 7):
        f12df['abundance_'+str(i)] = f12df['reads_brian_'+str(i)] / f12df['reads_brian_'+str(i)].sum()
    for oof in other_orig_frags:
        odf = dataframes[oof]
        odf['abundance_pre'] = odf['count_pre'] / odf['count_pre'].sum()
        odf['abundance_post'] = odf['count_post'] / odf['count_post'].sum()

    # calculate log fold enrichment relative to library 1
    for i in range(2, 7):
        f12df['log_fold_enrichment_'+str(i)] = np.log(f12df['abundance_'+str(i)] / f12df['abundance_1'])
    # use point 4, post-1-selection, as "default" count & log fold enrichment
    f12df['count_pre'] = f12df['reads_brian_1']
    f12df['count_post'] = f12df['reads_brian_4']
    f12df['abundance_pre'] = f12df['abundance_1']
    f12df['abundance_post'] = f12df['abundance_4']
    f12df['log_fold_enrichment'] = f12df['log_fold_enrichment_4']
    for oof in other_orig_frags:
        odf = dataframes[oof]
        invalid = (float('inf'), -float('inf'), 0, float('NaN'))
        filter_odf = odf[~(odf.abundance_post.isin(invalid)
                     | odf.abundance_pre.isin(invalid))]
        filter_odf['log_fold_enrichment'] = np.log(filter_odf['abundance_post'] / filter_odf['abundance_pre'])
        print('------------FROM NGS DATA PARSING SCRIPT---------', max(filter_odf['log_fold_enrichment']))
        dataframes[oof] = filter_odf

    for orig_frag, df in dataframes.items():
        def diff_letters(seq):
            return sum ( seq[i] != orig_frag[i] for i in range(len(seq)) )
    # translate the DNA lib member into protein
        df['variant_aa_seq'] = [str(Seq(x).translate()) for x in df.trimmed]
        # calculate how many mutations away each library member is
        df['point_mut_num'] = df.variant_aa_seq.apply(diff_letters)
        df.to_csv(os.path.join(interim_file_dir, 'finished_' + orig_frag + '.csv'))
else:
    try:
        dataframes = {frag:pd.read_csv(os.path.join(interim_file_dir, 'finished_' + frag + '.csv')) for frag in all_frags}
    except FileNotFoundError as e:
        raise type(e)('\n' + str(e) + "\n\nif finished dataframe files don't exist yet, you may want to run this script (" + __file__ + "), first with the --make-interim-files switch if necessary, and then  with the --make-dataframes switch\n")
# the dict dataframes is the main result of running this script. it can be imported.


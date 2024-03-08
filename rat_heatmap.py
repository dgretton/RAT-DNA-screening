from variant_tools import build_fitness_dict, funtrp_datas, mhv_for_fragment, mhv_likelihood
import numpy as np 
import sys

# provide switch --combinatorial to run the combinatorial escapee fraction calculation

fitness_dicts = {}

# function to only print a couple choice entries from a dictionary
def print_dict(d, n=3):
    if not d:
        print("Empty dictionary")
        return
    i = 0
    for k, v in d.items():
        print(k, v)
        i += 1
        if i >= n:
            break

def get_fitness_dict(fragment, dataframe=None):
    if fragment in fitness_dicts:
        return fitness_dicts[fragment]
    if dataframe is None:
        raise ValueError(f"Original fragment {fragment} not found in fitness_dicts and no dataframe provided")
    fd = build_fitness_dict(fragment, dataframe)
    fitness_dicts[fragment] = fd
    return fd

def fitness_bin_key(fitness):
    if fitness in (float('inf'), float('-inf'), float('nan')):
        return fitness
    return int(fitness * 1000) // 10

# build a dictionary of fitness bins to sets of variants not in the db
def escapees_by_fitness_bin(fragment, db_set):
    """Build a dictionary of fitness bins to sets of variants."""
    fragment_fitness_dict = get_fitness_dict(fragment)
    binned_escapees = {}
    for variant, fitness in fragment_fitness_dict.items():
        if variant in db_set:
            continue
        variant_bin_key = fitness_bin_key(fitness)
        if variant_bin_key not in binned_escapees:
            binned_escapees[variant_bin_key] = set()
        binned_escapees[variant_bin_key].add(variant)
    # make sure all the fitness bins add up to the total number of variants not in the db
    print(f"sum of escapees: {sum(len(variants) for variants in binned_escapees.values())}")
    print(f"total variants: {len(fragment_fitness_dict)}")
    print(f"total variants not in db: {len(set(fragment_fitness_dict) - db_set)}")

    assert sum(len(variants) for variants in binned_escapees.values()) == len(set(fragment_fitness_dict) - db_set)
    return binned_escapees

#defined like get_fitness_dict
binned_escapees_dicts = {}
def get_binned_escapees(fragment, db_set=None):
    if fragment in binned_escapees_dicts:
        return binned_escapees_dicts[fragment]
    if db_set is None:
        raise ValueError(f"Original fragment {fragment} not found in binned_escapees_dicts and no db_set provided")
    binned_escapees_dicts[fragment] = escapees_by_fitness_bin(fragment, db_set)
    return binned_escapees_dicts[fragment]

cum_escapees = {}
def cumulative_escapees_above(fragment, fbin_key, binned_escapees=None):
    if (fragment, fbin_key) in cum_escapees:
        return cum_escapees[(fragment, fbin_key)]
    escapees = 0
    for bin_key, variants in binned_escapees.items():
        if bin_key > fbin_key:
            escapees += len(variants)
    cum_escapees[(fragment, fbin_key)] = escapees
    return escapees

fit_variants_above = {}
def total_fit_variants_above(fragment, fitness):
    bin_key = fitness_bin_key(fitness)
    if (fragment, bin_key) in fit_variants_above:
        return fit_variants_above[(fragment, bin_key)]
    fragment_fitness_dict = get_fitness_dict(fragment)
    num_fit_variants = 0
    for fit in fragment_fitness_dict.values():
        if fit > fitness:
            num_fit_variants += 1
    fit_variants_above[(fragment, bin_key)] = num_fit_variants
    return num_fit_variants

def escapee_frac(fragment, db_set, threshold):
    """Calculate the fraction of escapees in a given database."""
    fragment_fitness_dict = get_fitness_dict(fragment)
    # all of the variants in fitness bins strictly higher than the threshold bin are escapees
    threshold_bin = fitness_bin_key(threshold)
    binned_escapees = get_binned_escapees(fragment, db_set)
    escapees = cumulative_escapees_above(fragment, threshold_bin, binned_escapees)
    total = total_fit_variants_above(fragment, threshold)
    # now we need to calculate the fraction of the variants in the threshold bin that are escapees
    # first, we need to find the variants in the threshold bin
    if threshold_bin not in binned_escapees:
        # there are no variants in the threshold bin, so the extra escapee count is 0
        return escapees / total
    threshold_bin_variants = binned_escapees[threshold_bin]
    for variant in threshold_bin_variants:
        fitness = fragment_fitness_dict[variant]
        if fitness > threshold:
            escapees += 1
    return escapees / total

def escapee_frac_inefficient(fragment, db_set, threshold):
    """Calculate the fraction of escapees in a given database."""
    fragment_fitness_dict = get_fitness_dict(fragment)
    escapees = 0
    total_fit_variants = 0
    for variant in fragment_fitness_dict:
        fitness = fragment_fitness_dict[variant]
        if fitness >= threshold:
            total_fit_variants += 1
            if variant not in db_set:
                escapees += 1
    return escapees / total_fit_variants

# a function that gives the product of all the escape fractions for a particular fitness.
def product_escapee_frac(fragments, db_set, threshold):
    """Calculate the product of the escapee fractions for a given list of thresholds."""
    product = 1
    for fragment in fragments:
        ef = escapee_frac(fragment, db_set, threshold)
        product *= ef
    return product

# utility, convert an array of log numbers from natural log to base 10. if it's broadcast-able, do that, and if not, just iterate and do it elementwise.
def log_to_base_10(log_array, base_change=1/np.log(10)):
    try:
        return log_array * base_change # should never accidentally list-multiply because base_change is a float, not an int
    except (ValueError, TypeError):
        return np.array([log * base_change for log in log_array])
    
def combinatorial_escapees(db_set, num_samples=100000, threshold=.05):
    # this parallelizable method works similarly to escapee_frac_inefficient, but instead of calculating the escapee fraction for each fragment, it calculates the escapee fraction for each combination of fragments and compares it to the threshold to get the fraction of combinations that are escapees
    import random
    # come up with a unique identifier for this thread
    unique_id = str(random.randint(100000000, 1000000000))
    verbose = False

    def printuid(output=True, *args):
        if output:
            print(f"{unique_id}: ", *args)

    printuid(True, f"num_samples: {num_samples}, threshold: {threshold}")
    log_threshold = np.log(threshold)
    printuid(True, f"log_threshold: {log_threshold}")
    
    # lists for each fragment
    variant_lists = []
    filtered_fitness_dicts = {}
    for fragment in funtrp_datas:
        # no point sampling fitnesses already smaller than the threshold because they'll never be functional in combination with any other variants, so we'll just sample from the fitnesses above the threshold
        # also prevent fitnesses from being larger than wild type, i.e. greater than 0
        filtered_fitness_dicts[fragment] = {variant: min(0, fitness) for variant, fitness in get_fitness_dict(fragment).items() if fitness >= log_threshold}
        variant_lists.append(list(filtered_fitness_dicts[fragment].keys()))
        # print the length and the first few values
        printuid(verbose, f"length for {fragment}: {len(variant_lists[-1])}")
        printuid(verbose, f"first few values for {fragment}: {variant_lists[-1][:5]}")
    # only brief summary if not verbose
    printuid(not verbose, f"lengths: {[len(variant_list) for variant_list in variant_lists]}")
    # make overall fitness dictionary for fast lookup
    fitness_dict = {}
    for fragment in funtrp_datas:
        fitness_dict.update(filtered_fitness_dicts[fragment])
    escapees = 0
    fit_variants = 0
    for i in range(num_samples):
        # benchmark progress
        if i % 1000000 == 0:
            printuid(True, f"[Sample {i}]")
            printuid(not verbose, f"escapees so far: {escapees}, fit variants so far: {fit_variants}")
        # choose a random fragment from each set of fragments
        random_variants = [random.choice(variant_list) for variant_list in variant_lists]
        random_lfes = [fitness_dict[variant] for variant in random_variants]
        # calculate the product of the log random numbers
        product = sum(random_lfes)
        # if the product is greater than the threshold, increment the escapee count
        if product > log_threshold:
            # this is a fit variant but not an escapee yet; it's an escapee if it isn't in the db set
            printuid(verbose, f"fit variant: {random_variants} ({random_lfes})")
            fit_variants += 1
            printuid(verbose, f"fit variants so far: {fit_variants}")
            if not any(variant in db_set for variant in random_variants):
                printuid(verbose, f"ESCAPEE: {random_variants} ({random_lfes})")
                escapees += 1
                printuid(verbose, f"escapees so far: {escapees}")
            printuid(verbose, f"(iteration {i})")
    printuid(True, f"escapees: {escapees}, fit variants: {fit_variants}, num samples: {num_samples}")
    return escapees

if __name__ == '__main__':
    from analyze_ngs import (
        dataframes,
        )
    import matplotlib.pyplot as plt
    from matplotlib.cm import get_cmap
    plot_resolution = 500
    fragments = list(funtrp_datas)
    fragment_nums = {f: i for i, f in enumerate(fragments)}
    # every time we use log scale we'll convert to base 10 before plotting
    # initialize a (num fragments) x (resolution) grid named escape_fracs where the resolution dimension is populated with (resolution) steps between -8 and 0
    escape_fracs = np.zeros((len(fragments), plot_resolution))
    thresholds = np.linspace(-8, 0, plot_resolution)
    plot_thresholds = log_to_base_10(thresholds)
    def setcolors():
        plt.gca().set_prop_cycle(color=get_cmap('tab20').colors)
    plt.figure()
    setcolors()
    db_size = 50000
    combined_db_set = set()
    import pandas as pd
    only_saving_csvs = '--csv' in sys.argv
    if only_saving_csvs:
        simple_fitness_csv_df = pd.DataFrame() # this will be a dataframe with the measured fitnesses of all the variants in the combined database, and also the calculated/estimated fitness from mhv_likelihood. there should be a column for the original fragment that each variant is derived from.
        escapees_by_fitness_threshold_df = pd.DataFrame() # this will be a dataframe with the escapee fraction for each fragment at each fitness threshold

    for fragment in fragments:
        # Construct the fitness dictionary
        df = dataframes.get(fragment)
        if df is None:
            raise ValueError(f"Original fragment {fragment} not found in data")
        fragment_fitness_dict = get_fitness_dict(fragment, df) # this will build the fitness dictionary if it doesn't already exist
        # files named like "[fragment]_sorted_variants.txt" (e.g. YSYLTPYLSHGRYFKPLNL_sorted_variants.txt) can be found in the directory resources/screening_databases/ 
        # its entries look like this:
        # YSYLTPYLSHGRYFKPLNL	0.30454089050927763
        # YTYLTPYLSHGRYFKPLNL	0.3013848028344842
        # YNYLTPYLSHGRYFKPLNL	0.3013848028344842
        # etc.
        # the numbers are "coverage," not relevant to this test
        # read the first column (variants) into a set
        db_set = set()
        with open(f"resources/screening_databases/{fragment}_sorted_variants.txt") as f:
            lines = list(f)
            for _, line in zip(range(db_size), lines):
                db_set.add(line.split()[0])
        print(f"Number of variants in db: {len(db_set)}")
        combined_db_set.update(db_set)

        if only_saving_csvs:
            # add the fitnesses and fragment names to the simple_fitness_csv_df dataframe
            mhv = mhv_for_fragment(fragment)
            variants_in_db_or_measured = list(set(fragment_fitness_dict.keys()) | db_set) # no duplicates
            def variant_fitness_log10(variant): # fitnesses are in natural logs, meaning they're not easily interpretable. also, they're mostly negative.
                # for securedna manuscript we want to plot them in base 10, so we'll convert them here
                if variant in fragment_fitness_dict:
                    return log_to_base_10(fragment_fitness_dict[variant])
                return 'no measurement'
            simple_fitness_csv_df = simple_fitness_csv_df.append(
                    pd.DataFrame({'variant': variants_in_db_or_measured,
                                  'fitness_log10': [variant_fitness_log10(variant) for variant in variants_in_db_or_measured],
                                  'window': fragment, # original fragments are windows by the terminology in the securedna manuscript
                                  'predicted_fitness': [mhv_likelihood(variant, mhv) for variant in variants_in_db_or_measured],
                                  'in_database': [variant in db_set for variant in variants_in_db_or_measured]
                                  }), ignore_index=True)
        
        x = []
        y = []
        for threshold in np.linspace(-8, 0, 10):
            x.append(threshold)
            y.append(escapee_frac_inefficient(fragment, db_set, threshold))
        x = log_to_base_10(x) # for securedna manuscript we want to plot in base 10

        for i, threshold in enumerate(thresholds):
            escape_fracs[fragment_nums[fragment]][i] = escapee_frac(fragment, db_set, threshold)
        plt.plot(plot_thresholds, escape_fracs[fragment_nums[fragment], :], label=fragment)
        if only_saving_csvs:
            escapees_by_fitness_threshold_df = escapees_by_fitness_threshold_df.append(
                pd.DataFrame({'window': fragment,
                              'fitness_threshold': plot_thresholds,
                              'escapee_fraction': escape_fracs[fragment_nums[fragment], :]
                              }), ignore_index=True)
    
    if only_saving_csvs:
        # reorder columns: window, variant, predicted_fitness, fitness
        simple_fitness_csv_df = simple_fitness_csv_df[['window', 'variant', 'predicted_fitness', 'fitness_log10', 'in_database']]
        simple_fitness_csv_df.to_csv(f'resources/screening_databases/predicted_and_measured_fitnesses_db{db_size}.csv', index=False)
        escapees_by_fitness_threshold_df.to_csv(f'resources/screening_databases/escapees_by_fitness_threshold_db{db_size}.csv', index=False)
        exit()

    if '--combinatorial' in sys.argv:
        from multiprocessing import Pool, cpu_count
        cpu_count = cpu_count()
        print(f"cpu_count: {cpu_count}")
        def combinatorial_escapees_func(num_samples):
            return combinatorial_escapees(combined_db_set, num_samples=num_samples, threshold=.05)
        # progressively larger sample numbers for combinatorial escape fraction
        for sample_num in [1000, 10000, 100000, 1000000000, 5000000000, 10000000000]:
            with Pool(cpu_count) as p:
                samples_per_thread = sample_num // cpu_count
                print(f"samples_per_thread: {samples_per_thread}")
                total_combinatorial_escapees = sum(p.map(combinatorial_escapees_func, [samples_per_thread] * cpu_count))
            print(f"FINISHED POOL MAP CALL FOR COMBINATORIAL ESCAPEES WITH SAMPLE_NUM {sample_num}")
            # use 3 significant figures for the escapee fraction
            print(f"sample num: {sample_num}, combinatorial escapees: {total_combinatorial_escapees}, fraction: {total_combinatorial_escapees / sample_num:.3g}\n\n\n")

    plt.plot(plot_thresholds, np.product(escape_fracs, axis=0), label='Product', linewidth=3)
    # the following two lines "plot" the minimum and maximum of the product of the escapee fractions, just one point each
    plt.scatter(plot_thresholds[np.argmin(np.product(escape_fracs, axis=0))], min(np.product(escape_fracs, axis=0)), color='red', label=f'Product minimum ({min(np.product(escape_fracs, axis=0)):.2E})') # show the product minimum as a number formatted in scientific notation with 3 digits past the decimal point
    plt.scatter(plot_thresholds[np.argmax(np.product(escape_fracs, axis=0))], max(np.product(escape_fracs, axis=0)), color='blue', label=f'Product maximum ({max(np.product(escape_fracs, axis=0)):.2g})')
    # also plot the geometric mean across all of the fragments for each of plot_resolution divisions of the range of fitness -8 to 0
    geom_mean = np.product(escape_fracs, axis=0) ** (1/len(fragments))
    # also calculate geometric mean excluding VEIKA, which may be an outlier
    geom_mean_no_veika = np.product(np.delete(escape_fracs, fragment_nums['VEIKASPAKVLQGHNVFGT'], axis=0), axis=0) ** (1/(len(fragments) - 1))
    # bold the geometric mean so it stands out
    plt.plot(plot_thresholds, geom_mean, label='Geometric mean', linewidth=4, linestyle='--')
    plt.plot(plot_thresholds, geom_mean_no_veika, label='Geometric mean (no VEIKA)', linewidth=4, linestyle='--')
    plt.legend()
    plt.xlabel('Fitness threshold')
    plt.ylabel('Fraction of escapees')
    plt.title('Escapee fraction vs. fitness threshold')
    plt.gca().xaxis.set_major_locator(plt.MultipleLocator(1))
    # format y scale labels as 10.0, 1.0, 0.1, 0.01, i.e. 10^(actual number on the plot) showing exactly 2 digits past the decimal point
    plt.gca().xaxis.set_major_formatter(plt.FuncFormatter(lambda x, pos: f"{10**x:.2g}"))
    # New figure. Contour plot showing powers of the geometric mean (or whatever other means of combining the various escapees vs fitness) horizontally and escapee fractions vertically. Include powers up to n_windows. a general term for things like geom_mean is "aggregated escapee fraction"
    def plot_contour(aggregated_escapee_fraction, aggregation_method_label):
        ax = plt.figure().add_subplot()

        # Nines is a list of floats .9, .99, .999, etc.
        n_levels = 6
        nines = [float('0.' + '9' * i) for i in range(1, n_levels + 1)]
        # calculate number of windows necessary to hit highest escapee fraction in nines
        n_windows = int(np.ceil(np.log(1 - nines[-1]) / np.log(aggregated_escapee_fraction.max()))) # make the width of the plot approximately match the expected extent of the rightmost countour line. typically between 30 and 50.
        powers = np.arange(0, n_windows + 1)
        # Create a meshgrid to compute Z, the escapee fraction
        P, Y = np.meshgrid(powers, plot_thresholds)
        Z = 1 - aggregated_escapee_fraction[:, None] ** P

        fig = plt.gcf()
        fig.set_size_inches(1300 / fig.dpi, 800 / fig.dpi)
        n_minor_levels = 5
        # space n_minor_levels values between each major level (given by nines) log spaced. leave out the elements of nines themselves so we don't overdraw.
        # first, take the logs of the differences of the nines array from 1, so spacing can be linear
        log_nines = np.log(1 - np.array(nines))
        minor_levels = []
        minor_levels.extend(np.linspace(0, log_nines[0], n_minor_levels, endpoint=False)[1:])
        for i in range(len(nines) - 1):
            minor_levels.extend(np.linspace(log_nines[i], log_nines[i+1], n_minor_levels, endpoint=False)[1:])
        
        minor_levels = 1 - np.exp(minor_levels)

        # Choose a colormap
        # import matplotlib.cm as cm; cmap = cm.get_cmap('winter') # original for colored contour figure
        # >>> revision for no-veika comparison
        # define a function that takes in 2 colors and returns a cmap gradient
        import matplotlib.colors as colors
        def cmap(start_color, end_color):
            return colors.LinearSegmentedColormap.from_list('custom', [start_color, end_color], N=256)
        cmap = cmap('#DDDDDD', 'white')
        # revision for no-veika comparison <<<

        # Generate the list of colors
        colors = [cmap(pos) for pos in np.linspace(0, 1, len(nines) + len(minor_levels) + 2)]
        # Plot the contour
        CF = ax.contourf(P, Y, Z, levels=sorted([0, 1] + list(nines) + list(minor_levels)), colors=colors)

        # linecolor = 'white'
        # >>> revision for no-veika comparison
        linecolor = 'black'
        # revision for no-veika comparison <<<
        CS_major = ax.contour(P, Y, Z, levels=nines, colors=linecolor)
        # no labels on minor contours
        CS_minor = ax.contour(P, Y, Z, levels=minor_levels, colors=linecolor, linewidths=0.5)

        # This custom formatter removes trailing zeros, e.g. "1.0" becomes "1", and
        # then adds a percent sign.
        def fmt(x):
            s = f"{x*100:.9f}"
            while s.endswith("0"):
                s = s[:-1]
            if s.endswith("."):
                s = s[:-1]
            return rf"   {s}\%   " if plt.rcParams["text.usetex"] else f"    {s}%   "
        
        # dotted horizontal line at fitness = -3 labeled "Fitness .05"
        # reference viruses: 1918 influenza (R0~2.5), SARS-CoV-2 (R0~4), mumps (R0~10), and measles (R0~18)
        pathogen_cutoffs = [('1918 influenza', 1/2.5), ('SARS-CoV-2', 1/4), ('Mumps', 1/10), ('Measles', 1/18)]
        # use smallest cutoff to set label offset
        label_offset = log_to_base_10(np.log(min(pathogen_cutoffs, key=lambda x: x[1])[1])) * .08 # fyi, is negative
        for path_name, pathogen_cutoff in pathogen_cutoffs:
            line_pos = log_to_base_10(np.log(pathogen_cutoff))
            ax.axhline(y=line_pos, linestyle='--', color=linecolor)
            ax.text(.97*n_windows, line_pos + label_offset, f"{path_name} (fitness {pathogen_cutoff:.2f})", horizontalalignment='right', color='black')
        # bold text on contour lines
        ax.clabel(CS_major, inline=True, fmt=fmt, fontsize=10)
        for txt in CS_major.labelTexts:
            txt.set_weight('bold')
        ax.set_ylim(min(plot_thresholds), max(plot_thresholds))
        ax.set_xlim(1, n_windows)
        # 2:1 window aspect ratio (wide)
        ax.set_aspect(.5 * (n_windows - 1) / log_to_base_10(8))
        # use a "grid" (really just a bunch of vertical lines) to show every horizontal power (e.g. there would be 42 vertical lines if n_windows = 43)
        # add minor ticks at every integer
        ax.xaxis.set_minor_locator(plt.MultipleLocator(1))
        # add major ticks at every 5
        ax.xaxis.set_major_locator(plt.MultipleLocator(5))
        ax.grid(True, axis='x', linewidth=0.5, which='both')
        ax.set_xlabel('Database size in windows, 50,000 variants each')
        # only major ticks on y axis, only on integers, using null locator
        # ax.yaxis.set_major_locator(plt.NullLocator())
        ax.yaxis.set_major_locator(plt.MultipleLocator(1))
        # format y scale labels as 10.0, 1.0, 0.1, 0.01, i.e. 10^(actual number on the plot) showing exactly 2 digits past the decimal point
        ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, pos: f"{10**x:.2g}"))
        ax.set_ylabel('Fitness threshold')
        ax.set_title(f"{aggregation_method_label.title()}")


    plot_contour(geom_mean, 'geometric mean')
    plot_contour(geom_mean_no_veika, 'geometric mean (no VEIKA)')
    # plot_contour(np.mean(escape_fracs, axis=0), 'mean')
    # plot_contour(np.median(escape_fracs, axis=0), 'median')
    # plot_contour(np.max(escape_fracs, axis=0), 'max')

    plt.show()

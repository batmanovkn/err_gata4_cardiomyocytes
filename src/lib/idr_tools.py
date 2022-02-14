import os
import math
from collections import namedtuple, defaultdict, OrderedDict
from itertools import chain
import sys

import pandas as pd
import numpy as np
import idr
from idr.optimization import estimate_model_params, old_estimator
from idr.utility import calc_post_membership_prbs, compute_pseudo_values
from scipy.stats.stats import rankdata

import matplotlib.style
import matplotlib as mpl
# Say, "the default sans-serif font is Arial"
matplotlib.rcParams['font.sans-serif'] = "Arial"
# Then, "ALWAYS use sans-serif fonts"
matplotlib.rcParams['font.family'] = "sans-serif"
mpl.rcParams['font.size'] = 15
mpl.rcParams['axes.spines.top'] = False
mpl.rcParams['axes.spines.right'] = False
mpl.rcParams['legend.frameon'] = False
mpl.rcParams['legend.loc'] = "upper left"
mpl.rcParams['figure.facecolor'] = "white"

from . import misc, parallel_run, bedtools, homer_tools


def peaks_to_narrowpeak(peaks, out=None):
    if out is None:
        out = misc.replace_extension(peaks, "narrowPeak")

    peaks_data = pd.read_csv(peaks, sep="\t", header=None, comment="#")
    peaks2_data = peaks_data[[1, 2, 3, 4, 6, 4, 7, 6, 10, 5]].copy()

    peaks2_data.loc[:, 2] -= 1
    peaks2_data[5] = peaks2_data[5].round().astype(int)
    peaks2_data.to_csv(out, sep="\t", header=False, index=False)
    return out


def idr_analysis(peaks1, peaks2, output_file, plot_output=None, runner=parallel_run.SimpleRunner()):
    runner.add(f"python -m lib.cl peaks_to_narrowpeak {misc.join_quoted_paths(peaks1)}")
    runner.add(f"python -m lib.cl peaks_to_narrowpeak {misc.join_quoted_paths(peaks2)}")
    narrowpeak1 = misc.replace_extension(peaks1, "narrowPeak")
    narrowpeak2 = misc.replace_extension(peaks2, "narrowPeak")
    plot_command = "--plot" * (plot_output is not None)
    runner.add(f"idr --samples {misc.join_quoted_paths([narrowpeak1, narrowpeak2])} -o {output_file} {plot_command}")
    if plot_output is not None:
        misc.rename_file(output_file + ".png", plot_output, runner)


def pseudo_idr(tag_dir, peaks_folder, output_file, plot_output=None, peak_calling_params=None):
    if peak_calling_params is None:
        peak_calling_params = {}

    pseudo_dirs = [tag_dir + f".pseudo{i}" for i in [1, 2]]
    homer_tools.pseudoreplicate(tag_dir, *pseudo_dirs)
    pseudo_peaks = [os.path.join(peaks_folder, os.path.basename(pd) + ".peaks") for pd in pseudo_dirs]
    parallel_run.parallel_run(
        lambda: homer_tools.find_peaks(pseudo_dirs[0], pseudo_peaks[0], parallel_run.SimpleRunner(), **peak_calling_params),
        lambda: homer_tools.find_peaks(pseudo_dirs[1], pseudo_peaks[1], parallel_run.SimpleRunner(), **peak_calling_params))

    idr_analysis(*pseudo_peaks, output_file, plot_output)


Peak = namedtuple(
    'Peak', ['chrm', 'strand', 'start', 'stop', 'signal', 'summit', 'signalValue', 'pValue', 'qValue'])

MergedPeak = namedtuple(
    'Peak', ['chrm', 'strand', 'start', 'stop', 'summit',
             'merged_signal', 'signals', 'pks'])


def load_bed(fp, signal_index, peak_summit_index=None):
    grpd_peaks = defaultdict(list)
    for line in fp:
        if line.startswith("#"): continue
        if line.startswith("track"): continue
        data = line.split()
        #print(data, signal_index)
        signal = float(data[signal_index])
        if idr.ONLY_ALLOW_NON_NEGATIVE_VALUES and signal < 0:
            raise ValueError("Invalid Signal Value: {:e}".format(signal))
        if peak_summit_index == None or int(data[peak_summit_index]) == -1:
            summit = None
        else:
            summit = int(data[peak_summit_index]) + int(float(data[1]))
        assert summit == None or summit >= 0
        peak = Peak(data[0], data[5],
                    int(float(data[1])), int(float(data[2])),
                    signal, summit,
                    float(data[6]), float(data[7]), float(data[8])
        )
        grpd_peaks[(peak.chrm, peak.strand)].append(peak)
    return grpd_peaks


def mean(items):
    items = list(items)
    return sum(items)/float(len(items))


def iter_merge_grpd_intervals(
        intervals, n_samples, pk_agg_fn,
        use_oracle_pks, use_nonoverlapping_peaks):
    # grp peaks by their source, and calculate the merged
    # peak boundaries
    grpd_peaks = OrderedDict([(i + 1, []) for i in range(n_samples)])
    pk_start, pk_stop = 1e12, -1
    for interval, sample_id in intervals:
        # if we've provided a unified peak set, ignore any intervals that
        # don't contain it for the purposes of generating the merged list
        if (not use_oracle_pks) or sample_id == 0:
            pk_start = min(interval.start, pk_start)
            pk_stop = max(interval.stop, pk_stop)
        # if this is an actual sample (ie not a merged peaks)
        if sample_id > 0:
            grpd_peaks[sample_id].append(interval)

    # if there are no identified peaks, continue (this can happen if
    # we have a merged peak list but no merged peaks overlap sample peaks)
    if pk_stop == -1:
        return None

    # skip regions that dont have a peak in all replicates
    if not use_nonoverlapping_peaks:
        if any(0 == len(peaks) for peaks in grpd_peaks.values()):
            return None

    # find the merged peak summit
    # note that we can iterate through the values because
    # grpd_peaks is an ordered dict
    replicate_summits = []
    for sample_id, pks in grpd_peaks.items():
        # if an oracle peak set is specified, skip the replicates
        if use_oracle_pks and sample_id != 0:
            continue

        # initialize the summit to the first peak
        try:
            replicate_summit, summit_signal = pks[0].summit, pks[0].signal
        except IndexError:
            replicate_summit, summit_signal = None, -1e9
        # if there are more peaks, take the summit that corresponds to the
        # replicate peak with the highest signal value
        for pk in pks[1:]:
            if pk.summit != None and pk.signal > summit_signal:
                replicate_summit, summit_signal = pk.summit, pk.signal
        # make sure a peak summit was specified
        if replicate_summit != None:
            replicate_summits.append(replicate_summit)

    summit = (int(mean(replicate_summits))
              if len(replicate_summits) > 0 else None)

    # note that we can iterate through the values because
    # grpd_peaks is an ordered dict
    signals = [pk_agg_fn(pk.signal for pk in pks) if len(pks) > 0 else 0
               for pks in grpd_peaks.values()]
    merged_pk = (pk_start, pk_stop, summit,
                 pk_agg_fn(signals), signals, grpd_peaks)

    yield merged_pk
    return


def iter_matched_oracle_pks(
        pks, n_samples, pk_agg_fn, use_nonoverlapping_peaks=False):
    """Match each oracle peak to it nearest replicate peaks.
    """
    oracle_pks = [pk for pk, sample_id in pks
                  if sample_id == 0]
    # if there are zero oracle peaks in this
    if len(oracle_pks) == 0: return None
    # for each oracle peak, find score the replicate peaks
    for oracle_pk in oracle_pks:
        peaks_and_scores = OrderedDict([(i + 1, []) for i in range(n_samples)])
        for pk, sample_id in pks:
            # skip oracle peaks
            if sample_id == 0: continue

            # calculate the distance between summits, setting it to a large
            # value in case the peaks dont have summits
            summit_distance = sys.maxsize
            if oracle_pk.summit != None and pk.summit != None:
                summit_distance = abs(oracle_pk.summit - pk.summit)
            # calculate the fraction overlap witht he oracle peak
            overlap = (1 + min(oracle_pk.stop, pk.stop)
                       - max(oracle_pk.start, pk.start))
            overlap_frac = overlap / (oracle_pk.stop - oracle_pk.start + 1)

            peaks_and_scores[sample_id].append(
                ((summit_distance, -overlap_frac, -pk.signal), pk))

        # skip regions that dont have a peak in all replicates.
        if not use_nonoverlapping_peaks and any(
                0 == len(peaks) for peaks in peaks_and_scores.values()):
            continue

        # build the aggregated signal value, which is jsut the signal value
        # of the replicate peak witgh the closest match
        signals = []
        rep_pks = []
        for rep_id, scored_pks in peaks_and_scores.items():
            scored_pks.sort()
            if len(scored_pks) == 0:
                assert use_nonoverlapping_peaks
                signals.append(0)
                rep_pks.append(None)
            else:
                signals.append(scored_pks[0][1].signal)
                rep_pks.append([scored_pks[0][1], ])
        all_peaks = [oracle_pk, ] + rep_pks
        new_pk = (oracle_pk.start, oracle_pk.stop, oracle_pk.summit,
                  pk_agg_fn(signals),
                  signals,
                  OrderedDict(zip(range(len(all_peaks)), all_peaks)))
        yield new_pk

    return


def merge_peaks_in_contig(all_s_peaks, pk_agg_fn, oracle_pks=None,
                          use_nonoverlapping_peaks=False):
    """Merge peaks in a single contig/strand.

    returns: The merged peaks.
    """
    # merge and sort all peaks, keeping track of which sample they originated in
    oracle_pks_iter = []
    if oracle_pks != None:
        oracle_pks_iter = oracle_pks

    # merge and sort all of the intervals, leeping track of their source
    all_intervals = []
    for sample_id, peaks in enumerate([oracle_pks_iter, ] + all_s_peaks):
        all_intervals.extend((pk, sample_id) for pk in peaks)
    all_intervals.sort()

    # grp overlapping intervals. Since they're already sorted, all we need
    # to do is check if the current interval overlaps the previous interval
    grpd_intervals = [[], ]
    curr_start, curr_stop = all_intervals[0][:2]
    for pk, sample_id in all_intervals:
        if pk.start < curr_stop:
            curr_stop = max(pk.stop, curr_stop)
            grpd_intervals[-1].append((pk, sample_id))
        else:
            curr_start, curr_stop = pk.start, pk.stop
            grpd_intervals.append([(pk, sample_id), ])

    # build the unified peak list, setting the score to
    # zero if it doesn't exist in both replicates
    merged_pks = []
    if oracle_pks == None:
        for intervals in grpd_intervals:
            for merged_pk in iter_merge_grpd_intervals(
                    intervals, len(all_s_peaks), pk_agg_fn,
                    use_oracle_pks=(oracle_pks != None),
                    use_nonoverlapping_peaks=use_nonoverlapping_peaks):
                merged_pks.append(merged_pk)
    else:
        for intervals in grpd_intervals:
            for merged_pk in iter_matched_oracle_pks(
                    intervals, len(all_s_peaks), pk_agg_fn,
                    use_nonoverlapping_peaks=use_nonoverlapping_peaks):
                merged_pks.append(merged_pk)

    return merged_pks


def merge_peaks(all_s_peaks, pk_agg_fn, oracle_pks=None,
                use_nonoverlapping_peaks=False):
    """Merge peaks over all contig/strands

    """
    # if we have reference peaks, use its contigs: otherwise use
    # the union of the replicates contigs
    if oracle_pks != None:
        contigs = sorted(oracle_pks.keys())
    else:
        contigs = sorted(set(chain(*[list(s_peaks.keys()) for s_peaks in all_s_peaks])))

    merged_peaks = []
    for key in contigs:
        # check to see if we've been provided a peak list and, if so,
        # pass it down. If not, set the oracle peaks to None so that
        # the callee knows not to use them
        if oracle_pks != None:
            contig_oracle_pks = oracle_pks[key]
        else:
            contig_oracle_pks = None

        # since s*_peaks are default dicts, it will never raise a key error,
        # but instead return an empty list which is what we want
        merged_contig_peaks = merge_peaks_in_contig(
            [s_peaks[key] for s_peaks in all_s_peaks],
            pk_agg_fn, contig_oracle_pks,
            use_nonoverlapping_peaks=use_nonoverlapping_peaks)
        merged_peaks.extend(
            MergedPeak(*(key + pk)) for pk in merged_contig_peaks)

    merged_peaks.sort(key=lambda x: x.merged_signal, reverse=True)
    return merged_peaks


def load_idr_samples(samples):
    # decide what aggregation function to use for peaks that need to be merged
    signal_type = 'signal.value'

    signal_index = {"score": 4, "signal.value": 6,
                    "p.value": 7, "q.value": 8}[signal_type]

    if signal_index in (4, 6):
        peak_merge_fn = sum
    else:
        peak_merge_fn = min

    summit_index = 9
    f1, f2 = [load_bed(fp, signal_index, summit_index)
              for fp in samples]
    oracle_pks = None
    # build a unified peak set
    idr.log("Merging peaks", 'VERBOSE')
    merged_peaks = merge_peaks([f1, f2], peak_merge_fn,
                               oracle_pks, False)
    return merged_peaks, signal_type


def plot_idr(scores, ranks, IDRs, output_plot_file):
    soft_idr_threshold = idr.DEFAULT_SOFT_IDR_THRESH
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot

    colors = np.zeros(len(ranks[0]), dtype=str)
    colors[:] = 'k'
    colors[IDRs > soft_idr_threshold] = 'r'

    # matplotlib.rc('font', family='normal', weight='bold', size=10)

    fig = matplotlib.pyplot.figure(num=None, figsize=(12, 12))

    matplotlib.pyplot.subplot(221)
    matplotlib.pyplot.axis([0, 1, 0, 1])
    matplotlib.pyplot.xlabel("Sample 1 Rank")
    matplotlib.pyplot.ylabel("Sample 2 Rank")
    matplotlib.pyplot.title(
        "Ranks - (red >= %.2f IDR)" % soft_idr_threshold)
    matplotlib.pyplot.scatter((ranks[0] + 1) / float(max(ranks[0]) + 1),
                              (ranks[1] + 1) / float(max(ranks[1]) + 1),
                              edgecolor=colors,
                              c=colors,
                              alpha=0.05)

    matplotlib.pyplot.subplot(222)
    matplotlib.pyplot.xlabel("Sample 1 log10 Score")
    matplotlib.pyplot.ylabel("Sample 2 log10 Score")
    matplotlib.pyplot.title(
        "Log10 Scores - (red >= %.2f IDR)" % soft_idr_threshold)
    matplotlib.pyplot.scatter(np.log10(scores[0] + 1),
                              np.log10(scores[1] + 1),
                              edgecolor=colors,
                              c=colors,
                              alpha=0.05)

    def make_boxplots(sample_i):
        groups = defaultdict(list)
        norm_ranks = (ranks[sample_i] + 1) / float(max(ranks[sample_i]) + 1)
        for rank, idr_val in zip(norm_ranks, -np.log10(IDRs)):
            groups[int(20 * rank)].append(float(idr_val))
        group_labels = sorted((x + 2.5) / 20 for x in groups.keys())
        groups = [x[1] for x in sorted(groups.items())]

        matplotlib.pyplot.title(
            "Sample %i Ranks vs IDR Values" % (sample_i + 1))
        matplotlib.pyplot.axis([0, 21,
                                0, 0.5 - math.log10(idr.CONVERGENCE_EPS_DEFAULT)])
        matplotlib.pyplot.xlabel("Sample %i Peak Rank" % (sample_i + 1))
        matplotlib.pyplot.ylabel("-log10 IDR")
        matplotlib.pyplot.xticks([], [])

        matplotlib.pyplot.boxplot(groups, sym="")

        matplotlib.pyplot.axis([0, 21,
                                0, 0.5 - math.log10(idr.CONVERGENCE_EPS_DEFAULT)])
        matplotlib.pyplot.scatter(20 * norm_ranks,
                                  -np.log10(IDRs),
                                  alpha=0.01, c='black')

    matplotlib.pyplot.subplot(223)
    make_boxplots(0)

    matplotlib.pyplot.subplot(224)
    make_boxplots(1)

    matplotlib.pyplot.savefig(output_plot_file)


def build_rank_vectors(merged_peaks):
    # allocate memory for the ranks vector
    s1 = np.zeros(len(merged_peaks))
    s2 = np.zeros(len(merged_peaks))
    # add the signal
    for i, x in enumerate(merged_peaks):
        s1[i], s2[i] = x.signals

    rank1 = np.lexsort((np.random.random(len(s1)), s1)).argsort()
    rank2 = np.lexsort((np.random.random(len(s2)), s2)).argsort()

    return (np.array(rank1, dtype=np.int),
            np.array(rank2, dtype=np.int))


def calc_local_IDR(theta, r1, r2):
    """
    idr <- 1 - e.z
    o <- order(idr)
    idr.o <- idr[o]
    idr.rank <- rank(idr.o, ties.method = "max")
    top.mean <- function(index, x) {
        mean(x[1:index])
    }
    IDR.o <- sapply(idr.rank, top.mean, idr.o)
    IDR <- idr
    IDR[o] <- IDR.o
    """
    mu, sigma, rho, p = theta
    z1 = compute_pseudo_values(r1, mu, sigma, p, EPS=1e-12)
    z2 = compute_pseudo_values(r2, mu, sigma, p, EPS=1e-12)
    localIDR = 1-calc_post_membership_prbs(np.array(theta), z1, z2)
    if idr.FILTER_PEAKS_BELOW_NOISE_MEAN:
        localIDR[z1 + z2 < 0] = 1

    # it doesn't make sense for the IDR values to be smaller than the
    # optimization tolerance
    localIDR = np.clip(localIDR, idr.CONVERGENCE_EPS_DEFAULT, 1)
    return localIDR


def fit_model_and_calc_local_idr(r1, r2,
                                 starting_point=None,
                                 max_iter=idr.MAX_ITER_DEFAULT,
                                 convergence_eps=idr.CONVERGENCE_EPS_DEFAULT,
                                 fix_mu=False, fix_sigma=False):
    # in theory we would try to find good starting point here,
    # but for now just set it to somethign reasonable
    if starting_point is None:
        starting_point = (idr.DEFAULT_MU, idr.DEFAULT_SIGMA,
                          idr.DEFAULT_RHO, idr.DEFAULT_MIX_PARAM)

    idr.log("Initial parameter values: [%s]" % " ".join(
        "%.2f" % x for x in starting_point))

    # fit the model parameters
    idr.log("Fitting the model parameters", 'VERBOSE');

    theta, loss = estimate_model_params(
        r1, r2,
        starting_point,
        max_iter=max_iter,
        convergence_eps=convergence_eps,
        fix_mu=fix_mu, fix_sigma=fix_sigma)

    idr.log("Finished running IDR on the datasets", 'VERBOSE')
    idr.log("Final parameter values: [%s]" % " ".join("%.2f" % x for x in theta))

    # calculate the global IDR
    localIDRs = calc_local_IDR(np.array(theta), r1, r2)
    return localIDRs


def calc_global_IDR(localIDR):
    local_idr_order = localIDR.argsort()
    ordered_local_idr = localIDR[local_idr_order]
    ordered_local_idr_ranks = rankdata( ordered_local_idr, method='max' )
    IDR = []
    for i, rank in enumerate(ordered_local_idr_ranks):
        IDR.append(ordered_local_idr[:rank].mean())
    IDR = np.array(IDR)[local_idr_order.argsort()]
    return IDR


def idr_plot(samples, plot_file):
    with open(samples[0], "r") as fp1:
        with open(samples[1], "r") as fp2:
            merged_peaks, signal_type = load_idr_samples([fp1, fp2])

    s1 = np.array([pk.signals[0] for pk in merged_peaks])
    s2 = np.array([pk.signals[1] for pk in merged_peaks])

    # build the ranks vector
    idr.log("Ranking peaks", 'VERBOSE')
    r1, r2 = build_rank_vectors(merged_peaks)

    localIDRs = fit_model_and_calc_local_idr(
        r1, r2,
        starting_point=(
            idr.DEFAULT_MU, idr.DEFAULT_SIGMA,
            idr.DEFAULT_RHO, idr.DEFAULT_MIX_PARAM),
        max_iter=idr.MAX_ITER_DEFAULT,
        convergence_eps=idr.CONVERGENCE_EPS_DEFAULT,
        fix_mu=False, fix_sigma=False)

    IDRs = calc_global_IDR(localIDRs)
    plot_idr([s1, s2], [r1, r2], IDRs, plot_file, )

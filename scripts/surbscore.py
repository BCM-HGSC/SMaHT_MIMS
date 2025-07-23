"""
Automatic scoring of stratifications
"""
import os
import sys
import logging
import argparse
import joblib
import truvari
from truvari.vcf2df import vcf2df_main
import numpy as np
import pandas as pd
# requires scipy==1.10.1 and statsmodels==0.13.5
from scipy.stats import chi2_contingency
from statsmodels.stats.multitest import multipletests


def get_scores(df, state="state", features=["svtype", "szbin"], min_obs=10):
    """
    Actively filter 0 in either - has to be something
    """
    tot_counts = []
    for feat in features:
        cnt = df[feat].value_counts()
        cnt.name = "count"
        cnt = pd.DataFrame(cnt)
        cnt['feat'] = feat
        tot_counts.append(cnt)
    tot_counts = pd.concat(tot_counts)
    tot_counts = tot_counts.reset_index().set_index(['feat', 'index'])

    logging.debug("Full Counts\n%s", str(tot_counts))
    all_cont_obs = pd.crosstab(df[state], [df[i] for i in features])
    logging.info("%d stratifications had ≥1 observed variant",
                 all_cont_obs.shape[1])
    view = pd.DataFrame(all_cont_obs.sum(axis=0))
    filtered = view[0] < min_obs
    if filtered.any():
        logging.warning(
            "%d stratifications with < %d observations ignored", filtered.sum(), min_obs)
        with pd.option_context('display.max_rows', 100):
            logging.debug("Fitered Bins\n%s", str(all_cont_obs.T[filtered]))

    cont_obs = all_cont_obs.T[~filtered].T

    chi, pval, dof, exp = chi2_contingency(cont_obs)
    score = (np.sign(cont_obs.values[1] - exp[1])) * \
            (cont_obs.values[1] - exp[1])**2 / exp[1]

    df_rank = pd.DataFrame(list(zip(list(cont_obs.columns),
                                    list(cont_obs.values[1]), list(exp[1]), list(score))),
                           columns=["Values", "observed", "expected", "score"])

    df_rank.sort_values(["score"], ascending=True, inplace=True)
    df_rank.reset_index(inplace=True, drop=True)
    df_rank.index.names = ['rank']
    df_rank[features] = pd.DataFrame(
        df_rank["Values"].tolist(), index=df_rank.index)
    df_rank = df_rank.drop("Values", axis=1)
    return df_rank, tot_counts, all_cont_obs


def permutation_test(a_values, b_values, n_samps=10000, tailed="two"):
    """
    a_values - b_values mean difference from permutations' mean differences
    tail can be "two", "left", or "right"
    """
    delta = a_values.mean() - b_values.mean()
    if tailed == "two":
        delta = abs(delta)
    logging.debug("Case %.2f (%d)", a_values.mean(), len(a_values))
    logging.debug("Control %.2f (%d)", b_values.mean(), len(b_values))
    logging.debug("delta %.2f", delta)

    def run_t(pooled, sizeA, sizeB):
        np.random.shuffle(pooled)
        starA = pooled[:sizeA]
        starB = pooled[-sizeB:]
        return starA.mean() - starB.mean()
    pooled = np.hstack([a_values, b_values])
    estimates = np.array(
        [run_t(pooled, a_values.size, b_values.size) for _ in range(n_samps)])
    if tailed == "left":
        diffCount = len(np.where(estimates <= delta)[0])
    elif tailed == "right":
        diffCount = len(np.where(estimates >= delta)[0])
    else:
        estimates = np.abs(estimates)
        diffCount = len(np.where(estimates >= delta)[0])
    pct_extreme = float(diffCount) / float(n_samps)
    one, fif, nine = np.quantile(estimates, [0.01, 0.5, 0.99])
    logging.debug(
        "Quantiles: 1% {:.10f} - 50% {:.10f} - 99% {:.10f}".format(one, fif, nine))
    logging.debug("Pct at least as extreme %.2f", pct_extreme)
    return pct_extreme, delta, one, fif, nine


def surb_score(m_data, state="state", features=["svtype", "szbin"],  min_obs=10, nperm=10000, tail="left", alpha=0.01):
    """
    create the surb score data-frame that's plottable
    """
    rank, counts, all_counts = get_scores(m_data, state, features, min_obs=10)
    rows = []
    # Also need the pvalue distributions (normalized?)
    for feat in features:
        for val in rank[feat].unique():
            sub = rank[feat] == val
            case = np.array(rank[sub].index)
            control = np.array(rank[~sub].index)
            if len(case) < min_obs or len(control) < min_obs:
                logging.warning("%s : %s with < %d observations ignored (%d case, %d control)", feat, val, min_obs, len(
                    case), len(control))
                continue
            logging.debug("Testing %s : %s", feat, val)
            pval, delta, one, fif, nin = permutation_test(
                case, control, tailed=tail)
            acc = m_data[m_data[feat] == val][state].mean()
            rows.append([feat, val, acc, delta, one, fif, nin,
                        counts.loc[feat, val]["count"], pval])

    plt_rank = pd.DataFrame(rows, columns=[
                            "feature", "value", "acc", "delta", "q1", "q50", "q99", "obs", "pval"])
    reject, adj_pval, x, y = multipletests(
        pvals=plt_rank['pval'], alpha=alpha, method="fdr_bh")
    adj = pd.concat([plt_rank,
                     pd.Series(adj_pval, name="adj_pval"),
                     pd.Series(reject, name="reject"),
                     ], axis=1)
    return adj, all_counts

def bin_equal_freq(series, n_bins) -> pd.Series:
    """
    Bin a continuous series into equal-frequency (quantile) bins.
    
    Parameters:
    - series: The continuous variable as a pandas Series.
    - n_bins: The number of bins to create.
    
    Returns:
    - A Series with bin labels (0 to n_bins-1).
    """
    return pd.qcut(series, q=n_bins, duplicates='drop')

def bin_equal_width(series, n_bins) -> pd.Series:
    """
    Bin a continuous series into equal-width bins.
    
    Parameters:
    - series: The continuous variable as a pandas Series.
    - n_bins: The number of bins to create.
    
    Returns:
    - A Series with bin labels (0 to n_bins-1).
    """
    return pd.cut(series, bins=n_bins, include_lowest=True)

# PresetName: (FeatureList, StatesStr)
known_features = {'mims': (['svtype', 'szbin', 'TRF', 'isolated', 'VAF_bin'], 'base')}


def parse_args(args):
    parser = argparse.ArgumentParser(prog="surbscore", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("bench_dir", metavar="DIR", type=str,
                        help="Truvari bench result directory")
    parser.add_argument("--features", type=str, default=None,
                        help="Comma-separated list of categorical features to score (e.g. svtype,szbin)")
    parser.add_argument("--states", type=str, default='base', choices=['base', 'comp', 'both'],
                        help="Score VCF entries from (base)line, (comp)arison, or (both) (default:base)")
    parser.add_argument("--preset", type=str, default=None, choices=known_features.keys(),
                        help="Preset features to score")
    parser.add_argument("-o", "--output", type=str, default="/dev/stdout",
                        help="Score output tsv file (stdout)")
    parser.add_argument("-c", "--count-output", type=str, default=None,
                        help="Stratification variant count tsv file")
    parser.add_argument("--debug", action='store_true',
                        help="Verbose logging")

    thresg = parser.add_argument_group("statistics arguments")
    thresg.add_argument("--min-obs", type=int, default=10,
                        help="Minimum observation for a bin (%(default)s)")
    thresg.add_argument("--nperm", type=int, default=10000,
                        help="Number of permutations to perform (%(default)s)")
    thresg.add_argument("--tail", type=str, default="left", choices=['left', 'right', 'two'],
                        help="Tail to test (%(default)s)")
    thresg.add_argument("--alpha", type=truvari.restricted_float, default=0.01,
                        help="Significance threshold (%(default)s)")
    thresg.add_argument("--seed", type=int, default=42,
                        help="Random seed for permutation tests (%(default)s)")

    args = parser.parse_args(args)
    truvari.setup_logging(args.debug, show_version=False)
    if not args.features:
        if not args.preset:
            logging.error("Either `--features` or `--preset` must be provided")
            sys.exit(1)
        args.features, args.states = known_features[args.preset]
    else:
        args.features = args.features.split(',')
    return args


def set_mims_bins(data):
    """
    Categorization of DataFrame columns based when they contain SMaHT MIMS annotations
    """
    data['is_tp'] = data['state'].str.startswith('tp')
    data['isolated'] = data['NumNeighbors'] == 0
    bins = [0, 0.05, 0.3, 0.7, 1]
    labels = ["SomaticLow(<5%)", "Somatic", "GermlineHet", "GermlineHom"]
    data['VAF_bin'] = pd.cut(data['VAF_alt'], bins=bins, labels=labels)


if __name__ == '__main__':
    """
    Run SurbScore on a truvari directory with benchmarking results
    """
    args = parse_args(sys.argv[1:])
    np.random.seed(args.seed)
    data_file = os.path.join(args.bench_dir, "data.jl")
    if not os.path.exists(data_file):
        logging.info("Constructing DataFrame in %s", data_file)
        vcf2df_main(f"-i -f -b {args.bench_dir} {data_file}".split(' '))

    data = joblib.load(data_file)
    if args.preset == "mims":
        set_mims_bins(data)
    else:
        data['is_tp'] = data['state'].str.startswith('tp')

    if args.states == 'base':
        data = data[data['state'].isin(['tpbase', 'fn'])]
    elif args.states == "comp":
        data = data[data['state'].isin(['tp', 'fp'])]

    for i in range(len(args.features)):
        col = args.features[i]
        if ':' in col:
            name, key = col.split(':')
            method, cnt = key[0], key[1:]
            cnt = int(cnt)
            if method == 'f':
                data[name] = bin_equal_freq(data[name], cnt)
            elif method == 'f':
                data[name] = bin_equal_width(data[name], cnt)
            args.features[i] = name

    n_strats = 1
    n_values = 0
    for i in args.features:
        c = data[i].nunique()
        n_strats *= c
        n_values += c
    n_feat = len(args.features)

    logging.info("Scoring %d features with %d values (max %d stratifications)",
                 n_feat, n_values, n_strats)


    scores, counts = surb_score(data, state="is_tp",
                                features=args.features,
                                min_obs=args.min_obs,
                                nperm=args.nperm,
                                tail=args.tail,
                                alpha=args.alpha,
                                )

    if args.count_output:
        counts.T.to_csv(args.count_output, sep='\t')
    scores.round(4).sort_values(by=['adj_pval', 'pval', 'acc'], ascending=True).to_csv(
        args.output, sep='\t', index=False)
    logging.info("Finished")

import os
import sys
import pysam
import argparse
import joblib
import truvari
import numpy as np
import pandas as pd
import seaborn as sb
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
from truvari.vcf2df import bench_dir_to_df

def size_type_metrics(df):
    """
    Makes a dataframe of precision/recall/f1 by sv type and size
    """
    rows = []
    for i, sub in df.groupby(['szbin', 'svtype'], observed=False):
        cnts = sub['state'].value_counts()
        rows.append((*i, * truvari.performance_metrics(**cnts)))
    ret = pd.DataFrame(rows, columns=['szbin', 'svtype', 'precision', 'recall', 'f1'])
    return ret.dropna()

def size_type_barplot(metrics, y="recall", title="Size/Type Performance"):
    """
    Simple barplot of y=recall/precision/f1 by svsize hue svtype
    metrics = output from `size_type_metrics`
    Returns a seaborn plot
    """
    p = sb.barplot(data=metrics, x='szbin', y=y, hue='svtype')
    p.set(title=title, xlabel="SizeBin")
    _ = plt.xticks(rotation=45, ha='right')
    return p

def split_base_comp(data):
    """
    Separates dataframe rows into those from the baseline and from the comparison VCF.
    Extracts Truvari's MatchId annotation so the subsets have a key to potentially joine one another
    """
    # Separate out the variants from the base VCF and add new columns of the base/comp ids
    base = data[data['state'].isin(['tpbase', 'fn'])].copy()
    base['base_id'] = base['MatchId'].apply(lambda x: x[0])
    base['comp_id'] = base['MatchId'].apply(lambda x: x[1])

    # Separate out the variants from the comparison VCF and add new columns of the base/comp ids
    comp = data[data['state'].isin(['tp', 'fp'])].copy()
    comp['base_id'] = comp['MatchId'].apply(lambda x: x[0])
    comp['comp_id'] = comp['MatchId'].apply(lambda x: x[1])
    base.reset_index(inplace=True)
    comp.reset_index(inplace=True)
    return base, comp

def severus_vaf(df, sample):
    """
    Severus already calculates it
    """
    return df[f'{sample}_VAF']

def sniffles_vaf(df, sample):
    """
    Calculate sniffles' predicted vaf
    """
    dv = df[f'{sample}_DV']
    dr = df[f'{sample}_DR']
    tot = dv + dr
    return dv / tot

def sawfish_vaf(df, sample):
    """
    Calculate sawfish' predicted vaf
    """
    dr = df[f'{sample}_AD_ref']
    dv = df[f'{sample}_AD_alt']
    tot = dr + dv
    return dv / (dv+dr)

def svdss_vaf(df, sample):
    """
    Calculate sawfish' predicted vaf
    """
    return df['WEIGHT'] / df['COV']

def pbsv_vaf(df, sample):
    return sawfish_vaf(df, sample)

def delly_vaf(df, sample):
    return sniffles_vaf(df, sample)

def svaba_vaf(df, sample):
    return df[f'{sample}_AD'] / df[f'{sample}_DP']

def gridss_vaf(df, sample):
    vf = df[f'{sample}_VF']
    ref = df[f'{sample}_REF']
    refpair = df[f'{sample}_REFPAIR']
    tot = vf + ref
    big_tot = tot + refpair
    vaf = vf / tot.where(df['svlen'] < 1000, big_tot)
    return vaf

tool_vaf = {'sniffles': sniffles_vaf,
            'svdss': svdss_vaf,
            'sawfish': sawfish_vaf,
            'severus': severus_vaf,
            'pbsv': pbsv_vaf,
            'delly': delly_vaf,
            'svaba': svaba_vaf,
            'gridss': gridss_vaf,
            }

def get_predicted_vaf(comp, sample="", tool=None):
    """
    Calculate the predicted_vaf annotation for a dataframe 
    Supported tools: svdss, sniffles, sawfish
    """
    if tool not in tool_vaf.keys():
        raise Exception(f"Unexpected tool name {tool}")
    return tool_vaf[tool](comp, sample)

def get_vaf_bin(vaf, bins=[0, 0.05, 0.3, 0.7, 1], labels=["Ultra-Low", "Low", "Medium", "High"]):
    """
    Bin VAFs and apply labels
    """
    return pd.cut(vaf, bins=bins, labels=labels) 

def fill_in_vafs(base, comp):
    """
    Sets predicted and expected vaf columns for base and comp dataframes
    """
    sum_vaf = base[base['state'] == 'tpbase'].groupby(['comp_id'])['expected_vaf'].sum()
    a = comp[(comp['state'] == 'tp') & ~comp['comp_id'].isin(sum_vaf.index)][['base_id', 'comp_id']].copy()
    b = base[base['base_id'].isin(a['base_id'])]
    a.set_index('base_id', inplace=True)
    b.set_index('base_id', inplace=True)
    a['expected_vaf'] = b['expected_vaf']
    a = a.reset_index(drop=True).set_index('comp_id')
    sum_vaf = pd.concat([sum_vaf, a])

    comp.set_index('comp_id', inplace=True)
    comp['expected_vaf'] = sum_vaf.where(sum_vaf <= 1, 1)
    comp.reset_index(inplace=True)

    sum_vaf = comp[comp['state'] == 'tp'].groupby(['base_id'])['predicted_vaf'].sum()
    a = base[(base['state'] == 'tpbase') & ~base['base_id'].isin(sum_vaf.index)][['base_id', 'comp_id']].copy()
    b = comp[comp['comp_id'].isin(a['comp_id'])]
    a.set_index('comp_id', inplace=True)
    b.set_index('comp_id', inplace=True)

    a['predicted_vaf'] = b['predicted_vaf']
    a = a.reset_index(drop=True).set_index('base_id')
    sum_vaf = pd.concat([sum_vaf, a])

    base.set_index('base_id', inplace=True)
    base['predicted_vaf'] = sum_vaf.where(sum_vaf <= 1, 1)
    base.reset_index(inplace=True)

def vaf_performance(predicted, expected):
    """
    Summarize the difference of predicted and expected vaf
    Returns a DataFrame with the description of the predicted-expected distribution
    As well as PearsonR coefficient and Pval, and mean-square error (MSE)
    """
    o1 = (predicted - expected).describe().to_frame()
    o2 = pearsonr(predicted, expected)
    o1.loc["PearsonR"] = o2.statistic
    o1.loc["Pval"] = o2.pvalue

    MSE = ((predicted - expected) ** 2).sum()  / len(predicted)
    o1.loc['MSE'] = MSE
    return o1

def comp_add_base_anno(comp, base):
    """
    Add the NumNeighbors and TRF annotations from base to the comp dataframe
    Edits them in place
    """
    comp.set_index('base_id', inplace=True)
    base.set_index('base_id', inplace=True)

    comp['NumNeighbors'] = base['NumNeighbors']
    comp['TRF'] = base['TRF']
    comp.reset_index(inplace=True)
    base.reset_index(inplace=True)

def make_mosaic_df(base, comp):
    """
    Given a truvari vcf2df dataframe
    Unite the predicted/expected vaf and baseline specific annotations 
    Returns a new, cleaned up datafraome
    """
    base['expected_vaf'] = base['VAF_alt']

    fill_in_vafs(base, comp)
    comp_add_base_anno(comp, base)

    columns = ['svtype', 'svlen', 'szbin', 'TRF', 'NumNeighbors', 'NeighId', 
               'expected_vaf', 'predicted_vaf', 'state', 'hash', 'base_id', 'comp_id']

    base = base[columns]
    comp = comp[columns]

    data = pd.concat([base, comp])

    data['expected_vaf_bin'] = get_vaf_bin(data['expected_vaf'])
    data['predicted_vaf_bin'] = get_vaf_bin(data['predicted_vaf'])
    data['matching_vaf_bin'] = data['expected_vaf_bin'] == data['predicted_vaf_bin']
    data['vaf-delta'] = data['predicted_vaf'] - data['expected_vaf']
    return data

def collapse_mosaic(data):
    """
    Returns a 'collapsed' set of variants
    """
    # Filter base and comp states once
    base = data[data['state'].isin(['tpbase', 'fn'])]
    comp = data[data['state'].isin(['tp', 'fp'])].copy()
    comp['comp_id'] = comp['comp_id'].astype(float)

    # Process tpbase entries
    tpbase = base[base['state'] == 'tpbase']
    grouped = tpbase.groupby('comp_id')

    # Identify the best matches and collect kept and collapsed variants
    best_match_ids = comp.set_index('comp_id')['base_id']
    kept = []
    collapsed = []

    for cid, subset in grouped:
        if len(subset) == 1:
            kept.append(subset)
        else:
            best_ids = subset['base_id'].isin(best_match_ids.get(cid, []))
            kept.append(subset[best_ids])
            collapsed.append(subset[~best_ids])

    # Combine kept variants
    m_kept = pd.concat(kept, ignore_index=True)

    # Add back 'fn' entries
    new_base = pd.concat([m_kept, base[base['state'] == 'fn']], ignore_index=True)

    # Handle missing base_ids in new_base
    missing_base_ids = comp.loc[~comp['base_id'].isin(new_base['base_id']), 'base_id']
    missing_base = base.loc[base['base_id'].isin(missing_base_ids)]
    new_base = pd.concat([new_base, missing_base], ignore_index=True)

    # Return combined result
    return pd.concat([new_base, comp], ignore_index=True)

def collapse_mosaic_old(data):
    """
    Returns a 'collapsed' set of variants
    Collapsing is performed by saying: if multiple baseline variants match to a single comparison variant,
    we say all these baseline variants are similar enough that they look like a single variant to the caller
    We then keep the single baseline variant that best matches to the comp, and remove the rest 
    But then we put back in any missing base which are a comp is missing
    """
    base = data[data['state'].isin(['tpbase', 'fn'])]
    comp = data[data['state'].isin(['tp', 'fp'])]
    kept = []
    collapsed = []
    for cid, subset in base[base['state'] == 'tpbase'].groupby("comp_id"):
        if len(subset) == 1:
            kept.append(subset)
            continue
        other = comp[comp['comp_id'].astype(float) == cid]
        best = other['base_id']
        k = subset['base_id'].isin(best)
        kept.append(subset[k])
        collapsed.append(subset[~k])
    m_kept = pd.concat(kept)
    # Need to put back in the base calls which are pointed at by comp
    new_base = pd.concat([m_kept, base[base['state'] == 'fn']])
    missing_base_ids = comp[~comp['base_id'].isin(new_base['base_id'])]['base_id']
    missing_base = base[base['base_id'].isin(missing_base_ids)]
    new_base = pd.concat([new_base, missing_base])
    return pd.concat([new_base, comp])

def parse_args(args):
    """
    """
    parser = argparse.ArgumentParser(prog="mims_summary", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--input", type=str, required=True,
                        help="Truvari bench result directory")
    parser.add_argument("-p", "--program", type=str, required=True, choices=tool_vaf.keys(),
                        help="Program which produced the comparison VCF")
    parser.add_argument("-o", "--output", type=str, default="mims.summary.txt",
                        help="Output tsv")
    args = parser.parse_args(args)
    return args
    # write a tsv of the summarized table

if __name__ == '__main__':
    args = parse_args(sys.argv[1:])
    tru_data = bench_dir_to_df(args.input, with_info=True, with_format=True)
    base, comp = split_base_comp(tru_data)
    
    vcf = pysam.VariantFile(os.path.join(args.input, "tp-comp.vcf.gz"))
    sample = vcf.header.samples[0]
    comp['predicted_vaf'] = get_predicted_vaf(comp, sample, args.program)

    # Now, we can recombine and simplify the annotations
    non_dedup = make_mosaic_df(base, comp)
    data = collapse_mosaic(non_dedup)
    data.to_csv(args.output, sep='\t', index=False)

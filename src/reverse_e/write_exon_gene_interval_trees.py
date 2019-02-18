import pandas as pd
import pickle
from intervaltree import IntervalTree


def rebuild_gene_interval_tree():
    ts = {}
    # I just realized that dynamically building this results in an overhead
    # TODO: build the Interval tree statically
    for index, row in genes.iterrows():
        l = row["start"]
        r = row["end"]
        g = row["chromosome"]
        if g not in ts:
            ts[g] = IntervalTree()
        ts[g][l:(r + 1)] = ["gene", index]
        if row['strand'] == '+':
            ts[g][(l - FLANK_RANGE):l] = ["promoter", index]
        else:
            ts[g][(r + 1):(r + FLANK_RANGE + 1)] = ["promoter", index]
    return ts


def rebuild_exon_interval_tree():
    # TODO: build the Interval tree statically
    name2id = {}
    for index, row in genes.iterrows():
        name2id[row["name"]] = index
    exon_ts = {}
    for index, row in exons.iterrows():
        l = row["start"]
        r = row["end"]
        gid = row["Gene"]
        if gid not in exon_ts:
            exon_ts[gid] = IntervalTree()
        exon_ts[gid][l:(r + 1)] = name2id[gid]
    return exon_ts


def rebuild_all_intervals_trees():
    # rebuild inverval trees and serialize them
    pickle.dump(rebuild_gene_interval_tree(), open("ts.p", "wb"))
    pickle.dump(rebuild_exon_interval_tree(), open("exon_ts.p", "wb"))


GENE_PATH = "GeneData.csv"
EXON_PATH = "whitelist_exon.csv"
FLANK_RANGE = 1000
exons = pd.read_csv(EXON_PATH)
genes = pd.read_csv(GENE_PATH)
rebuild_all_intervals_trees

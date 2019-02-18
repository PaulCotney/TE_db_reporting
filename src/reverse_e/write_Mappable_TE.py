import pandas as pd
import pickle
from intervaltree import IntervalTree


# helper function for retrieving gene location information
def lookup_(i, loc):
    t = ts[i]
    retrieved = t[loc]
    retrieved = list(map(lambda x: x[2], retrieved))
    r = []
    for e in retrieved:
        if e[0] == 'gene' or e[0] == 'intron':
            e[0] = 'intron'
            g = genes.iloc[e[1]]
            exons_hit = []
            if g["name"] in exon_ts:
                exons_hit = map(lambda x: x[2], exon_ts[g["name"]][loc])
            if len(exons_hit) > 0:
                for ex in exons_hit:
                    r.append(["exon", ex])
            else:
                r.append(e)
        else:
            r.append(e)
    return r


def lookup(i, loc):
    if isinstance(i, int):
        i = str(i)
    r = lookup_(i, loc)
    return r

    # modified = collapsed.copy()


def mapping(v, d):
    g = genes.iloc[v[1]]
    return [v[0], g['name'], g['ensemblID'], g['strand'] == d['strand']]


def joiner(pn):
    return '||'.join(['|'.join(p) for p in pn])


def append_region_info(tbl):
    for i, d in tbl.iterrows():
        chromo = d['chromo']
        pos = d['pos']
        r = lookup(chromo, pos)
        r_ = map(lambda x: mapping(x, d), r)
        r_.sort(key=lambda x: x[0])
        gns = genes.iloc[list(map(lambda x: x[1], r))]
        if len(r) == 0:
            tbl.loc[i, 'region'] = 'intergenic'
            continue
        pos = filter(lambda x: x[3] == True, r_)
        neg = filter(lambda x: x[3] == False, r_)
        p = [pos, neg]
        if len(list(neg)) == 0:
            p = [pos]
        tbl.loc[i, 'region'] = joiner([[ps[0] for ps in pn] for pn in p])
        tbl.loc[i, 'gene_id'] = joiner([[ps[1] for ps in pn] for pn in p])
        tbl.loc[i, 'gene_name'] = joiner([[ps[2] for ps in pn] for pn in p])


def append_orientation_info(tbl):
    for i, d in tbl.iterrows():
        chromo = d['chromo'] # get TE's chromosome
        pos = d['pos'] # get TE's position
        r = lookup(chromo, pos) # lookup the chromosome and position, to the result
        if len(r) == 0:
            continue
        if r[0][0] in ['gene', 'intron']:
            if len(r) > 1:
                tbl.loc[i, 'same_orientation'] = -1
            else:
                loc = r[0][1]
                gene = genes.iloc[loc]
                b = gene['strand'] == d['strand']
                tbl.loc[i, 'same_orientation'] = int(b)


# OUT_PATH will be set to be in_path if empty
def append_region_info_to_file(in_path, out_path=None):
    out_path = out_path or in_path
    tbl = pd.read_csv(in_path, dtype = {'chromo':str, 'pos':int})
    append_region_info(tbl)
    tbl.to_csv(out_path, index=False)


# Loading datasets
GENE_PATH="GeneData.csv"
genes = pd.read_csv(GENE_PATH)

IN_PATH="FinalMatrix_ALL_collapsed.csv"
OUT_PATH = "Mappable_TEs.csv"


#Two dictionaries
ts = pickle.load( open( "ts.p", "rb" ) )
exon_ts = pickle.load( open( "exon_ts.p", "rb" ) )
append_region_info_to_file(IN_PATH, OUT_PATH)
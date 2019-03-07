from datetime import datetime
import csv
import gc
import glob
import numpy as np
import os
import pandas as pd
import sqlite3
import time

def getPopulation():
    founder_bwt = '/csbiodata/BWTs/sangerMsbwt'
    cc_bwt = '/csbiodata/perlegen/CC_bwt'
    cc_avoid = ['CC039M3730_UNC_NYGC', 'CC019M6839_UNC_UNC', 'CC081M0332_UNC_UNC']
    sanger_founder = ["C57BL6JN", "129S1", "NOD", "NZOdna", "CASTdna", "PWKdna", "WSBdna"]
    unc_founder = ["AJ_F321", "C57BL6J_m003636A", "129S1_M157", "NOD_M146", "NZO_FNNN", "CAST_M090", "PWK_F6114",
                   "WSB_FNNN"]
    cc_samples = []
    bwt_list = sorted(glob.glob("%s/*" % cc_bwt))
    for bwtfile in bwt_list:
        sample = bwtfile.split("/")[-1]
        if sample in cc_avoid:
            continue
        cc_samples.append(sample)

    sample_list = []
    sample_population = {}
    for sample in sanger_founder:
        sample_population[sample] = 'sanger_founder'
    for sample in unc_founder:
        sample_population[sample] = 'unc_founder'
    for sample in cc_samples:
        sample_population[sample] = 'cc'

    return sample_population


def ref_seq(side):
    filename = 'tmp/bowtie_data/input_file_%s_final.csv' % (side)
    df = pd.read_csv(filepath_or_buffer=filename, sep=',', dtype={'chromo': str, 'pos': int})
    data = df.loc[:, ].values
    ref_dict = {}
    for d in data:
        my_id, context, te, ref_te, ref_prefix, ref_suffix, chromo, pos, strand = d[0:9]
        ref_prefix = ref_prefix[-25:]
        ref_suffix = ref_suffix[:25]
        ref_dict[(chromo, pos, strand)] = [ref_prefix, ref_suffix]

    return ref_dict


def get_ref_te():
    filename = "tmp/Table/MERGED_v1.csv"
    df = pd.read_csv(filepath_or_buffer=filename, dtype= {'chromo': str,'pos': int}, sep=',')
    df = df.replace(np.nan, '', regex=True)
    df = df[(df['TE'] == 1) & (df['REF'] == 1)]
    data = df.iloc[:,].values
    ref_te = set()
    for d in data:
        [my_id,side,chromo,pos,strand] = d[0:5]
        ref_te.add((chromo,pos,strand))
    return ref_te

def get_counts_te():
    filename = "MERGED_v1.csv"
    dtype = {'chromo': str, 'pos': int}
    df = pd.read_csv(filepath_or_buffer=filename, dtype=dtype, sep=',')
    df = df.replace(np.nan, '', regex=True)
    df = df[(df['TE'] == 1)]
    data = df.iloc[:, ].values
    ref_dict = ref_seq('start')
    ref_te = get_ref_te()
    header = list(df)
    sample_list = header[11:]
    pos_dict = {}
    all_te_context = {}

    pos_to_TE_id = {}
    for d in data:
        [my_id, side, chromo, pos, strand] = d[0:5]
        pos_to_TE_id[(chromo, pos, strand)] = my_id
        [r, t] = d[9:11]
        [p, ts, te, s] = d[5:9]
        mystring = p + ts + te + s
        # other_context = True if len([i for i in mystring if i.islower()]) > 0 else False

        pid = pos_dict.get((chromo, pos, strand), len(pos_dict))
        pos_dict[(chromo, pos, strand)] = pid

        if (chromo, pos, strand) in ref_te:
            ref_prefix, ref_suffix = '', ''
            [ref_prefix, _] = ref_dict[(chromo, pos, strand)]
            all_te_context[(chromo, pos, strand)] = [ref_prefix, ref_suffix]
        else:
            [ref_prefix, ref_suffix] = ref_dict[(chromo, pos, strand)]
            all_te_context[(chromo, pos, strand)] = [ref_prefix, ref_suffix]

    m, n = len(pos_dict), len(sample_list)
    with_te = np.zeros((m, n), dtype=int)
    without_te = np.zeros((m, n), dtype=int)

    for d in data:
        [my_id, side, chromo, pos, strand] = d[0:5]
        [r, t] = d[9:11]
        pid = pos_dict[(chromo, pos, strand)]
        sample_count = d[11:]
        for sid, c in enumerate(sample_count):
            with_te[pid][sid] += c

    labels = []
    pid_to_chromo = {v: k for k, v in pos_dict.iteritems()}
    for pid in pid_to_chromo:
        [chromo, pos, strand] = pid_to_chromo[pid]
        l = [chromo, pos, strand]
        labels.append(l)

    filename = "tmp/Table/MERGED_v1.csv"
    dtype = {'chromo': str, 'pos': int, 'old_chromo': str, 'old_pos': int}
    df = pd.read_csv(filepath_or_buffer=filename, dtype=dtype, sep=',')
    df = df.replace(np.nan, '', regex=True)
    df = df[(df['TE'] == 0)]
    data = df.iloc[:, ].values

    for d in data:
        [my_id, side, chromo, pos, strand] = d[0:5]
        [r, t] = d[9:11]
        if (chromo, pos, strand) not in pos_dict.keys():
            continue
        pid = pos_dict[(chromo, pos, strand)]
        sample_count = d[11:]
        for sid, c in enumerate(sample_count):
            without_te[pid][sid] += c

    return with_te, without_te, labels, sample_list, all_te_context, pos_to_TE_id


def make_string(a):
    a = "'%s'" % a
    return a


def get_genotype(high_count_strain, chromo, pos, state):
    DATABASE = '/home/mcmillan/Notebooks/Anwica/CCIntervalDatabase/CCIntervalsAK.db'
    db = sqlite3.connect(DATABASE)
    db.text_factory = str
    cursor = db.cursor()

    pos_start, pos_end = pos, pos
    chromo = make_string(chromo)
    query = """select S.name,I.genotype from Intervals I, Samples S where S.sid=I.sid 
                and I.chromosome = %s and I.start <= %d and I.end >= %d""" % (chromo, pos_start, pos_end)
    cursor.execute(query)
    slist = [i for i in cursor.fetchall() if i[0].find('CC081') < 0]

    founders = ['AA', 'BB', 'CC', 'DD', 'EE', 'FF', 'GG', 'HH']
    genotype_to_strain, strain_to_genotype = {}, {}
    for s, g in slist:
        strain = s[:5]
        if g not in founders:
            continue
        genotype_to_strain[g] = genotype_to_strain.get(g, []) + [strain]
        strain_to_genotype[strain] = g
    db.close()

    geno_with_TE = set()
    for strain in high_count_strain:
        g = strain_to_genotype.get(strain, 'XX')
        geno_with_TE.add(g)

    fully_consistent_geno = []
    if state == 'private':
        return "%s(%s)" % (list(high_count_strain)[0], list(geno_with_TE)[0])

    geno_with_TE = geno_with_TE - set(['XX'])
    for g in geno_with_TE:
        strain_list = genotype_to_strain[g]
        clist = set(strain_list) - set(high_count_strain)
        if len(clist) == 0:
            fully_consistent_geno.append(g)
        else:
            consistent = len(strain_list) - len(clist)
            total = len(strain_list)
            g = "%s(%d/%d)" % (g, consistent, total)
            fully_consistent_geno.append(g)

    return "|".join(sorted(fully_consistent_geno))

def collapse():
    # with_te, without_te, labels, sample_list, all_te_context, pos_to_TE_id
    with_te, without_te, labels, sample_list, all_te_context, pos_to_TE_id = get_counts_te()
    (m, n) = with_te.shape

    collapsed_count = np.zeros((m, n), dtype=int) - 1
    sample_population = getPopulation()
    new_data = []
    # m = 10
    ### change above line later
    for pid in range(m):
        [chromo, pos, strand] = labels[pid]
        my_id = pos_to_TE_id[(chromo, pos, strand)]
        side = 'start'
        [ref_prefix, ref_suffix] = all_te_context[(chromo, pos, strand)]
        sample_count = with_te[pid]

        # print sample_count
        high_count = np.where(with_te[pid] > 4)[0]
        high_count_sample = [sample_list[sid] for sid in high_count]

        high_count_strain = set()
        non_cc = False
        for sample in high_count_sample:
            if sample_population[sample] == "cc":
                strain = sample[:5]
            else:
                strain = sample
                non_cc = True
            high_count_strain.add(strain)

        state = ''
        if len(high_count_strain) == 1:
            state = 'questionable' if non_cc else 'private'
        elif len(high_count_sample) == len(sample_list):
            state = 'fixed'
        else:
            state = 'shared'
        geno = get_genotype(high_count_strain, chromo, pos, state)

        for sid in high_count:
            collapsed_count[pid][sid] = 1

        high_count_nonTE = np.where(without_te[pid] > 4)[0]

        for sid in high_count_nonTE:
            collapsed_count[pid][sid] += 1

        row_data = [my_id, side, chromo, pos, strand, state, ref_prefix, ref_suffix, geno] + list(collapsed_count[pid])
        new_data.append(row_data)

    print(collapsed_count)

    header = ['my_id', 'side', 'chromo', 'pos', 'strand', 'state', 'before', 'after', 'genotype']
    header = header + sample_list

    filename = "FinalMatrix_ALL_collapsed.csv"
    with open(filename, 'wb') as fp:
        a = csv.writer(fp, delimiter=',')
        a.writerows([header])
        for row_data in new_data:
            a.writerows([row_data])
            # print row_data
    print("Wrote file %s [%d lines]" % (filename, len(new_data)))


collapse()

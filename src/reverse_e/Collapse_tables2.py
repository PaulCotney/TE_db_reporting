import csv
import numpy as np
import pandas as pd
import sqlite3
import datetime


# Copied from ./Anwica/TransposableElement/TEprobeCounts/Step03_CollapseTable-allState.ipynb
database_file = '../TE_db.sqlite'
connection = sqlite3.connect(database_file)
connection.text_factory = str
source_df = None


# This is stupidly inefficient, but is temporary
def build_source_df():
    cursor = connection.cursor()
    sample_name_loc = dict()
    loc_sample_name = dict()
    pos = 0
    for sample_name in cursor.execute('select distinct name from sample'):
        sample_name_loc[sample_name[0]] = pos
        loc_sample_name[pos] = sample_name[0]
        pos = pos + 1
    row_counts = cursor.execute('select count(*) from tmp_summary').fetchone()[0]
    column_counts = len(sample_name_loc.keys())
    calls = np.zeros((row_counts, column_counts), dtype=int)
    summary = []
    row_pos = {}
    row_count = 0
    for row in cursor.execute('SELECT * FROM tmp_summary'):
        line = []
        for i in range(1, len(row)):
            line.append(row[i])
        summary.append(line)
        row_pos[row[0]] = row_count
        row_count += 1

    for row in cursor.execute('SELECT a.raw_id, b.name, a.value FROM sample_call a JOIN sample b WHERE a.sample_id = b.pid'):
        calls[row_pos[row[0]]][sample_name_loc[row[1]]] = row[2]
    for summary_line_on in range(0, len(summary)):
        for sample_on in range(0, len(sample_name_loc)):
            summary[summary_line_on].append(calls[summary_line_on, sample_on])
    columns = ['my_id','side','chrom','pos','strand','ref_like_prefix','insertion_after_prefix',
               'insertion_before_suffix','ref_like_suffix','REF','TE']
    for i in range(len(loc_sample_name)):
        columns.append(loc_sample_name[i])

    return pd.DataFrame(summary, columns=columns)


def getPopulation():
    founder_bwt = '/csbiodata/BWTs/sangerMsbwt'
    cc_bwt = '/csbiodata/perlegen/CC_bwt'
    cc_avoid = ['CC039M3730_UNC_NYGC', 'CC019M6839_UNC_UNC', 'CC081M0332_UNC_UNC']
    sanger_founder = ["C57BL6JN", "129S1", "NOD", "NZOdna", "CASTdna", "PWKdna", "WSBdna"]
    unc_founder = ["AJ_F321", "C57BL6J_m003636A", "129S1_M157", "NOD_M146", "NZO_FNNN", "CAST_M090", "PWK_F6114",
                   "WSB_FNNN"]
    cc_samples = []
    #TODO fix this.

    for bwtfile in open('cc_bwt.txt', 'r'):
        sample = bwtfile.split("/")[-1]
        if sample in cc_avoid:
            continue
        cc_samples.append(sample)

    sample_population = {}
    for sample in sanger_founder:
        sample_population[sample] = 'sanger_founder'
    for sample in unc_founder:
        sample_population[sample] = 'unc_founder'
    for sample in cc_samples:
        sample_population[sample] = 'cc'

    return sample_population


# Builds dictionary with the values for chrom, pos, strand, ref_prefix and ref_suffix from Mappable_TEs and
# FinalMatrix_collapsed2.csv. The data for Mappable is in the db and is all that is needed.
def get_ref_dict():
    file1 = "Mappable_TEs.csv"
    file2 = "FinalMatrix_collapsed2.csv"
    ref_dict = {}

    df = pd.read_csv(filepath_or_buffer=file1, sep=',', dtype={'chromo': str, 'pos': int})
    df = df.replace(np.nan, '', regex=True)
    data = df.loc[:, ].values
    for d in data:
        [chromo, pos, strand] = d[0:3]
        [ref_prefix, ref_suffix] = d[6:8]
        ref_dict[(chromo, pos, strand)] = [ref_prefix, ref_suffix]

    # chromo,pos,strand,ref_te,state,before
    df = pd.read_csv(filepath_or_buffer=file2, sep=',', dtype={'chromo': str, 'pos': int})
    df = df.replace(np.nan, '', regex=True)
    data = df.loc[:, ].values
    for d in data:
        [chromo, pos, strand, _, _, ref_prefix, ref_suffix] = d[0:7]
        ref_dict[(chromo, pos, strand)] = [ref_prefix, ref_suffix]
    return ref_dict


# Removes the sequences without ref locations from the ref dictionary produced by get_ref_dict into the scratch file
def remove_not_found_refseq():
    connection.execute('DROP TABLE IF EXISTS tmp_summary')
    connection.execute('CREATE TABLE tmp_summary AS select * FROM raw')
    cursor = connection.cursor()
    ref_dict = get_ref_dict()
    delete_count = 0
    for key in ref_dict.keys():
        delete_count += 1
        cursor.execute('DELETE FROM tmp_summary WHERE chrom = ? AND pos = ? AND strand =?', (key[0], key[1], key[2]))
    print ("Removed [{} locations]".format(delete_count))


def get_counts_te(df):
    dtype = {'chromo': str, 'pos': int}
#    df_te = pd.read_sql_query("", connection)
    # df_te = df_te.replace(np.nan, '', regex=True)
    df_te = df[(df['TE'] == 1)]
    data = df_te.iloc[:, ].values
    ref_dict = get_ref_dict()
    # print ref_dict
    header = list(df_te)

    all_te_context = {}
    sample_list = header[11:]

    # We are building the following dictionaries below
    pos_dict = {}
    pos_to_TE_id = {}
    pos_to_side = {}
    for d in data:
        [my_id, side, chromo, pos, strand] = d[0:5]
        [prefix, te_start, te_end, suffix] = d[5:9]
        original_context = False if len(
            [i for i in prefix + te_start + te_end + suffix if i.islower()]) > 0 else True
        pos_to_TE_id[(chromo, pos, strand)] = my_id
        if original_context:
            pos_to_side[(chromo, pos, strand)] = side

        [ref_prefix, ref_suffix] = ref_dict[(chromo, pos, strand)]
        pid = pos_dict.get((chromo, pos, strand), len(pos_dict))
        pos_dict[(chromo, pos, strand)] = pid
        all_te_context[(chromo, pos, strand)] = [ref_prefix, ref_suffix]

    #We now initialize the
    m, n = len(pos_dict), len(sample_list)
    with_te = np.zeros((m, n), dtype=int)
    without_te = np.zeros((m, n), dtype=int)

    # Populate the with_te with the counts from the samples
    for d in data:
        [chromo, pos, strand] = d[2:5]
        pid = pos_dict[(chromo, pos, strand)]
        sample_count = d[13:]
        for sid, c in enumerate(sample_count):
            with_te[pid][sid] += c

    labels = []
    pid_to_chromo = {v: k for k, v in pos_dict.iteritems()}
    for pid in pid_to_chromo:
        [chromo, pos, strand] = pid_to_chromo[pid]
        l = [chromo, pos, strand]
        labels.append(l)

    dtype = {'chromo': str, 'pos': int}
#    df_te = pd.read_csv(filepath_or_buffer="FinalMatrix_v3.csv", dtype=dtype, sep=',')
#    df_te = df_te.replace(np.nan, '', regex=True)
    df_no_te = df[(df['TE'] == 0)]
    data = df_no_te.iloc[:, ].values

    for d in data:
        [chromo, pos, strand] = d[2:5]
        sample_count = d[13:]
        pid = pos_dict[(chromo, pos, strand)]
        for sid, c in enumerate(sample_count):
            without_te[pid][sid] += c

    return with_te, without_te, labels, sample_list, all_te_context, pos_to_TE_id, pos_to_side


# Uses the interval and sample tables to determine the genotype.
def get_genotype(high_count_strain, chrom, pos, state):
    cursor = connection.cursor()

    chrom = "'%s'" % chrom
    query = """select S.name,I.genotype from Intervals I, Samples S where S.sid=I.sid 
                and I.chromosome = %s and I.start <= %d and I.end >= %d""" % (chrom, pos, pos)
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

    geno_with_TE = set()
    for strain in high_count_strain:
        g = strain_to_genotype.get(strain, 'XX')
        geno_with_TE.add(g)

    fully_consistent_geno = []
    if state == 'private':
        return "%s(%s)" % (list(high_count_strain)[0], list(geno_with_TE)[0])

    geno_with_TE = geno_with_TE - {'XX'}
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


def collapse(df):
    with_te, without_te, labels, sample_list, all_te_context, pos_to_TE_id, pos_to_side = get_counts_te(df)
    (m, n) = with_te.shape

    collapsed_count = np.zeros((m, n), dtype=int) - 1
    sample_population = getPopulation()
    new_data = []
    # m = 10
    # change above line later
    for pid in range(m):
        [chromo, pos, strand] = labels[pid]
        my_id = pos_to_TE_id[(chromo, pos, strand)]
        side = pos_to_side[(chromo, pos, strand)]
        [ref_prefix, ref_suffix] = all_te_context[(chromo, pos, strand)]

        # print sample_count
        high_count = np.where(with_te[pid] > 4)[0]
        high_count_sample = [sample_list[sid] for sid in high_count]

        high_count_strain = set()
        non_cc = False
        for sample in high_count_sample:
            # Paul added as bandaid
            if sample not in sample_population.keys():
                sample_population[sample] = "Unknown"
            if sample_population[sample] == "cc":
                strain = sample[:5]
            else:
                strain = sample
                non_cc = True
            high_count_strain.add(strain)

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

        row_data = [my_id, side, chromo, pos, strand, state, ref_prefix, ref_suffix, geno] + list(
            collapsed_count[pid])
        new_data.append(row_data)

    print collapsed_count

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


if __name__ == '__main__':
    remove_not_found_refseq()
    timestart = datetime.datetime.now()
    print "Begin\t%02d:%02d:%02d" % (timestart.time().hour, timestart.time().minute, timestart.time().second)
    timeloadfin = datetime.datetime.now()
    source_df = build_source_df()
    difference = timeloadfin - timestart
    source_df.to_csv('test.csv')
    collapse(source_df)

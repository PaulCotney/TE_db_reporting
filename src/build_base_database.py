import sqlite3
import itertools
import re
import sys


# This takes a dictionary of values and uses that to create the extension105mer / junct_105_hits / hit tables from
# load tables.
def convert_tmp_full_105_table(args):
    call_columns = dict()
    call105_columns = dict()
    row_columns = cursor_read.execute("SELECT * FROM {}".format(args['table_name'])).fetchone()
    call_columns['start'] = row_columns.keys().index('TE') + 1
    try:
        call_columns['end'] = row_columns.keys().index('tes_i')-1
    except ValueError:
        call_columns['end'] = row_columns.keys().index('tee_i') - 1
    call105_columns['start'] = row_columns.keys().index('side:1') + 1
    call105_columns['end'] = len(row_columns) - 1

    for row in cursor_read.execute("SELECT * FROM {}".format(args['table_name'])):
        pid_105 = cursor_write.execute('INSERT INTO extension105mer(start_or_end, side, chromo, pos, strand, '
                                       'ref_like_prefix_seq, insertion_after_prefix_seq, insertion_before_suffix_seq, '
                                       'ref_like_suffix_seq, TELike_seq, REF) '
                                       'VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)',
                                       (args['start_pos'], row['side'], row['chromo'], row['pos'], row['strand'],
                                        row['ref_like_prefix'], row['insertion_after_prefix'],
                                        row['insertion_before_suffix'], row['ref_like_suffix'],
                                        row['TELike'], row['REF'])).lastrowid
        for sample in itertools.islice(row.keys(), call_columns['start'] + 1, call_columns['end']):
            if int(row[sample]) > 0:
                pid_hit = cursor_write.execute('INSERT INTO hit (target, value) VALUES (?,?)',
                                               (sample, row[sample])).lastrowid
                cursor_write.execute('INSERT INTO junct_105_hits(extension105mer_id, hit_id) VALUES (?,?)',
                                     (pid_105, pid_hit))
        for sample in itertools.islice(row.keys(), call105_columns['start'] + 1, call105_columns['end']):
            if int(row[sample].replace('.0', '')) > 0:
                pid_hit = cursor_write.execute('INSERT INTO hit (target, value) VALUES (?,?)',
                                               (sample.replace(':1', ''), row[sample])).lastrowid
                cursor_write.execute('INSERT INTO junct_105_hits105(extension105mer_id, hit_id) VALUES (?,?)',
                                     (pid_105, pid_hit))


genotype_re = re.compile('\((.*)\)')
cc_re = re.compile('CC...')
founder_exists = re.compile('(AA|BB|CC|DD|EE|FF|HH)')

connection = sqlite3.connect('TE_db.sqlite')
connection.row_factory = sqlite3.Row
cursor_read = connection.cursor()
cursor_read2 = connection.cursor()
cursor_write = connection.cursor()

with open('build_base_database.sql', 'r') as myfile:
    commands = myfile.read().replace('\n', '').split(";")
    for command in commands:
        if len(command.replace(" ", "")) != 0:
            try:
                cursor_write.execute(command + ";")
            except sqlite3.OperationalError as e:
                print("Fail on command: {} with {}", (command, e))
                raise

# Remove items from mappable_TEs where no founder strain is in the genotypes field
rows_to_delete = []
for row in cursor_read.execute('SELECT rowid, genotype FROM load_mappable_tes'):
    if not founder_exists.search(row[1]):
        rows_to_delete.append(row[0])
for row_id in rows_to_delete:
    cursor_write.execute('DELETE FROM load_mappable_tes WHERE rowid = ? ', (int(row_id),))

# Load sample table
strain_codes = {'A': 'A/J',
                'B': 'C57BL/6J',
                'C': '129S1/SvImJ',
                'D': 'NOD/ShiLtJ',
                'E': 'NZO/HlLtJ',
                'F': 'CAST/EiJ',
                'G': 'PWK/PhJ',
                'H': 'WSB/EiJ'}
row = cursor_read.execute('SELECT * FROM load_finalMatrix_ALL').fetchone()
row_columns_start = row.keys().index('TE') + 1
row_columns_end = len(row.keys())
sample_table = {}
iCol = row_columns_start + 1
for sample in itertools.islice(row.keys(), row_columns_start + 1, row_columns_end):
    strain_value = None
    strain_code = None
    long_strain_code = None
    founder = False
    for code, strain_full in strain_codes.iteritems():
        strain_short = strain_full.split("/")[0]
        if code == 'A':
            if sample.find('AJ') != -1:
                strain_value = strain_full
                strain_code = code
                long_strain_code = code + code
                founder = True
                break
        elif sample.find(strain_short) != -1:
            strain_value = strain_full
            strain_code = code
            long_strain_code = code + code
            founder = True
            break
    if strain_value is None:
        strain_value = cc_re.match(sample).group() if cc_re.match(sample) is not None else None
        strain_code = ''
        long_strain_code = ''
    rowid = cursor_write.execute('INSERT INTO SAMPLE (name, strain, founder, strain_code, long_strain_code) values '
                                 '(?,?,?,?,?)',
                                 (sample, strain_value, founder, strain_code, long_strain_code)).lastrowid
    sample_table[sample] = rowid

connection.commit()

# Now with the sample table loaded we create the raw_matrix and the sample_calls table
cursor_read.close()
cursor_read = connection.cursor()
for row in cursor_read.execute('SELECT * FROM load_FinalMatrix_ALL'):
    raw_rowid = cursor_write.execute('INSERT INTO raw (my_id, side, chrom, pos , strand, ref_like_prefix, '
                                     'insertion_after_prefix, insertion_before_suffix, ref_like_suffix, REF, TE) '
                                     'values (?,?,?,?,?,?,?,?,?,?,?)',
                                     (row['my_id'], row['side'], row['chromo'], row['pos'], row['strand'],
                                      row['ref_like_prefix'], row['insertion_after_prefix'],
                                      row['insertion_before_suffix'], row['ref_like_suffix'], row['REF'],
                                      row['TE'])).lastrowid
    for sample in itertools.islice(row.keys(), row_columns_start + 1, row_columns_end):
        if int(row[sample]) > 0:
            cursor_write.execute('INSERT INTO sample_call (raw_id, sample_id, value) VALUES (?,?,?)',
                                 (raw_rowid, sample_table[sample], row[sample]))
connection.commit()
cursor_read.close()
cursor_read = connection.cursor()

# Build genotype helper dictionary (params['potential_genotype'][genotype]) and the genotype table
params = dict()

params['table_name'] = 'tmp_full_105mer_start'
params['start_pos'] = 'start'
convert_tmp_full_105_table(params)
params['table_name'] = 'tmp_full_105mer_end'
params['start_pos'] = 'end'
convert_tmp_full_105_table(params)

# Build Mappable Extensions
for row in cursor_read.execute("SELECT *, rowid FROM load_mappable_tes"):
    totals = {}
    calls = []
    for i in xrange(9, 152):
        calls.append(int(row[i]))
    for call in calls:
        if call in totals.keys():
            totals[call] = totals[call] + 1
        else:
            totals[call] = 1
    for key, value in totals.iteritems():
        cursor_write.execute('INSERT INTO load_map_calls (call_value, number_of_hits, load_mappable_tes_id) '
                             'VALUES (?,?,?)', (key, value, row['rowid']))

# Adding to the report table (TE_in_CC and gene_detail) the values from the derived table, load_mappable_tes.
for row in cursor_read.execute('SELECT * FROM load_mappable_tes'):
    regions = row['region'].split('|')
    gene_ids = row['gene_id'].split('|')
    gene_names = row['gene_name'].split('|')
    assert (len(regions) == len(gene_ids) == len(gene_names))
    gene_detail_rowids = list()
    for i in range(0, len(regions)):
        gene_detail_rowids.append(
            cursor_write.execute('INSERT INTO gene_detail(region, gene_id, gene_name) VALUES (?,?,?)',
                                 (regions[i], gene_ids[i], gene_names[i])).lastrowid)
    genotype_in_founders = None
    SDP = None
    count_in_founders = None
    count_in_CC = None
    te_in_cc_id = cursor_write.execute('INSERT INTO TE_in_CC (my_id, chrom, pos, strand) VALUES (?,?,?,?)',
                                       (row['my_id'], row['chromo'], row['pos'], row['strand'])).lastrowid
    # Write out the join table
    for gene_detail_id in gene_detail_rowids:
        cursor_write.execute('INSERT INTO junct_te_in_cc_gene(gene_detail_id, te_in_cc_id) '
                             'VALUES(?,?)', (gene_detail_id, te_in_cc_id))

rows_to_delete = []
for row_driver in cursor_read.execute('SELECT rowid, my_id, chrom, pos, strand FROM TE_in_CC '):
    if cursor_read2.execute('SELECT count(*) FROM raw where my_id = ? AND chrom = ? AND pos = ? AND strand = ? '
                            'AND TE = 1 ',
                            (row_driver['my_id'], row_driver['chrom'], row_driver['pos'], row_driver['strand'])
                            ).fetchone()[0] == 0:
        # print('{} {} {} {}'.format(row_driver['my_id'], row_driver['chrom'], row_driver['pos'], row_driver['strand']))
        rows_to_delete.append(row_driver[0])
for rowid in rows_to_delete:
    cursor_write.execute('DELETE FROM TE_in_CC where rowid = ?', (rowid,))
connection.commit()

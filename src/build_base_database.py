import sqlite3
import itertools
import re

cc_re = re.compile('CC...')
founder_exists = re.compile('(AA|BB|CC|DD|EE|FF|HH)')



connection = sqlite3.connect('TE_db.sqlite')
connection.row_factory = sqlite3.Row
cursor_read = connection.cursor()
cursor_write = connection.cursor()

# Remove items from mappable_TEs where no founder strain is in the genotypes field
rows_to_delete = []
for row in cursor_read.execute('SELECT rowid, genotype FROM load_mappable_tes'):
    if not founder_exists.search(row[1]):
        rows_to_delete.append(row[0])
for row_id in rows_to_delete:
    cursor_write.execute('DELETE FROM load_mappable_tes WHERE rowid = ? ', (int(row_id),))

# Create the sample table and populate it with information
creation_sql = ('CREATE TABLE sample (pid INTEGER PRIMARY KEY, name TEXT, strain TEXT, founder INTEGER, '
                'strain_code TEXT, long_strain_code TEXT )',)
strain_codes = {'A': 'A/J',
                'B': 'C57BL/6J',
                'C': '129S1/SvImJ',
                'D': 'NOD/ShiLtJ',
                'E': 'NZO/HlLtJ',
                'F': 'CAST/EiJ',
                'G': 'PWK/PhJ',
                'H': 'WSB/EiJ'}

# Load sample table
cursor_write.execute('DROP TABLE IF EXISTS sample')
for command in creation_sql:
    cursor_write.execute(command)
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
cursor_write.execute('DROP TABLE IF EXISTS raw')
cursor_write.execute('CREATE TABLE raw (pid INTEGER PRIMARY KEY, my_id STRING, side STRING, chrom STRING, '
                     'pos INTEGER, strand STRING, ref_like_prefix STRING, insertion_after_prefix STRING, '
                     'insertion_before_suffix STRING, ref_like_suffix STRING, REF STRING, TE STRING)')
cursor_write.execute('DROP TABLE IF EXISTS sample_call')
cursor_write.execute('CREATE TABLE sample_call (raw_id INTEGER, sample_id INTEGER, value INTEGER, FOREIGN KEY(raw_id) '
                     'REFERENCES raw(pid), FOREIGN KEY(sample_id) REFERENCES sample(pid)) ')

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

# We now create a test example with TE_in_CC
cursor_read.close()
cursor_read = connection.cursor()
cursor_read2 = connection.cursor()
cursor_write.execute('DROP TABLE IF EXISTS gene_detail')
cursor_write.execute('CREATE TABLE gene_detail (pid INTEGER PRIMARY KEY, region STRING, gene_id STRING, '
                     'gene_name STRING)')
cursor_write.execute('DROP TABLE IF EXISTS TE_in_CC')
cursor_write.execute('CREATE TABLE TE_in_CC (pid INTEGER PRIMARY KEY, my_id STRING, chrom STRING, pos INTEGER, '
                     'strand STRING, genotype_in_founders STRING, SDP STRING, Count_in_founders INTEGER, '
                     'Count_in_CC INTEGER, gene_matrix_id INTEGER)')
cursor_write.execute('DROP TABLE IF EXISTS gene_matrix')
cursor_write.execute('CREATE TABLE gene_matrix (gene_detail_id INTEGER, te_in_cc_id INTEGER, '
                     'FOREIGN KEY(gene_detail_id) REFERENCES gene_detail(pid), FOREIGN KEY(te_in_cc_id) REFERENCES '
                     'TE_in_CC(pid))')

# Adding to the report table the values from the derived table: load_mappable_tes and build the gene_detail join table
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
    #Write out the join table
    for gene_detail_id in gene_detail_rowids:
        cursor_write.execute('INSERT INTO gene_matrix(gene_detail_id, te_in_cc_id) VALUES(?,?)', (gene_detail_id,
                                                                                                  te_in_cc_id))
connection.commit()

# Delete load tables
delete_list = []
for row in cursor_read.execute("SELECT name FROM sqlite_master WHERE type='table' AND name like 'load_%' ORDER BY name"):
    delete_list.append(row[0])
for delete_item in delete_list:
    cursor_write.execute("DROP TABLE " + delete_item )
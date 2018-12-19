import sqlite3
import itertools
import re

cc_re = re.compile('CC...')

#Create the sample table and populate it with information
creation_sql = ('CREATE TABLE sample (pid INTEGER PRIMARY KEY, name TEXT, strain TEXT, founder INTEGER, '
                'strain_code TEXT, long_strain_code TEXT )',)
strain_codes = {'A' : 'A/J',
                'B' : 'C57BL/6J',
                'C' : '129S1/SvImJ',
                'D' : 'NOD/ShiLtJ',
                'E' : 'NZO/HlLtJ',
                'F' : 'CAST/EiJ',
                'G' : 'PWK/PhJ',
                'H' : 'WSB/EiJ'}

connection = sqlite3.connect('TE_db.sqlite')
connection.row_factory = sqlite3.Row
cursor_read = connection.cursor()
cursor_write = connection.cursor()

#Load sample table
cursor_write.execute('DROP TABLE IF EXISTS sample')
for command in creation_sql:
    cursor_write.execute(command)
row = cursor_read.execute('SELECT * FROM load_finalMatrix_ALL').fetchone()
row_columns_start = row.keys().index('TE') + 1
row_columns_end = len(row.keys())
sample_table = {}
for sample in itertools.islice(row.keys(), row_columns_start + 1, row_columns_end):
    strain_value = None
    strain_code = None
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
        strain_value = cc_re.match(sample).group() if cc_re.match(sample) != None else None
        strain_code = ''
        long_strain_code = ''
    rowid = cursor_write.execute('INSERT INTO SAMPLE (name, strain, founder, strain_code, long_strain_code) values '
                                 '(?,?,?,?,?)', (sample, strain_value, founder, strain_code, long_strain_code)).lastrowid
    connection.commit()
    sample_table[sample] = rowid

#We now create a test example with TE_in_CC

cursor_read.close()
cursor_read = connection.cursor()
cursor_write.execute('DROP TABLE IF EXISTS GENE_DETAIL')
cursor_write.execute('CREATE TABLE GENE_DETAIL (pid INTEGER PRIMARY KEY, region STRING, gene_id STRING, '
                     'gene_name STRING)')
cursor_write.execute('DROP TABLE IF EXISTS TE_in_CC')
cursor_write.execute('CREATE TABLE TE_in_CC (pid INTEGER PRIMARY KEY, my_id STRING, chrom STRING, pos INTEGER, '
                     'strand STRING, genotype_in_founders STRING, SDP STRING, Count_in_founders INTEGER, '
                     'Count_in_CC INTEGER, gene_detail_id INTEGER, FOREIGN KEY (gene_detail_id) REFERENCES '
                     'gene_detail (pid))')
#Loop through the Mappable TEs to build the base. Mappable TEs is derived by what process I don't know
for row in cursor_read.execute('SELECT * FROM load_mappable_tes'):
    cursor_write.execute('INSERT INTO TE_in_CC (my_id, chrom, pos, strand, genotype_in_founders, SDP, count_in_founders,'
                         'count_in_CC values (?,?,?,?,?,?,?)',
                         (row['my_id'], row['chrom'], row['pos'], row['strand'], row['genotype_in_founders'], row['SDP'],
                          row['count_in_founders'], row['Count_in_CC']))
connection.close()
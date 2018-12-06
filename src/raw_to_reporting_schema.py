import sqlite3
import itertools

creation_sql = (
    'CREATE TABLE te ( te_id INTEGER PRIMARY KEY, next_te INTEGER, chromo TEXT, pos TEXT, strand TEXT, '
    'distance TEXT, ref_te TEXT, state TEXT, before TEXT, after TEXT, genotype TEXT, my_id TEXT, region TEXT, '
    'gene_id TEXT, gene_name TEXT, FOREIGN KEY (next_te) REFERENCES te (te_id) )',
    'CREATE TABLE sample (sample_id INTEGER PRIMARY KEY, name TEXT)',
    'CREATE TABLE algorithm_call(te_id INTEGER, call INTEGER, sample_id INTEGER, PRIMARY KEY (te_id, sample_id), '
    'FOREIGN KEY (te_id) REFERENCES te (te_id), FOREIGN KEY (sample_id) REFERENCES sample (sample_id))')

connection = sqlite3.connect('TE_db.sqlite')
connection.row_factory = sqlite3.Row
cursor = connection.cursor()
insert_cursor = connection.cursor()
cursor.execute('DROP TABLE IF EXISTS algorithm_call')
cursor.execute('DROP TABLE IF EXISTS te')
cursor.execute('DROP TABLE IF EXISTS sample')
for command in creation_sql:
    cursor.execute(command)
row = cursor.execute('SELECT * FROM load_mappable_tes').fetchone()
row_columns_start = row.keys().index('genotype') + 1
row_columns_end = row.keys().index('my_id') - 1
#Create the sample table and fill it with the correct values
sample_table = {}
for sample in itertools.islice(row.keys(), row_columns_start + 1, row_columns_end):
    rowid = insert_cursor.execute('INSERT INTO SAMPLE (name) values (?)', (sample,)).lastrowid
    sample_table[sample] = rowid
rowid = None
old_row_id = None
row_pairs = {}
for row in cursor.execute('SELECT * FROM load_mappable_tes'):
    te_insert_row = []
    algorithm_calls = {}
    for col in itertools.islice(row, 0, row_columns_start):
        te_insert_row.append(col)
    for col in itertools.islice(row, row_columns_end + 1, len(row.keys())):
        te_insert_row.append(col)
    for col in itertools.islice(row.keys(), row_columns_start + 1, row_columns_end):
        algorithm_calls[col] = row[col]
    if rowid:
        old_row_id = rowid
    rowid = insert_cursor.execute("INSERT INTO te (chromo,pos,strand,distance,ref_te,state,before,"
                                  "after,genotype,my_id,region,gene_id,gene_name ) VALUES "
                                  "(?,?,?,?,?,?,?,?,?,?,?,?,?) ", te_insert_row).lastrowid
    if old_row_id:
        row_pairs[old_row_id] = rowid
    for sample in itertools.islice(row.keys(), row_columns_start + 1, row_columns_end):
        insert_cursor.execute('INSERT INTO algorithm_call (te_id, sample_id, call) VALUES (?,?,?)', (rowid,
                                                                                              sample_table[sample],
                                                                                              algorithm_calls[sample]))
for key, value in row_pairs.items():
    insert_cursor.execute("UPDATE te SET next_te = ? WHERE te_id = ?", (value, key))

insert_cursor.execute('commit')

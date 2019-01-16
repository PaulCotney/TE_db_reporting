import sqlite3

connection = sqlite3.connect('TE_db.sqlite')
connection.row_factory = sqlite3.Row
cursor_read = connection.cursor()
cursor_read2 = connection.cursor()
cursor_write = connection.cursor()

connection.commit()
for row_driver in cursor_read.execute('SELECT my_id, chrom, pos, strand FROM TE_in_CC '):
    if cursor_read2.execute('SELECT count(*) FROM raw where my_id = ? AND chrom = ? AND pos = ? AND strand = ? '
                            'AND insertion_before_suffix = '' AND ref_like_suffix = '' ',
                            (row_driver['my_id'], row_driver['chrom'], row_driver['pos'], row_driver['strand'])
                            ).fetchone()[0] != 1:
        assert False


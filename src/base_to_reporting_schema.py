import sqlite3

connection = sqlite3.connect('TE_db.sqlite')
connection.row_factory = sqlite3.Row
cursor_read = connection.cursor()
cursor_read2 = connection.cursor()
cursor_write = connection.cursor()

# Remove items from TE_in_CC that have no raw value.
rows_to_delete = []
for row_driver in cursor_read.execute('SELECT rowid, my_id, chrom, pos, strand FROM TE_in_CC '):
    if cursor_read2.execute('SELECT count(*) FROM raw where my_id = ? AND chrom = ? AND pos = ? AND strand = ? '
                            'AND TE = 1 ',
                            (row_driver['my_id'], row_driver['chrom'], row_driver['pos'], row_driver['strand'])
                            ).fetchone()[0] == 0:
        print('{} {} {} {}'.format(row_driver['my_id'], row_driver['chrom'], row_driver['pos'], row_driver['strand']))
        rows_to_delete.append(row_driver[0])
for rowid in rows_to_delete:
    cursor_write.execute('DELETE FROM TE_in_CC where rowid = ?', (rowid,) )
connection.commit()
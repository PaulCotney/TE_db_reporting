import sqlite3

connection = sqlite3.connect('TE_db.sqlite')
connection.row_factory = sqlite3.Row
cursor_read = connection.cursor()
cursor_read2 = connection.cursor()
cursor_write = connection.cursor()

cursor_write.execute('DROP TABLE IF EXISTS load_map_calls')
cursor_write.execute('CREATE TABLE load_map_calls (call_value INTEGER, number_of_hits INTEGER, load_mappable_tes_id INTEGER, FOREIGN KEY(load_mappable_tes_id) REFERENCES load_mappable_tes(rowiid))')


for row in cursor_read.execute("SELECT *,rowid FROM load_mappable_tes"):
    totals = {}
    calls = []
    for i in xrange(9,152):
        calls.append(int(row[i]))
    for call in calls:
        if call in totals.keys():
            totals[call] = totals[call] + 1
        else:
            totals[call] = 1
    for key, value in totals.iteritems():
        cursor_write.execute('INSERT INTO load_map_calls (call_value, number_of_hits, load_mappable_tes_id) '
                             'VALUES (?,?,?)', (key, value, row['rowid']))
connection.commit()
#    for i in calls:

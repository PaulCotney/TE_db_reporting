import sqlite3

delete_list = []
connection = sqlite3.connect('TE_db.sqlite')
connection.row_factory = sqlite3.Row
cursor_read = connection.cursor()
cursor_read2 = connection.cursor()
cursor_write = connection.cursor()
for row in cursor_read.execute(
        "SELECT name FROM sqlite_master WHERE type='table' AND ( name like 'load_%' OR name like "
        "'tmp_%' ) ORDER BY name"):
    delete_list.append(row[0])
for delete_item in delete_list:
    cursor_write.execute("DROP TABLE " + delete_item)
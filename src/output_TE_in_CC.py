import sqlite3

connection = sqlite3.connect('TE_db.sqlite')
connection.row_factory = sqlite3.Row
cursor_read = connection.cursor()
#Write out table
out_file = open('TE_in_CC.csv','w+')
# write header
select_sql = 'SELECT my_id, chrom, pos, strand, genotype_in_founders, SDP, count_in_founders, count_in_CC ' \
             'FROM TE_in_CC'

row_columns = cursor_read.execute(select_sql) .fetchone()
row_text = []
for column in row_columns.keys():
    row_text.append(column)
for row_cell in row_text[:-1]:
    out_file.write("{},".format(row_cell))
out_file.write("{}\n".format(row_text[-1]))
for row in cursor_read.execute(select_sql):
    row_text = []
    for column in row:
        row_text.append(column)
    for row_cell in row_text[:-1]:
        out_file.write("{},".format(row_cell))
    out_file.write("{}\n".format(row_text[-1]))
out_file.close()
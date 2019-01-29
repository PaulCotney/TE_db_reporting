import sqlite3

connection = sqlite3.connect('TE_db.sqlite')
connection.row_factory = sqlite3.Row
cursor_read = connection.cursor()
cursor_read2 = connection.cursor()
cursor_write = connection.cursor()
csv_out_file = open("Duplicates.csv", "w+")

def printorder(columns, skipdistance = False):
    text_to_write = []
    out_line = ""
    if not skipdistance:
        for column_on in columns[0:9]:
            if column_on != '' and column_on != '-' and column_on[0] == '-' :
                text_to_write.append("'{},".format(column_on))
            else:
                text_to_write.append("{},".format(column_on))
    else:
        for column_on in columns[0:3]:
            if column_on != '' and column_on != '-'  and column_on[0] == '-':
                text_to_write.append("'{},".format(column_on))
            else:
                text_to_write.append("{},".format(column_on))
        text_to_write.append(",")
        for column_on in columns[4:9]:
            if column_on != '' and column_on != '-' and column_on[0] == '-':
                text_to_write.append("'{},".format(column_on))
            else:
                text_to_write.append("{},".format(column_on))
    for column_on in columns[152:160]:
        if column_on == 'None':
            text_to_write.append("0,")
        else:
            text_to_write.append("{},".format(column_on))
    for column_on in columns[9:151]:
        text_to_write.append("{},".format(column_on))
    text_to_write.append("{}\n".format(columns[151]))
    for cell in text_to_write:
        out_line = out_line + cell
    csv_out_file.write(out_line)
    print(out_line)

row_columns = cursor_read.execute("SELECT * FROM v_mappable_tes_ratios") .fetchone()
printorder(row_columns.keys())


row_previous = None
for row in cursor_read.execute("SELECT * FROM v_mappable_tes_ratios order by orig_order"):
    if not row_previous:
        row_previous = row
        continue

    if int(row_previous['distance']) < 4000 and \
            row_previous['strand'] == row['strand'] and row_previous['calls_1'] < 10 and row['calls_1'] < 10:
        row_text = []
        for column in row_previous:
            row_text.append(str(column))
        printorder(row_text)
        print(row_text)
        row_text = []
        for column in row:
            row_text.append(str(column))
        printorder(row_text, True)
        print(row_text)

        # spacer
        csv_out_file.write('\n')

    row_previous = row

import sqlite3
import itertools
import re


connection = sqlite3.connect('TE_db.sqlite')
connection.row_factory = sqlite3.Row
cursor_read = connection.cursor()
cursor_read2 = connection.cursor()

csv_out_file = open("collapse_master_detail.csv", "w+")
# Remove items from mappable_TEs where no founder strain is in the genotypes field

csv_out_file.write("collapsed,my_id,side,chrom,pos,strand,state,before,after,ref_like_prefix,insertion_after_prefix,"
                   "insertion_before_suffix,ref_like_suffix,REF,TE,genotype,region,gene_id,gene_name")
sample_columns = []
row_columns = cursor_read.execute("SELECT * FROM load_Final_Matrix_ALL_collapsed") .fetchone()
for column in row_columns.keys()[8:116]:
    csv_out_file.write("{},".format(column))
csv_out_file.write("{}\n".format(row_columns.keys()[116]))
for row in cursor_read.execute('SELECT * FROM load_Final_Matrix_ALL_collapsed'):
    row_text = []
    analysis_cols = []
    for column in row:
        row_text.append(column)
    analysis_cols.extend(row_text[:8])
    analysis_cols.append("")
    analysis_cols.append("")
    analysis_cols.append("")
    analysis_cols.append("")
    analysis_cols.append("")
    analysis_cols.append("")
    analysis_cols.append(row_text[8])
    analysis_cols.extend(row_text[117:])
    # THe collapsed boolean field
    analysis_cols.insert(0, "TRUE")
    analysis_cols.extend(row_text[9:117])
    for cell in analysis_cols[:-1]:
        csv_out_file.write("{},".format(cell))
    csv_out_file.write("{}\n".format(analysis_cols[-1]))
    for row_raw in cursor_read2.execute('SELECT * FROM load_finalmatrix_all WHERE my_id = ? AND chromo = ? AND'
                                       ' pos = ?', (row['my_id'], row['chromo'], row['pos'])):
        row_text =[]
        analysis_cols = []
        for column in row_raw:
            row_text.append(column)
        analysis_cols.extend(row_text[:5])
        analysis_cols.append("")
        analysis_cols.append("")
        analysis_cols.append("")
        analysis_cols.extend(row_text[5:11])
        analysis_cols.append("")
        analysis_cols.append("")
        analysis_cols.append("")
        analysis_cols.append("")
        analysis_cols.extend(row_text[11:])
        analysis_cols.insert(0, "FALSE")
        for cell in analysis_cols[:-1]:
            csv_out_file.write("{},".format(cell))
        csv_out_file.write("{}\n".format(analysis_cols[-1]))


import sqlite3
import sys
from itertools import islice
from string import ascii_uppercase

connection = sqlite3.connect('TE_db.sqlite')
connection.row_factory = sqlite3.Row
cursor_read = connection.cursor()
cursor_read2 = connection.cursor()

out_file = open('test_105mer.csv','w+')

# Build a list of 'AA' 'BB' ..
target_genotype = list()
for c in islice(ascii_uppercase, 0, 8):
    target_genotype.append('{0}{0}'.format(c))
# Choose some targets
target_sample = ('AJ_F321', 'C57BL6J_m003636A', '129S1_M157', 'NZO_FNNN', 'CAST_M090', 'PWK_F6114', 'NOD_M146')

# Main driver sql statement for the starting 105mer
# TODO get rid of limit
start_driver_sql = 'SELECT a.pid, a.side, a.chromo, a.pos, a.strand, a.ref_like_prefix_seq, ' \
             'a.insertion_after_prefix_seq, a.insertion_before_suffix_seq, a.ref_like_suffix_seq, a.TELike_seq, a.REF, ' \
             'a.state FROM extension105mer a WHERE a.start_or_end = \'start\' limit 50 '
# Sql to return the genotype values for the driver row
genotype_sql = 'SELECT a.value, a.short_value FROM genotype a INNER JOIN junct_105_genotype b ON a.pid = ' \
                 'b.genotype_id WHERE b.extension105mer_id = ?'

where_in_clause = "("
for target in target_sample[:-1]:
    where_in_clause = where_in_clause + '\'{}\', '.format(target)
where_in_clause = where_in_clause + '\'{}\') '.format(target)

hits_105_sql = 'SELECT a.target, a.value FROM hit a INNER JOIN junct_105_hits b ON a.pid = b.hit_id WHERE ' \
           'b.extension105mer_id = ? AND a.target IN {}'.format(where_in_clause)
hits_105ext_sql = 'SELECT a.target, a.value FROM hit a INNER JOIN junct_105_hits105 b ON a.pid = b.hit_id WHERE ' \
              'b.extension105mer_id = ? AND a.target IN {}'.format(where_in_clause)
row_columns = cursor_read.execute(start_driver_sql).fetchone()
row_text = []

# write header
# write columns to list and then list to outfile
for column in row_columns.keys()[1:]:
    row_text.append(column)
for strain in target_genotype:
    row_text.append(strain)
for sample in target_sample:
    row_text.append(sample)
for row_cell in row_text:
    out_file.write("{},".format(row_cell))
out_file.write(("{}\n".format(row_text[-1])))
# Loop through each line in extensions105mer selecting genotypes we want
for row in cursor_read.execute(start_driver_sql):
    line_out = list()
    row_text = list()
    for column in row:
        row_text.append(column)
    for row_cell in row_text[1:]:
        line_out.append("{}".format(row_cell))

    # Add the strain columns for filtering
    row_has_strain = list()
    for subrow in cursor_read2.execute(genotype_sql, (row['pid'],)):
        if subrow['short_value'] in target_genotype:
            row_has_strain.append(subrow['short_value'])
    row_has_any_strain = False
    for target_strain in target_genotype:
        if target_strain in row_has_strain:
            line_out.append("TRUE".format(row_cell))
            row_has_any_strain = True
        else:
            line_out.append("FALSE".format(row_cell))

    if not row_has_any_strain:
        # The row has no strain from target strain list so exit
        continue

    # Write hits
    # TODO failing here
    hits_values = dict()
    for subrow in cursor_read2.execute(hits_105_sql, (row['pid'],)):
        hits_values[subrow['target']] = subrow['value']
    for sample in target_sample:
        line_out.append(hits_values.get(sample, 0))
    for row_cell in line_out:
        out_file.write("{},".format(row_cell))
    out_file.write(("{}\n".format(line_out[-1])))
out_file.close()
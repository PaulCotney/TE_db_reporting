import sqlite3
import itertools
import re

genotype_re = re.compile('\((.*)\)')
connection = sqlite3.connect('TE_db.sqlite')
connection.row_factory = sqlite3.Row
cursor_read = connection.cursor()
cursor_write = connection.cursor()


def convert_tmp_full_105_table(args):
    call_columns = dict()
    call105_columns = dict()
    row_columns = cursor_read.execute("SELECT * FROM {}".format(args['table_name'])).fetchone()
    call_columns['start'] = row_columns.keys().index('genotype') + 1
    call_columns['end'] = row_columns.keys().index('my_id') - 1
    call105_columns['start'] = row_columns.keys().index('side:1') + 1
    call105_columns['end'] = len(row_columns) - 1

    for row in cursor_read.execute("SELECT * FROM {}".format(args['table_name'])):
        pid_105 = cursor_write.execute('INSERT INTO extension105mer(start_or_end, side, chromo, pos, strand, '
                                      'ref_like_prefix_seq, insertion_after_prefix_seq, insertion_before_suffix_seq, '
                                      'ref_like_suffix_seq, TELike_seq, REF, state) '
                                      'VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)',
                                      (args['start_pos'], row['side'], row['chromo'], row['pos'], row['strand'],
                                       row['ref_like_prefix'], row['insertion_after_prefix'],
                                       row['insertion_before_suffix'], row['ref_like_suffix'],
                                       row['TELike'], row['REF'], row['state'])).lastrowid
        for sample in itertools.islice(row.keys(), call_columns['start'] + 1, call_columns['end']):
            if int(row[sample]) > 0:
                pid_hit = cursor_write.execute('INSERT INTO hit (target, value) VALUES (?,?)',
                                     (row_columns[sample], row[sample])).lastrowid
                cursor_write.execute('INSERT INTO junct_105_hits(extension105mer_id, hit_id) VALUES (?,?)',
                                     (pid_105, pid_hit))
        for sample in itertools.islice(row.keys(), call105_columns['start'] + 1, call105_columns['end']):
            if int(row[sample].replace('.0', '')) > 0:
                pid_hit = cursor_write.execute('INSERT INTO hit (target, value) VALUES (?,?)',
                                     (row_columns[sample].replace(':0', ''), row[sample])).lastrowid
                cursor_write.execute('INSERT INTO junct_105_hits105(extension105mer_id, hit_id) VALUES (?,?)',
                                     (pid_105, pid_hit))
        for genotype_call in row['genotype'].split('|'):
            pid_hit = args['potential_genotype'][genotype_call]
            cursor_write.execute('INSERT INTO junct_105_genotype (extension105mer_id, genotype_id) VALUES'
                                 ' (?,?)', (pid_105, pid_hit))



with open('build_clean_105mers.sql', 'r') as myfile:
    commands = myfile.read().replace('\n', '').split(";")
    for command in commands:
        if len(command.replace(" ","")) <> 0:
            try :
                cursor_write.execute(command + ";")
            except sqlite3.OperationalError as e:
                print("Fail on command: {} with {}", (command, e ))
                raise

# Build genotype and load helper dictionary
params = dict()
params['potential_genotype'] = dict()
genotype_raw = []
for genotype in cursor_read.execute("SELECT DISTINCT genotype FROM (SELECT genotype FROM tmp_full_105mer_start UNION SELECT "
                               "genotype FROM tmp_full_105mer_end) "):
    genotype_raw.append(genotype)
genotypes = set()
for raw in genotype_raw:
    for item in raw['genotype'].split('|'):
        genotypes.add(item)
for genotype in genotypes:
    instance = genotype_re.search(genotype)
    if instance is None:
        pid = cursor_write.execute('INSERT INTO genotype ( value, short_value) VALUES (?, ?) ',
                                   (genotype, genotype.replace('*', ''))).lastrowid
    else:
        pid = cursor_write.execute('INSERT INTO genotype ( value, short_value) VALUES (?, ?) ',
                                   (genotype, instance.group(1).replace('*', ''))).lastrowid

    params['potential_genotype'][genotype] = pid


params['table_name'] = 'tmp_full_105mer_start'
params['start_pos'] = 'start'
convert_tmp_full_105_table(params)
params['table_name'] = 'tmp_full_105mer_end'
params['start_pos'] = 'end'
convert_tmp_full_105_table(params)
connection.commit()
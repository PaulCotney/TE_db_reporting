import sqlite3


def write_row_to_file(*argv):
    pass


def write_row_to_db(row):
    cursor_write.execute('INSERT INTO summary (my_id, chrom, pos, strand, before, after) VALUES (?, ?, ?, ?, ?, ?)',
                         (row['my_id'], row['chrom'], row['pos'], row['strand'], row['before'], row['after']))


def process_summary_line_with_detail(**kwargs):
    if kwargs['summary_field'] is not None and kwargs['summary_field'] != '' and kwargs['summary_field'] == \
            kwargs['detail_field']:
        pass
    elif kwargs['summary_field'] is not None and (
            kwargs['summary_field'] == '' and kwargs['summary_field'] != kwargs['detail_field']):
        return ''
    elif kwargs['summary_field'] is None and kwargs['detail_field'] != '':
        return kwargs['detail_field']
    elif kwargs['summary_field'] is None and kwargs['detail_field'] == '':
        pass
    else:
        assert True


connection = sqlite3.connect('TE_db.sqlite')
connection.row_factory = sqlite3.Row
cursor_read = connection.cursor()
cursor_read2 = connection.cursor()
cursor_write = connection.cursor()

cursor_write.execute('DELETE FROM summary')
connection.commit()
for row in cursor_read.execute('SELECT my_id, chrom, pos, strand FROM raw GROUP BY my_id, chrom, pos, strand ORDER BY '
                               'chrom, pos limit 5'):
    detail_matrix = list()
    result = {'my_id' : row['my_id'], 'chrom' : row['chrom'], 'pos' : row['pos'], 'strand' : row['strand'],
              'before' : None, 'after' : None}
    for detail_row in cursor_read2.execute('SELECT side, ref_like_prefix, insertion_after_prefix, '
                                           'insertion_before_suffix, ref_like_suffix, REF FROM raw WHERE my_id = ? AND '
                                           'chrom = ? AND pos = ? AND strand = ?',
                                           (row['my_id'], row['chrom'], row['pos'], row['strand'])):
        if detail_row['REF'] == 1:
            result['before' ] = process_summary_line_with_detail(summary_field=result['before'], detail_field=detail_row['ref_like_prefix'])
            result['after'] = process_summary_line_with_detail(summary_field=result['after'], detail_field=detail_row['ref_like_suffix'])
    write_row_to_db(result)
connection.commit()
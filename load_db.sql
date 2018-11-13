.mode csv
.import anwica_table_start.clean.csv  load_start_anwica
UPDATE load_start_anwica SET genotype = replace( genotype, '-', '|' );
.import join_table_start.csv load_start_join
.import indexed_merged_te_cnts_start.csv load_start_indexed
.import anwica_table_end.clean.csv  load_end_anwica
UPDATE load_end_anwica SET genotype = replace( genotype, '-', '|' );
.import join_table_end.csv load_end_join
.import indexed_merged_te_cnts_end.csv load_end_indexed
.mode csv
.import ./raw_csv/anwica_table_start.clean.csv  load_start_anwica
UPDATE load_start_anwica SET genotype = replace( genotype, '-', '|' );
.import ./raw_csv/join_table_start.csv load_start_join
.import ./raw_csv/indexed_merged_te_cnts_start.csv load_start_indexed
.import ./raw_csv/anwica_table_end.clean.csv  load_end_anwica
UPDATE load_end_anwica SET genotype = replace( genotype, '-', '|' );
.import ./raw_csv/join_table_end.csv load_end_join
.import ./raw_csv/indexed_merged_te_cnts_end.csv load_end_indexed
.import ./raw_csv/Mappable_TEs.clean.csv load_mappable_tes
.import ./raw_csv/denovoTE105_start.clean.csv load_denovo_te105_start
UPDATE load_denovo_te105_start SET genotype = replace( genotype, '-', '|' );
.import ./raw_csv/denovoTE105_end.clean.csv load_denovo_te105_end
UPDATE load_denovo_te105_end SET genotype = replace( genotype, '-', '|' );
.import ./raw_csv/FinalMatrix_ALL.csv load_FinalMatrix_ALL
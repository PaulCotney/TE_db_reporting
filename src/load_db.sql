.mode csv
.import ./raw_csv/105mers_start.csv load_105mers_start
.import ./raw_csv/105mers_end.csv load_105mers_end
.import ./raw_csv/anwica_table_end.clean.csv load_anwica_table_end
.import ./raw_csv/anwica_table_start.clean.csv load_anwica_table_start
.import ./raw_csv/join_table_end.csv load_join_table_end
.import ./raw_csv/join_table_start.csv load_join_table_start

.import ./raw_csv/Mappable_TEs.clean.csv load_mappable_tes
UPDATE load_mappable_tes SET genotype = replace( genotype, '-', '|' );
UPDATE load_anwica_table_end SET genotype = replace( genotype, '-', '|' );
UPDATE load_anwica_table_start SET genotype = replace( genotype, '-', '|' );
.import ./raw_csv/FinalMatrix_ALL.csv load_FinalMatrix_ALL
CREATE TABLE tmp_full_105mer_start
AS SELECT *
FROM load_anwica_table_start a, load_join_table_start b, load_105mers_start c
WHERE a.as_i = b.as_i AND b.tes_i = c.tes_i;
CREATE TABLE tmp_full_105mer_end
AS SELECT *
FROM load_anwica_table_end a, load_join_table_end b, load_105mers_end c
WHERE a.ae_i = b.ae_i AND b.tee_i = c.tee_i;

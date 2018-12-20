.mode csv
.import ./raw_csv/Mappable_TEs.clean.csv load_mappable_tes
UPDATE load_mappable_tes SET genotype = replace( genotype, '-', '|' );
.import ./raw_csv/FinalMatrix_ALL.csv load_FinalMatrix_ALL
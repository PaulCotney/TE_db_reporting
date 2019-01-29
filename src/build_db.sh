 #!/bin/bash
rm TE_db.sqlite
touch TE_db.sqlite
# Fic the column that has pipes in it for separation of items in an array.
perl -p -e 's/\|/-/g' raw_csv/Mappable_TEs.csv > raw_csv/Mappable_TEs.clean.csv
perl -p -e 's/\|/-/g' raw_csv/anwica_table_start.csv > raw_csv/anwica_table_start.clean.csv
perl -p -e 's/\|/-/g' raw_csv/anwica_table_end.csv > raw_csv/anwica_table_end.clean.csv

# Loads the tables above into the database.
sqlite3 TE_db.sqlite < load_db.sql

# Cleanup of intermediate  files no longer
rm raw_csv/*.clean.csv

#Note Join works
#select a.*, c.* from start_anwica a, start_join b, start_indexed c WHERE a.as_i = b.as_i AND b.tes_i = c.tes_i; 

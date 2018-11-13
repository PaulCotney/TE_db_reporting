#!/bin/bash

# Fic the column that has pipes in it for seperation of items in an array.
perl -p -e 's/\|/-/g' anwica_table_start.csv > anwica_table_start.clean.csv
perl -p -e 's/\|/-/g' anwica_table_end.csv > anwica_table_end.clean.csv

# Loads the tables above into the database.
sqlite3 TE_db.sqlite < load_db.sql

# Cleanup of intermediate  files no longer
rm *.clean.csv

#Note Join works
#select a.*, c.* from start_anwica a, start_join b, start_indexed c WHERE a.as_i = b.as_i AND b.tes_i = c.tes_i; 

DROP TABLE IF EXISTS gene_matrix;
DROP TABLE IF EXISTS TE_in_CC;
DROP TABLE IF EXISTS gene_detail;
DROP TABLE IF EXISTS sample_call;
DROP TABLE IF EXISTS raw;
DROP TABLE IF EXISTS sample;
CREATE TABLE sample (pid INTEGER PRIMARY KEY, name TEXT, strain TEXT, founder INTEGER, strain_code TEXT, long_strain_code TEXT );
CREATE TABLE raw (pid INTEGER PRIMARY KEY, my_id STRING, side STRING, chrom STRING,
pos INTEGER, strand STRING, ref_like_prefix STRING, insertion_after_prefix STRING, 
insertion_before_suffix STRING, ref_like_suffix STRING, REF STRING, TE STRING);
CREATE TABLE sample_call (raw_id INTEGER, sample_id INTEGER, value INTEGER, FOREIGN KEY(raw_id)
REFERENCES raw(pid), FOREIGN KEY(sample_id) REFERENCES sample(pid));
CREATE TABLE gene_detail (pid INTEGER PRIMARY KEY, region STRING, gene_id STRING, 
gene_name STRING);
CREATE TABLE TE_in_CC (pid INTEGER PRIMARY KEY, my_id STRING, chrom STRING, pos INTEGER,
strand STRING, genotype_in_founders STRING, SDP STRING, Count_in_founders INTEGER,
Count_in_CC INTEGER, gene_matrix_id INTEGER);
CREATE TABLE gene_matrix (gene_detail_id INTEGER, te_in_cc_id INTEGER,
FOREIGN KEY(gene_detail_id) REFERENCES gene_detail(pid), FOREIGN KEY(te_in_cc_id) REFERENCES 
TE_in_CC(pid));
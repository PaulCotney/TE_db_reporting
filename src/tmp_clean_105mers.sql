DROP TABLE IF EXISTS junct_105_genotype;
DROP TABLE IF EXISTS junct_105_hits105;
DROP TABLE IF EXISTS junct_105_hits;
DROP TABLE IF EXISTS hit;
DROP TABLE IF EXISTS genotype;
DROP TABLE IF EXISTS extension105mer;
CREATE TABLE extension105mer ( 
pid INTEGER PRIMARY KEY, 
start_or_end STRING,
side STRING,
chromo STRING,
pos INTEGER,
strand STRING,
ref_like_prefix_seq STRING,
insertion_after_prefix_seq STRING,
insertion_before_suffix_seq STRING,
ref_like_suffix_seq STRING,
TELike_seq STRING,
REF INTEGER,
state STRING);
CREATE TABLE hit (pid INTEGER PRIMARY KEY, target STRING, value INTEGER);
CREATE TABLE genotype (pid INTEGER PRIMARY KEY, value STRING, short_value STRING);
CREATE TABLE junct_105_hits (extension105mer_id INTEGER, hit_id INTEGER,
	FOREIGN KEY(extension105mer_id) REFERENCES extension105mer(pid), 
	FOREIGN KEY(hit_id) REFERENCES hit(pid)
	);
CREATE TABLE junct_105_hits105 (extension105mer_id INTEGER, hit_id INTEGER,
	FOREIGN KEY(extension105mer_id) REFERENCES extension105mer(pid), 
	FOREIGN KEY(hit_id) REFERENCES hit(pid)
	);	
CREATE TABLE junct_105_genotype (extension105mer_id INTEGER, genotype_id INTEGER,
	FOREIGN KEY(extension105mer_id) REFERENCES extension105mer(pid), 
	FOREIGN KEY(genotype_id) REFERENCES genotype(pid)
	)

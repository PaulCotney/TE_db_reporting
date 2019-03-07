DROP TABLE IF EXISTS junct_te_in_cc_gene;
DROP TABLE IF EXISTS TE_in_CC;
DROP TABLE IF EXISTS gene_detail;
DROP TABLE IF EXISTS junct_sample_TE;
DROP TABLE IF EXISTS TE;
DROP TABLE IF EXISTS sample;
DROP TABLE IF EXISTS junct_105_genotype;
DROP TABLE IF EXISTS junct_105_hits105;
DROP TABLE IF EXISTS junct_105_hits;
DROP TABLE IF EXISTS hit;
DROP TABLE IF EXISTS genotype;
DROP TABLE IF EXISTS extension105mer;
DROP TABLE IF EXISTS load_map_calls;
DROP TABLE IF EXISTS collapsed_TE;
DROP TABLE IF EXISTS junct_collapsed_TE_junct_sample_TE;
DROP TABLE IF EXISTS junct_collapsed_TE_sample;
CREATE TABLE sample (pid INTEGER PRIMARY KEY, name TEXT, strain TEXT, founder INTEGER,
    strain_code TEXT, long_strain_code TEXT, bits_per_byte FLOAT, bases INTEGER, bytes_compressed INTEGER );
CREATE TABLE TE (pid INTEGER PRIMARY KEY, my_id STRING, side STRING, chrom STRING,
    pos INTEGER, strand STRING, ref_like_prefix STRING, insertion_after_prefix STRING,
    insertion_before_suffix STRING, ref_like_suffix STRING, REF STRING, TE STRING);
CREATE TABLE junct_sample_TE (pid INTEGER PRIMARY KEY, TE_id INTEGER, sample_id INTEGER, value INTEGER,
    FOREIGN KEY(TE_id) REFERENCES TE(pid), FOREIGN KEY(sample_id) REFERENCES sample(pid));
CREATE TABLE gene_detail (pid INTEGER PRIMARY KEY, region STRING, gene_id STRING, 
    gene_name STRING);
CREATE TABLE TE_in_CC (pid INTEGER PRIMARY KEY, my_id STRING, chrom STRING, pos INTEGER,
    strand STRING, genotype_in_founders STRING, SDP STRING, Count_in_founders INTEGER,
    Count_in_CC INTEGER, junct_te_in_cc_gene_id INTEGER);
CREATE TABLE junct_te_in_cc_gene (gene_detail_id INTEGER, te_in_cc_id INTEGER,
    FOREIGN KEY(gene_detail_id) REFERENCES gene_detail(pid), FOREIGN KEY(te_in_cc_id) REFERENCES
    TE_in_CC(pid));
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
	);
CREATE TABLE load_map_calls (call_value INTEGER, number_of_hits INTEGER, load_mappable_tes_id INTEGER, FOREIGN KEY(load_mappable_tes_id) REFERENCES load_mappable_tes(rowiid));
CREATE TABLE collapsed_TE (
    pid  INTEGER PRIMARY KEY,
    state TEXT,
    before TEXT,
    after TEXT,
    genotype TEXT,
    region TEXT,
    gene_id TEXT,
    gene_namegenotype TEXT);
CREATE TABLE junct_collapsed_TE_junct_sample_TE (collapsed_TE_id INTEGER, junct_sample_TE_id INTEGER, value INTEGER,
    FOREIGN KEY(collapsed_TE_id) REFERENCES collapsed_TE(pid),
    FOREIGN KEY(junct_sample_TE_id) REFERENCES junct_sample_TE(pid));
CREATE TABLE junct_collapsed_TE_sample (collapsed_TE_id INTEGER, sample_id INTEGER, value INTEGER,
    FOREIGN KEY(collapsed_TE_id) REFERENCES collapsed_TE(pid),
    FOREIGN KEY(sample_id) REFERENCES sample(pid));

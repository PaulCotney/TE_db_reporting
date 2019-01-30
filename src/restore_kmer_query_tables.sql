CREATE TABLE kmer_query (raw_id INTEGER, dataset STRING, forward INTEGER, reverse INTEGER, FOREIGN KEY(raw_id) REFERENCES raw(pid));
INSERT INTO kmer_query (raw_id,dataset,forward,reverse)
SELECT b.pid, a.dataset, a.forward, a.reverse
FROM kmer_query_export a INNER JOIN raw b ON a.ref_like_prefix = b.ref_like_prefix
 AND a.insertion_after_prefix = b.insertion_after_prefix AND a.insertion_before_suffix = b.insertion_before_suffix
 AND a.ref_like_suffix = b.ref_like_suffix AND a.pos = b.pos
 GROUP BY b.pid, a.dataset, a.forward, a.reverse;

SELECT DISTINCT cnt FROM
 ( select raw_id, dataset, count(*) cnt FROM kmer_query GROUP BY raw_id, dataset ) a
 WHERE a.cnt > 1
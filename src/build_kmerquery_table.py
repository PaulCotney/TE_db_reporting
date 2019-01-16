from MUSCython import MultiStringBWTCython as MSBWT
import time
import sqlite3

# Yanked from import msSharedUtil
dirLabels = ['CEGS 3x3 diallele',
             'CC Founders',
             'WashU HR Mice',
             'Wild Testis RNA-seq (PRN \'14)',
             'Maternal Nutrition Experiment',
             'Seisure Study',
             'Miscellaneous kwargs["dataset"]s',
             'Unclassified Samples',
             'CEGS2',
             'wild mice',
             'CC Genome',
             'Sister strains',
             'C57BL6 pedigree',
             'CC027 RNA Seq']
MSBWTdirs = ['/csbiodata/CEGS3x3BWT',
             '/csbiodata/BWTs/sangerMsbwt',  # '/playpen/sangerMsbwt',
             '/csbiodata/PompWU-nobackup/msbwts',
             '/csbiodata/BWTs/wild_testis_RNAseq',  # '/playpen/wild_testis_RNAseq',
             '/csbiodata/BWTs/matnut_full',  # '/csbiohome/holtjma/auxiliaryBWTSpace/matnut_full',
             '/csbiodataxv/SeizureStudy',
             '/csbiodata/BWTs/unclassified',  # '/playpen',
             '/csbiodata/BWTs/temporaryWebAccess',  # '/csbiohome/holtjma/auxiliaryBWTSpace/temporaryWebAccess',
             '/csbiodata/BWTs/cegs2/msbwt',  # '/csbiohome/holtjma/auxiliaryBWTSpace/cegs2/msbwt',
             '/csbiodata/BWTs/minibwt',  # '/csbiohome/holtjma/auxiliaryBWTSpace/minibwt',
             '/csbiodata/perlegen/CC_bwt',
             '/csbiodata/perlegen/sisterStrains',
             '/csbiodataxw/C57BL6',
             '/csbiodata/perlegen/CC027_RNASeq']

#Runs multiple kmers against 1 dataset at a time
def runQuery(**kwargs):
    pieces = kwargs["dataset"].split('-')    
    directory = MSBWTdirs[int(pieces[0])] + '/' + '-'.join(pieces[1:])
    # load the MSBWT
    msbwt = MSBWT.loadBWT(directory)
    if kwargs['forward'] == "true":
        forwardResults = [msbwt.countOccurrencesOfSeq(str(kmer)) for kmer in kwargs['kmerQueries']]
    else:
        forwardResults = []
    if kwargs['revComp'] == "true":
        rcResults = [msbwt.countOccurrencesOfSeq(MSBWT.reverseComplement(str(kmer))) for kmer in kwargs['kmerQueries']]
    else:
        rcResults = []
    return [forwardResults, rcResults]

itmp = 0
kmer_list = []
kmer_id_list = []
results_by_kmer = {}

connection = sqlite3.connect('TE_db.sqlite')
connection.row_factory = sqlite3.Row
cursor_read = connection.cursor()
cursor_write = connection.cursor()

#cursor_write.execute('DROP TABLE IF EXISTS kmer_query')
#cursor_write.execute('CREATE TABLE kmer_query (raw_id INTEGER, dataset STRING, forward INTEGER, reverse INTEGER, FOREIGN KEY(raw_id) REFERENCES raw(pid))')

# Loads the sample_list from the text file where each line is a dataset. I then remove datasets 
# that have already been loaded.
start_time = time.time()
sample_list = set()
already_loaded_list = set()
for sample in open('sample_list.txt','r'):
    sample_list.add(sample.replace("\n",""))
for sample in cursor_read.execute('SELECT dataset FROM kmer_query GROUP BY dataset'):
    already_loaded_list.add(sample[0])
print("Samples found already: {} / {}".format( len(already_loaded_list),len(sample_list)))
sample_list = sample_list.difference(already_loaded_list)
print("Loading sample_list took - {0:.2f}s seconds".format(time.time() - start_time))
start_time = time.time()

# I will now load the list of kmers we are running the query with from the database.
kmer_list = []
pid_list = []
for row in cursor_read.execute('SELECT * FROM raw WHERE TE = 1'):
    if row['ref_like_prefix'] != '' and row['insertion_after_prefix'] != '':
        kmer = row['ref_like_prefix'] + row['insertion_after_prefix']
    elif row['insertion_before_suffix'] != '' and row['ref_like_suffix'] != '':
        kmer = row['insertion_before_suffix'] + row['ref_like_suffix']    
    elif row['ref_like_prefix'] != '' and row['ref_like_suffix'] != '':
        kmer = row['ref_like_prefix'] + row['ref_like_suffix']
    else:
        kmer = ''

    if kmer != '':
        kmer_list.append(kmer.upper())
        pid_list.append(row['pid'])
print("Loading kmer_list took - {0:.2f}s for {1} entries".format((time.time() - start_time),len(kmer_list)))
start_time = time.time()

# For each sample run all of the kmers against it and write the results to the database. The write operation
# should be atomic so that no partial information is written for a sample.
for sample in sample_list:
    itmp = itmp + 1
    results = runQuery(dataset=sample.replace("\n",""),
                   kmerQueries=kmer_list,
                   forward='true',
                   revComp='true')
    for i in xrange(1,len(results[0])):
        if not(results[0][i] == 0 and results[1][i] == 0):
            cursor_write.execute('INSERT INTO kmer_query(raw_id, dataset, forward, reverse) VALUES (?,?,?,?)', (pid_list[i], sample.replace("\n",""), results[0][i], results[1][i]))
    # Commit after each sample
    connection.commit()
    print("Returning results from sample query for sample {0} took - {1:.2f}s".format(sample,time.time() - start_time))    
    start_time = time.time()


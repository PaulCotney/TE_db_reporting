import sqlite3


def findstartandstopforsample(self, rowon):
	self.samplestart = rowon.keys().index("genotype")

def loadrow(rowon):
	te_general = {}

connection = sqlite3.connect('TE_db.sqlite')
connection.row_factory = sqlite3.Row
cursor = connection.cursor()
first_row = {}
second_row = {}
for row in cursor.execute('SELECT * FROM load_mappable_tes'):
	findstartandstopforsample(row)
	print (type(row))
	print(row.keys())
	break



import sqlite3
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('gene')
args = parser.parse_args()

conn = sqlite3.connect('protein2ipr-paired.db')
c = conn.cursor()

output = open('../tmp/domain_out.txt', 'w')

row1 = c.execute('SELECT item_id_b FROM protein_actions WHERE item_id_a=?', (args.gene,))
#row2 = c.execute('SELECT start FROM protein_actions WHERE item_id_a=?', (args.gene,))
#row3 = c.execute('SELECT end FROM protein_actions WHERE item_id_a=?', (args.gene,))

data1 = []
data2 = []
data3 = []

for i in row1:
	data1.append(i[0])

row2 = c.execute('SELECT start FROM protein_actions WHERE item_id_a=?', (args.gene,))

for i in row2:
	data2.append(i[0])

row3 = c.execute('SELECT end FROM protein_actions WHERE item_id_a=?', (args.gene,))

for i in row3:
	data3.append(i[0])

#print data1

#for row in c.execute('SELECT end FROM protein_actions WHERE item_id_a=?', (args.gene,)):
for i in range(len(data1)):
    #print row[0][1]
    #out_list = list(row[0])
    #for item in row[0]:
    output.write("%s\n" % data1[i])
    #print ("%s" % data1[i])

for row in c.execute('SELECT item_id_a, start, end FROM protein_actions WHERE item_id_b=?', (args.gene,)):
    #print row[0][1]
    #output.write(row[0])
    #out_list = list(row[0])
    #for item in row[0]:
    output.write("%s\n" % row[0])

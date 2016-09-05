import sqlite3
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('gene')
args = parser.parse_args()

conn = sqlite3.connect('pdbtosp_paired.db')
c = conn.cursor()

output = open('../tmp/uni_id_out.txt', 'w')

for row in c.execute('SELECT item_id_b FROM protein_actions WHERE item_id_a=?', (args.gene,)):
    #print row[0][1]
    #out_list = list(row[0])
    #for item in row[0]:
    output.write("%s\n" % row[0])

for row in c.execute('SELECT item_id_a FROM protein_actions WHERE item_id_b=?', (args.gene,)):
    #print row[0][1]
    #output.write(row[0])
    #out_list = list(row[0])
    #for item in row[0]:
    output.write("%s\n" % row[0])

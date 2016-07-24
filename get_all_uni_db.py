import sqlite3
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('gene')
args = parser.parse_args()

conn = sqlite3.connect('all_go_knowledge_full_paired.db')
c = conn.cursor()

output = open('../tmp/all_out.txt', 'w')

for row in c.execute('SELECT item_id_b FROM protein_actions WHERE item_id_a=?', (args.gene,)):
    #print row[0]
    output.write("%s\n" % row[0])

for row in c.execute('SELECT item_id_a FROM protein_actions WHERE item_id_b=?', (args.gene,)):
    #print row[0]
    output.write("%s\n" % row[0])

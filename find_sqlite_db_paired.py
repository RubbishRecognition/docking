import sqlite3
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('gene')
args = parser.parse_args()

output = open('../tmp/' + str(args.gene), 'w')

conn = sqlite3.connect('protein.actions.v10_paired.db')
c = conn.cursor()

for row in c.execute('SELECT item_id_b FROM protein_actions WHERE item_id_a=?', (args.gene,)):
    #print row[0]
    output.write("%s\n" % row[0])

for row in c.execute('SELECT item_id_a FROM protein_actions WHERE item_id_b=?', (args.gene,)):
    #print row[0]
    output.write("%s\n" % row[0])

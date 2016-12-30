import re
import sqlite3

conn = sqlite3.connect('pdbtosp_paired.db')
c = conn.cursor()

c.execute('pragma main.page_size=4096;')
c.execute('pragma main.synchronous=0;')
conn.commit()
c.execute('CREATE TABLE IF NOT EXISTS protein_actions (item_id_a TEXT, item_id_b TEXT)')
c.execute('DELETE FROM protein_actions')
conn.commit()

interactions = []

with open('../id_mapping/_pdbtosp.txt') as f:
        f.readline()
        for line in f:

		l = re.sub (",", " ", line)
		l = (' '.join(l[28:].split())).split(' ')

		int1 = []
		int2 = []

		n = 0
		for curr in l:
			if (n % 2 == 0):
				int1.append(curr)
			else:
				print curr
				if (curr != 'see'):
					int2.append(curr.split('(')[1].split(')')[0])	
			n += 1

		n = 0
		#print int1
		#print int2
		for n in range(len(int1)):
			if (len(int1) == len(int2)):
				interactions.append( (int1[n], int2[n]) )
                #interaction = line.split('\t')
                #interaction1 = interaction[0] + '_' + interaction[1]
                #interaction2 = interaction[5]

                #interactions.append( (interaction1, interaction2) )

c.executemany('INSERT INTO protein_actions VALUES (?, ?)', interactions)

conn.commit()
c.execute('CREATE INDEX IF NOT EXISTS indexa ON protein_actions (item_id_a)')
c.execute('CREATE INDEX IF NOT EXISTS indexb ON protein_actions (item_id_b)')
conn.commit()
c.execute('pragma main.synchronous=2;')
conn.commit()
conn.close()

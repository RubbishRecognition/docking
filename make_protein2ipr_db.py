import sqlite3

conn = sqlite3.connect('protein2ipr-paired.db')
c = conn.cursor()

c.execute('pragma main.page_size=4096;')
c.execute('pragma main.synchronous=0;')
conn.commit()
c.execute('CREATE TABLE IF NOT EXISTS protein_actions (item_id_a TEXT, item_id_b TEXT, start TEXT, end TEXT)')
c.execute('DELETE FROM protein_actions')
conn.commit()

#interactions = []

with open('protein2ipr.dat') as f:
        f.readline()
        for line in f:
                interaction = line.split('\t')
                #interaction1 = interaction[1].split('|')[1].split('_')[0]
                #interaction2 = interaction[2]
		interactions = []
                interactions.append( (interaction[0], interaction[1], interaction[4], interaction[5].split('\n')[0]) )
		#print interactions
		#k = input()
		c.executemany('INSERT INTO protein_actions VALUES (?, ?, ?, ?)', interactions)
#c.executemany('INSERT INTO protein_actions VALUES (?, ?, ?, ?)', interactions)

conn.commit()
		
c.execute('CREATE INDEX IF NOT EXISTS indexa ON protein_actions (item_id_a)')
c.execute('CREATE INDEX IF NOT EXISTS indexb ON protein_actions (item_id_b)')
c.execute('CREATE INDEX IF NOT EXISTS indexb ON protein_actions (start)')
c.execute('CREATE INDEX IF NOT EXISTS indexb ON protein_actions (end)')
conn.commit()
c.execute('pragma main.synchronous=4;')
conn.commit()
conn.close()

#k = input()

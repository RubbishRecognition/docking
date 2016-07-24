import sqlite3

conn = sqlite3.connect('all_go_knowledge_full_paired.db')
c = conn.cursor()

c.execute('pragma main.page_size=4096;')
c.execute('pragma main.synchronous=0;')
conn.commit()
c.execute('CREATE TABLE IF NOT EXISTS protein_actions (item_id_a TEXT, item_id_b TEXT)')
c.execute('DELETE FROM protein_actions')
conn.commit()

interactions = []

with open('../id_mapping/all_go_knowledge_full.tsv') as f:
	f.readline()
	for line in f:
		interaction = line.split('\t')
		interaction1 = interaction[1]
		interaction2 = interaction[2]
		
		interactions.append( (interaction1, interaction2) )	

c.executemany('INSERT INTO protein_actions VALUES (?, ?)', interactions)

conn.commit()
c.execute('CREATE INDEX IF NOT EXISTS indexa ON protein_actions (item_id_a)')
c.execute('CREATE INDEX IF NOT EXISTS indexb ON protein_actions (item_id_b)')
conn.commit()
c.execute('pragma main.synchronous=2;')
conn.commit()
conn.close()

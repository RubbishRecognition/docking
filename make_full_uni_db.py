import sqlite3

conn = sqlite3.connect('full_uniprot_2_string.04_2015_paired.db')
c = conn.cursor()

c.execute('pragma main.page_size=4096;')
c.execute('pragma main.synchronous=0;')
conn.commit()
c.execute('CREATE TABLE IF NOT EXISTS protein_actions (item_id_a TEXT, item_id_b TEXT)')
c.execute('DELETE FROM protein_actions')
conn.commit()

interactions = []

with open('../id_mapping/full_uniprot_2_string.04_2015.tsv') as f:
        f.readline()
        for line in f:
                interaction = line.split('\t')
                interaction1 = interaction[1].split('|')[1]#.split('_')[0]
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

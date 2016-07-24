import sqlite3

conn = sqlite3.connect('protein.actions.v10_paired.db')
c = conn.cursor()

c.execute('pragma main.page_size=4096;')
c.execute('pragma main.synchronous=0;')
conn.commit()
c.execute('CREATE TABLE IF NOT EXISTS protein_actions (item_id_a TEXT, item_id_b TEXT)')
c.execute('DELETE FROM protein_actions')
conn.commit()

interactions = []

with open('../interaction_data/protein.actions.v10.txt') as f:
    f.readline()
    for line in f:
        interaction = line.split('\t')
        gene_start = interaction[0].find('.')+1
        interaction1 = interaction[0][gene_start:]
        gene_start = interaction[1].find('.')+1
        interaction2 = interaction[1][gene_start:]

        #interactions are repeated in the file as both a,b and b,a; let's #keep only
        #the copy where a<=b to save some space
        if interaction1 <= interaction2:
            interactions.append( (interaction1, interaction2) )

        #reasonable number of interactions to store in memory at once, delete a zero
        #if it is too many for your computer
        if len(interactions) > 10000000:
            c.executemany('INSERT INTO protein_actions VALUES (?, ?)', interactions)
            interactions = []

c.executemany('INSERT INTO protein_actions VALUES (?, ?)', interactions)

conn.commit()
c.execute('CREATE INDEX IF NOT EXISTS indexa ON protein_actions (item_id_a)')
c.execute('CREATE INDEX IF NOT EXISTS indexb ON protein_actions (item_id_b)')
conn.commit()
c.execute('pragma main.synchronous=2;')
conn.commit()
conn.close()

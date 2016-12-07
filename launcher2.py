import os

queries = open('queries.txt', 'r')

for q in queries:
        print q
        os.system(q)

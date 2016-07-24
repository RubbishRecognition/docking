import os

x = open ('files.txt', 'r')

dbg = False

if dbg:
	for query in x:
		print query
	k = input()

for query in x:
	#print query
	y = query.translate(None, '\n')
	y = y.lower()
	y1 = y.split(' ')[0]
	y2 = y.split(' ')[1].translate(None, '\n')
	y1 = y1.split('_')[0]
	y2 = y2.split('_')[0]
	#string = ('python grid.py ' + query.translate(None, '\n') + ' out_' + str(query.translate(None, '.fasta.txt')).translate('_', ' ') + '.txt')
	string = 'python grid.py ' + y1 + '.fasta.txt ' + y2 + '.fasta.txt ' + 'out_biogrid_' + y1 + '_' + y2 + '.txt ' + 'out_string_0_' + y1 + '_' + y2 + '.txt'
	print string
	#k = input()
	os.system (string)

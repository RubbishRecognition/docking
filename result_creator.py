import numpy as np
import pandas as pd
from os import listdir
from string import maketrans
import argparse
import ConfigParser
from id_converters import *

def check (a_s, a_e, b_s, b_e):
	res = 0

	if (a_s <= b_s) & (b_e <= a_e):
		res = 1
	if (b_s <= a_s) & (a_e <= b_e):
		res = 1
	if (a_s <= b_s) & (a_e <= b_e) & (b_s <= a_e):
		res = 1
	if (b_s <= a_s) & (b_e <= a_e) & (a_s <= b_e):
		res = 1

	return res

def main ():
	parser = argparse.ArgumentParser (description='Result maker')
	parser.add_argument ('lig', metavar='lig', help='ligand file')
	parser.add_argument ('rec', metavar='rec', help='receptor file')
	parser.add_argument ('z_lig', metavar='z_lig', help='zlab_ligand file')
	parser.add_argument ('z_rec', metavar='z_rec', help='zlab_receptor file')
	parser.add_argument ('out', metavar='out', help='output file')
	parser.add_argument ('prot1', metavar='prot1', help='sequence1')
	parser.add_argument ('prot2', metavar='prot2', help='sequence2')	

	args = parser.parse_args ()

	ligand = np.array(pd.read_csv(args.lig, header=0))
	receptor = np.array(pd.read_csv(args.rec, header=0))
	zlab_ligand = np.array(pd.read_csv(args.z_lig, header=0))
	zlab_receptor = np.array(pd.read_csv(args.z_rec, header=0))	
	output = open('buf.csv', 'w')
	
	prot = args.lig.split('_')[3] + '_' + args.lig.split('_')[4]
	
	len1 = len(open (args.prot1, 'r').readlines()[0])
        len2 = len(open (args.prot2, 'r').readlines()[0])
	
	thebest = []
	rank = 0
	bestlen = 0
	trantab = maketrans('\'',' ')
	for x in zlab_ligand:
		a = str(x[6])
		b = str(x[7])
		b = b.translate(trantab)
		a = int(a)
		b = int(b)
		if ((x[4] != 'None') & (x[4] != '\'0') & ((b - a) >= 60)):
			if (rank < int(x[3])):
				rank = int(x[3])
				thebest = x
				bestlen = b - a
	#print files[i]
	#print thebest

	true = [1 for j in range(len(ligand))]
	used = [0 for j in range(len(ligand))]

	for j in range(len(ligand)):
		if used[j] == 0:
			used[j] = 1
			for k in range(len(ligand)):
				if used[k] == 0:
					if (check(ligand[j][3], ligand[j][4], ligand[k][3], ligand[k][4])):
						used[k] == 1
						true[k] = 0

	#print true

	score = 1
	for j in range(len(ligand)):
		x = ligand[j]
		if (x[4] - x[3] >= 60) & ((x[1] + 1) > rank) & (true[j] == 1):
			score += 1
	if (len(thebest) > 0):
		print thebest[0] + ',ligand,' + uniprot_pdb_chain_converter(prot) + ',' + thebest[1] + ',' + str(score) + ',' + str(thebest[2]) + ',' + str(thebest[3]) + ',' + thebest[4].split('\'')[1]  + ',' + str(bestlen)
		output.write(thebest[0] + ',ligand,' + uniprot_pdb_chain_converter(prot) + ',' + thebest[1] + ',' + str(score) + ',' + str(thebest[2]) + ',' + str(thebest[3]) + ',' + thebest[4].split('\'')[1]  + ',' + str(bestlen) + ',' + str(len1) + '\n')
	else:
		print None
		output.write('None\n')
	#print score

	## Receptors

	prot = args.lig.split('_')[1] + '_' + args.lig.split('_')[2]

	thebest = []
	rank = 0
	bestlen = 0
	trantab = maketrans('\'',' ')
	for x in zlab_receptor:
		a = str(x[6])
		b = str(x[7])
		b = b.translate(trantab)
		a = int(a)
		b = int(b)
		if ((x[4] != 'None') & (x[4] != '\'0') & ((b - a) >= 60)):
			if (rank < int(x[3])):
				rank = int(x[3])
				thebest = x
				bestlen = b - a
	#print files[i]
	#print thebest

	true = [1 for j in range(len(receptor))]
	used = [0 for j in range(len(receptor))]

	for j in range(len(receptor)):
		if used[j] == 0:
			used[j] = 1
			for k in range(len(receptor)):
				if used[k] == 0:
					if (check(receptor[j][3], receptor[j][4], receptor[k][3], receptor[k][4])):
						used[k] == 1
						true[k] = 0

	#print true

	score = 1
	for j in range(len(receptor)):
		x = receptor[j]
		if (x[4] - x[3] >= 60) & ((x[1] + 1) > rank) & (true[j] == 1):
			score += 1
	#print score
	if (len(thebest) > 0):
		print thebest[0] + ',receptor,' + uniprot_pdb_chain_converter(prot) + ',' + thebest[1] + ',' + str(score) + ',' + str(thebest[2]) + ',' + str(thebest[3]) + ',' + thebest[4].split('\'')[1]  + ',' + str(bestlen)
		output.write(thebest[0] + ',receptor,' + uniprot_pdb_chain_converter(prot) + ',' + thebest[1] + ',' + str(score) + ',' + str(thebest[2]) + ',' + str(thebest[3]) + ',' + thebest[4].split('\'')[1]  + ',' + str(bestlen) + ',' + str(len2) + '\n')
	else:
		print None
		output.write('None\n')

	output.close()
	
	uniprot = np.array(pd.read_csv('buf.csv', sep=',', header = None))
	output = open (args.out, 'w')

	output.write('complex\tchain\tuniprot_id\tdomain_id\trank\tall_matches\tn_matches\tn_int_atoms\tdomain_len\tfull_seq_len\n')
	for x in uniprot:

        	line = [0 for i in range(2000)]
        	#domains = get_domain_id_list(x[2])
        	aim = x[3]
        	prot = x[2]
        	num = 0
        	domains = get_domain_id_list(x[2])
        	#print domains
        	for rec in domains:
                	if aim == rec[0].split(' ')[0]:
                        	i = int(rec[0].split(' ')[1])
                        	j = int(rec[0].split(' ')[2])
                        	for k in range(j - i):
                                	line[k + i] = 1
        	output.write(str(x[0]) + '\t' + str(x[1]) + '\t' + x[2] + '\t' + str(x[3]) + '\t' + str(x[4])  + '\t' + str(x[5]) + '\t' + str(x[6]) + '\t' + str(x[7]) + '\t' + str(sum(line)) + '\t' + str(x[8]) + '\n')
	
	os.system('rm buf.csv')


if __name__ == '__main__':
        main()

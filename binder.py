import os
import argparse
import numpy as np
import pandas as pd
from Bio import SeqIO
from multiprocessing import Pool
from Bio import SeqIO
import ConfigParser
from Bio.PDB import *
from domain_work import *
from get_ids_from_msa_files import *
from id_converters import *
from interaction_search import *
from script_running import *

#pdbtosp location
pdbtosp = '../id_mapping/pdbtosp.txt'

#BIOGRID database location
biogrid = '../interaction_data/BIOGRID-ALL-3.4.136.mitab.txt'

#Uniprot database location
uniprot = '../interaction_data/uniprot_sprot_varsplic.fasta'

#String database location
string = '../interaction_data/protein.actions.v10.txt'

#Uniprot20 database location
uniprot20 = '../seq_data/uniprot20_2016_02'

#pdb100 database location
pdb100 = '../seq_data/pdb100'

full_uniprot = ''

all_knowledge = ''

#config = open ('./config', 'r')
#for x in config.readlines():
#x = config.readlines()

def ConfigSectionMap(section):
    dict1 = {}
    options = Config.options(section)
    for option in options:
        try:
            dict1[option] = Config.get(section, option)
            if dict1[option] == -1:
                DebugPrint("skip: %s" % option)
        except:
            print("exception on %s!" % option)
            dict1[option] = None
    return dict1

Config = ConfigParser.ConfigParser()
Config.read("./config")

pdbtosp = ConfigSectionMap("PATHS")['pdbtosp']#x[0].split('|')[1].split('\n')[0]
biogrid = ConfigSectionMap("PATHS")['biogrid']#x[1].split('|')[1].split('\n')[0]
uniprot = ConfigSectionMap("PATHS")['uniprot']#x[2].split('|')[1].split('\n')[0]
string  = ConfigSectionMap("PATHS")['string']#x[3].split('|')[1].split('\n')[0]
uniprot20 = ConfigSectionMap("PATHS")['uniprot20']#x[4].split('|')[1].split('\n')[0]
pdb100 = ConfigSectionMap("PATHS")['pdb100']#x[5].split('|')[1].split('\n')[0]
full_uniprot = ConfigSectionMap("PATHS")['full_uniprot']#x[6].split('|')[1].split('\n')[0]
all_knowledge = ConfigSectionMap("PATHS")['all_knowledge']#x[7].split('|')[1].split('\n')[0]
interpro = ConfigSectionMap("PATHS")['interpro']

score = int(ConfigSectionMap("hhblits options")['score'])
cpu = int(ConfigSectionMap("hhblits options")['cpu'])
overlap = int (ConfigSectionMap("hhblits options")['overlap'])

def seq_test (seq, seq_len):

	## Right sequence - if each MSA block contains >= 60 letters

        count = 0
        flag = 0
        for i in seq:
                if (i == '-'):
                        if flag:
                                if (count < seq_len) & (count > 0):
                                        return 0
                        count = 0
                else:
                        flag = 1
                        count += 1
        if ((count > 0) & (count < seq_len)):
                return 0
        else:
                return 1



def config_parser (conf):
	
	## Config parser
	
	config = open (conf, 'r')
	#for x in config.readlines():
	x = config.readlines()

	pdbtosp = x[0].split('|')[1]
	biogrid = x[1].split('|')[1]			
	uniprot = x[2].split('|')[1]
	string  = x[3].split('|')[1]
	uniprot20 = x[4].split('|')[1]
	pdb100 = x[5].split('|')[1]
	full_uniprot = x[6].split('|')[1]
	all_knowledge = x[7].split('|')[1]



def find_align_pos (seq):

	## Find start and end pos
	
	start = 0
	end = 0
	flag_s = 1
	flag_e = 1
	for i in range(len(seq)):
		if ((seq[i] != '-') & (flag_s)):
			start = i
			flag_s = 0
		if ((seq[len(seq) - i - 1] != '-') & (flag_e)):
			end = len(seq) - i - 1
			flag_e = 0
	return start, end


def main ():

	#print "START"

	## Input sequences as a command line args
	
	parser = argparse.ArgumentParser (description='Finding protein interactions')
	parser.add_argument ('protein_A', metavar='prot_A',  help='file with sequence A')
	parser.add_argument ('protein_B', metavar='prot_B',  help='file with sequence B')
        parser.add_argument ('fout_biogrid', metavar='pdb_out',  help='output file for pdb interactions')
        parser.add_argument ('fout_string', metavar='uniprot_out',  help='output file for uniprot interactions')
	parser.add_argument ('fout_complex', metavar='complex_out',  help='output file for complexes')
	parser.add_argument ('fout_align_prot_A', metavar='align_out',  help='output file for multiple sequence alignment for sequence A')
	parser.add_argument ('fout_hhr_prot_A', metavar='hhr_out',  help='output file for hhalign result for sequence A')
	parser.add_argument ('fout_align_prot_B', metavar='align_out',  help='output file for multiple sequence alignment for sequence B')
        parser.add_argument ('fout_hhr_prot_B', metavar='hhr_out',  help='output file for hhalign result for sequence B')
	parser.add_argument ('fout_domain_A', metavar='domain_A',  help='output file for protein_A domain')
	parser.add_argument ('fout_domain_B', metavar='domain_B',  help='output file for protein_B domain')
	parser.add_argument ('fout_distance', metavar='distance_out',  help='output file for domains_distances')		

        args = parser.parse_args ()
	
	#Input reformat to original FASTA
	prot1 = open(args.protein_A + '.formatted', 'w')
	prot2 = open(args.protein_B + '.formatted', 'w')

	seq1 = open(args.protein_A, 'r').readlines()
	seq2 = open(args.protein_B, 'r').readlines()
		
	prot1.write('>1BQU:A|PDBID|CHAIN|SEQUENCE\n' + seq1[0])
	prot2.write('>1BQU:A|PDBID|CHAIN|SEQUENCE\n' + seq2[0])

	prot1.close()
	prot2.close()	

	print "START HHBLITS"

	run_hhblits (uniprot20, pdb100, cpu, score, args.protein_A + '.formatted', args.fout_align_prot_A, args.fout_hhr_prot_A)
	run_hhblits (uniprot20, pdb100, cpu, score, args.protein_B + '.formatted', args.fout_align_prot_B, args.fout_hhr_prot_B)
	
	print 'START ID MAPPING'
	prime_A, pdb_ids_A = get_pdb_ids (args.fout_align_prot_A, overlap)
	biogrid_ids_A = get_biogrid_chain_ids (args.fout_align_prot_A, overlap)
	string_ids_A = get_string_ids_sql (args.fout_align_prot_A, overlap)

	prime_B, pdb_ids_B = get_pdb_ids (args.fout_align_prot_B, overlap)
	biogrid_ids_B = get_biogrid_chain_ids (args.fout_align_prot_B, overlap)
        string_ids_B = get_string_ids_sql (args.fout_align_prot_B, overlap)

	## TO KEEP
	to_keep_a = biogrid_ids_A
	to_keep_b = biogrid_ids_B
	## Find interactions

	if (len (biogrid_ids_A) > len (biogrid_ids_B)):
		biogrid_interactions = find_biogrid_interactions (biogrid_ids_B, biogrid_ids_A, cpu, biogrid)
		string_interactions = find_string_interactions (string_ids_B, string_ids_A, cpu, string)
	else:
		biogrid_interactions = find_biogrid_interactions (biogrid_ids_A, biogrid_ids_B, cpu, biogrid)
                string_interactions = find_string_interactions (string_ids_A, string_ids_B, cpu, string)

	pdb_interactions = find_pdb_id_sql (biogrid_interactions, 'B') + find_pdb_id_sql (string_interactions, 'S')
	uniprot_interactions = biogrid_interactions + find_uniprot_id_sql (string_interactions, 'S')

	## Write down the results

	output = open (args.fout_biogrid, 'w')
	for x in pdb_interactions:
		output.write(str(x))

	output = open (args.fout_string, 'w')
	for x in uniprot_interactions:
                x_a = x.split(' ')[0]
		x_b = x.split(' ')[1].split('\n')[0]
		flag = 1
		if ((x_a in to_keep_a) & (x_b in to_keep_b)):
			output.write(x_a + ' ' + x_b + '\n')
		if ((x_a in to_keep_b) & (x_b in to_keep_a)):
			output.write(x_b + ' ' + x_a + '\n')
	#Finding complexes
	complexes = []
	test = list()
	for x in pdb_interactions:
		if (x.split(' ')[0].split('_')[0] == x.split(' ')[1].split('\n')[0].split('_')[0]):
			complexes.append(x.split(' ')[0] + ':' + x.split(' ')[1].split('\n')[0].split('_')[1] + '\n')
			test.append([(x.split(' ')[0] + ':' + x.split(' ')[1].split('\n')[0].split('_')[1] + '\n')])


	output = open (args.fout_complex, 'w')
	for x in test:
		for y in x:
			output.write(str(y))
	output.close()
	
	find_domains(args.fout_string, args.fout_domain_A, args.fout_domain_B)	
	domains_A_co, domains_B_co, domains_co = find_interacting_domains (interpro, args.fout_align_prot_A, args.fout_align_prot_B, args.fout_complex, args.protein_A + '.formatted', args.protein_B + '.formatted', pdb_ids_A, pdb_ids_B, args.fout_distance)	
	domains_A_contra, domains_B_contra, domains_contra = find_interacting_domains (interpro, args.fout_align_prot_B, args.fout_align_prot_A, args.fout_complex, args.protein_B + '.formatted', args.protein_A + '.formatted', pdb_ids_B, pdb_ids_A, args.fout_distance)
	dist_A, dist_B = prime_domains (args.fout_domain_A, args.fout_domain_B, domains_co, domains_contra, args.fout_distance)
	output_binder1 (interpro, prime_A, prime_B, args.fout_distance, domains_A_co, domains_A_contra, dist_A, dist_B, args.protein_A + '.formatted', args.protein_B + '.formatted', args.fout_domain_A, args.fout_domain_B)	

if __name__ == '__main__':
        main()


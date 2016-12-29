#from domain_work import *
#from get_ids_from_msa_files import *
from id_converters import *
#from interaction_search import *
#from script_running import *
#from binder import *
import os
import numpy as np
import pandas as pd
from Bio import SeqIO
from multiprocessing import Pool
from Bio import SeqIO
import ConfigParser
from Bio.PDB import *

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

##	Fantastic domains and where to find them
##	Output files will contain records in format: 'Domain: Uniprot ID's contains it'
##	Uniprot ID's are given from output file with uniprot interactions 

def find_domains (infile, outfile_A, outfile_B):

        ## Find domains from interaction output file with uniprot ID's

        pdb_a = []
        pdb_b = []

        pdb = open (infile, 'r')
        for i in pdb.readlines ():
                pdb_a.append (i.split (' ')[0])#.split('_')[0])
                pdb_b.append (i.split (' ')[1].split ('\n')[0])#.split('_')[0])

        #print pdb_a
        #print pdb_b

        uniprot_a = pdb_a
        uniprot_b = pdb_b

        #print uniprot_a
        #print uniprot_b

        for i in range(len(uniprot_a)):
                res = set()
                os.system ('python get_uni_id.py ' + str(uniprot_a[i]))
                tmp = open('../tmp/uni_id_out.txt', 'r')
                for j in tmp.readlines():
                        res.update([j.split('\n')[0]])
                if (len(list(res)) == 1):
                        uniprot_a[i] = list(res)[0]
                if (len(list(res)) > 1):
                        print 'Shit happened' + ' ' + uniprot_a[i]
                        print res
                if (len(list(res)) == 0):
                        uniprot_a[i] = uniprot_a[i].split('_')[0]

        for i in range(len(uniprot_b)):
                res = set()
                os.system ('python get_uni_id.py ' + str(uniprot_b[i]))
                tmp = open('../tmp/uni_id_out.txt', 'r')
                for j in tmp.readlines():
                        res.update([j.split('\n')[0]])
                if (len(list(res)) == 1):
                        uniprot_b[i] = list(res)[0]
                if (len(list(res)) > 1):
                        print 'Shit happened' + uniprot_b[i]
                if (len(list(res)) == 0):
                        uniprot_b[i] = uniprot_b[i].split('_')[0]

        #print uniprot_a
        #print uniprot_b

        domain_a = dict ()
        domain_b = dict ()

        flag = 1
        for i in range (len (uniprot_a)):
                curr = get_domain_id (uniprot_a[i])
                for rec in curr:
                        if (domain_a.get(rec) != None):
                                domain_a[rec] += [uniprot_a[i]]
                        else:
                                domain_a[rec] = [uniprot_a[i]]

                curr = get_domain_id (uniprot_b[i])
                for rec in curr:
                        if (domain_b.get(rec) != None):
                                domain_b[rec] += [uniprot_b[i]]
                        else:
                                domain_b[rec] = [uniprot_b[i]]

        out_A = open(outfile_A, 'w')
        out_B = open(outfile_B, 'w')
        #print(domain_a)
        #print(domain_b)
        for i in domain_a.items():
                out_A.write(i[0] + ' : ')
                for j in set(i[1]):
                        out_A.write(j + ' ')
                out_A.write('\n')
        for i in domain_b.items():
                out_B.write(i[0] + ' : ')
                for j in set(i[1]):
                        out_B.write(j + ' ')
                out_B.write('\n')


##	Caculates distaces between each domain in domain_list and chain_B
##	domain_list - list of domain ID's from chain_A sequence
##	structure1 and structure2 - pdb structures of chain_A * chain_B complex
##	chain_A and chain_B chain letters like 'A' or 'R'
##	Returns list, each record consist Domain ID, number of it's good atoms, all number of it's atoms, start and end positions in pdb strucure

def distance_calculator (domain_list, structure1, structure2, chain_A, chain_B):

	## Distance calculator from domains of chain_A to chain_B
	distance_list = []
	for model in structure1:
                chain = model[chain_A]        
		for domain in domain_list:
			place = 0
			distance = 1.0 * float(0)
			num = 0
			int_num = 0
			for residue in chain:
				if ((place >= domain[1]) & (place <= domain[2])):
					for atom in residue:
			         		min_dist = 99999999
						for model2 in structure2:
                                                        chain2 = model2[chain_B]
							for residue2 in chain2:
								for atom2 in residue2:
									if (abs(atom - atom2) < min_dist):
										min_dist = abs(atom - atom2)  		
						distance += min_dist
						## Means that atom exist in interation area
						if (min_dist < 5):
							int_num += 1
						num += 1
				place += 1
			distance = str(int_num) +  ' ' + str(num) + ' ' + str(domain[1]) + ' ' + str(domain[2])
			distance_list.append(domain + [distance])
	return distance_list	

##	Gets distance calculater output for each complex from complex output file
##	align_A, aling_B - pdb_msa output files, to find pdb sequence of Prime domains
##	complexes - complex output file, consist list of complex structures
##	prime_A, prime_B - PDB ID's of Prime proteins in format '1g0y_R'
##	pdb_ids_A, pdb_ids_B - lists of PDB ID's that identical to Prime_A and Prime_B
##	returns 3 lists - domains in Prime_A seq, domains in Prime_B seq, list of distances between Prime_A domains and Prime_B for each complex (same as distance calculater output)

def find_interacting_domains (interpro, align_A, align_B, complexes, prime_A, prime_B, pdb_ids_A, pdb_ids_B, outfile):

	## Calculate distances and so on
	result = []
	domains_A = get_domains_list (prime_A, interpro)
	domains_B = get_domains_list (prime_B, interpro)

	complex_list = open(complexes, 'r').readlines()
	complex_positions = list()
	for comp in complex_list:
		chain_A = comp.split(':')[0]
		chain_B = comp.split('_')[0] + '_' + comp.split(':')[1].split('\n')[0]
		comp_domains_A = list()
		comp_domains_B = list()
		if ((chain_A in pdb_ids_A) & (chain_B in pdb_ids_B) & (chain_A != chain_B)):
			## Find pdb struct for chain_A	
			handle = open('pdb_' + align_A, 'r')
			for record in SeqIO.parse(handle, "fasta"):
				if (len(record.id) > 10):
					record.id = record.id.split(':')[0].lower() + '_' + record.id.split(':')[1].split('|')[0]
				if (record.id == chain_A):
					start, end = find_align_pos (record.seq)	
					for domain in domains_A:
						if ((start <= int(domain[1])) & (end >= int(domain[2]))):
							comp_domains_A.append([domain[0], int(domain[1]) - start, int(domain[2]) - start])
			## Find pdb struct for chain_B
			handle = open('pdb_' + align_B, 'r')
                        for record in SeqIO.parse(handle, "fasta"):
				if (len(record.id) > 10):
                                        record.id = record.id.split(':')[0].lower() + '_' + record.id.split(':')[1].split('|')[0]
				if (record.id == chain_B):
                                        start, end = find_align_pos (record.seq)
                                        for domain in domains_B:
                                                if ((start <= int(domain[1])) & (end >= int(domain[2]))):
                                                        comp_domains_B.append([domain[0], int(domain[1]) - start, int(domain[2]) - start])	 

			## Distance calculating
			pdbl = PDBList()
			pdbl.retrieve_pdb_file(comp.split('_')[0].upper(), pdir = '../tmp/')
			parser = PDBParser()
			structure1 = parser.get_structure('X', '../tmp/pdb' + comp.split('_')[0] + '.ent')	
			structure2 = parser.get_structure('X', '../tmp/pdb' + comp.split('_')[0] + '.ent')
			distance_list_A = distance_calculator (comp_domains_B, structure1, structure2, chain_B.split('_')[1], chain_A.split('_')[1])
			structure1 = parser.get_structure('X', '../tmp/pdb' + comp.split('_')[0] + '.ent')
                        structure2 = parser.get_structure('X', '../tmp/pdb' + comp.split('_')[0] + '.ent')
			distance_list_B = distance_calculator (comp_domains_A, structure1, structure2, chain_A.split('_')[1], chain_B.split('_')[1])	
			dist_dict = dict()
			for i in distance_list_A:
				buf = dist_dict.get(i[0])
				if (buf == None):
					dist_dict[i[0]] = [i[3], i[2] - i[1]]
				else:
					if (float(i[3].split(' ')[0]) > float(buf[0].split(' ')[0])):
						dist_dict[i[0]] = [i[3], i[2] - i[1]]
			for i in distance_list_B:
                                buf = dist_dict.get(i[0])
                                if (buf == None):
					dist_dict[i[0]] = [i[3], i[2] - i[1]]
                                else:
                                        if (float(i[3].split(' ')[0]) > float(buf[0].split(' ')[0])):
						dist_dict[i[0]] = [i[3], i[2] - i[1]]		
			result.append([comp, dist_dict])

			output = open(outfile, 'a')
			output.write(comp + '\n')
			for i in distance_list_A:
				output.write(i[0] + ' ' + str(i[1]) + ' ' + str(i[2]) + ' ' + str(i[3]) + '\n')
			for i in distance_list_B:
                                output.write(i[0] + ' ' + str(i[1]) + ' ' + str(i[2]) + ' ' + str(i[3]) + '\n')
			output.write('\n\n')			
	return domains_A, domains_B, result


##	Creates output file, just for testing
##	domain_A, domain_B - domain* output files 
##	domains_co, domains_comtra - find_interacting_domains output, co for straight input (prime_A, prime_B), contra for reverce (prime_B, prime_A)

def prime_domains (domain_A, domain_B, domains_co, domains_contra, distance_out):

	infile = open(domain_A, 'r').readlines()
	outfile = open(distance_out, 'w')

	output_A = dict()
        output_B = dict()
	if (len(domains_contra) > len(domains_co)):
		right = domains_contra
	else:
		right = domains_co

	for rec in right:
		outfile.write(rec[0])
		output_A[rec[0]] = []
		for y in infile:
			x = str(y)
			domain = x.split(':')[0]
    			prot = x.split(':')[1].split(' ')

			dist = rec[1].get(str(domain.split(' ')[0].split('\n')[0]))
			output_A[rec[0]] += [domain + ' - ' + str(len(prot)) + ' matches, distance = ' + str(dist) + '\n']

	infile = open(domain_B, 'r').readlines()
        outfile = open(distance_out, 'a')

	for rec in right:
                outfile.write(rec[0])
		output_B[rec[0]] = []
                for y in infile:
                        x = str(y)
                        domain = x.split(':')[0]
                        prot = x.split(':')[1].split(' ')

                        dist = rec[1].get(str(domain.split(' ')[0].split('\n')[0]))
			output_B[rec[0]] += [domain + ' - ' + str(len(prot)) + ' matches, distance = ' + str(dist) + '\n']
	return output_A, output_B


##	Creates output files in csv format, files will contain find_interacting_domains output
##	pdb_A, pdb_B - PDB ID's of the Primes
##	domains_A, domains_B - lists of prime_A and prime_B domains, often given from find_interacting_domains output
##	dist_A, dist_B - prime_domains output, just for testing with zlab proteins
##	prime_A, prime_B - files with prime's formatted (with header) sequences

def output_binder1 (interpro, pdb_A, pdb_B, fout, domains_A, domains_B, dist_A, dist_B, prime_A, prime_B, domain_file_A, domain_file_B):
	
	tmp_prime_domains_A = get_domains_list (prime_A, interpro)
	tmp_prime_domains_B = get_domains_list (prime_B, interpro)
	
	out_file_A = open(fout + '_receptor.csv', 'w')
	out_file_A.write('domain,matches,all_matches,start,end\n')

	print prime_A
	print tmp_prime_domains_A
	
	prime_domains_A = list()
	prime_domains_B = list()

	for i in tmp_prime_domains_A:
		prime_domains_A.append(i[0])
	for i in tmp_prime_domains_B:
                prime_domains_B.append(i[0])

	print prime_domains_A	

	file_A = open(domain_file_A, 'r').readlines()
	uni_set = set()
        for rec in file_A:
                x = rec.split(':')[1].split(' ')
                for j in x:
                	uni_set.update([j])
	all_matches = len(uni_set)
	general_out_A = list()
	interpro_out = open('../tmp/interpro_' + prime_A, 'r').readlines()
	print ('General frequency information')
        x = open(domain_file_A, 'r').readlines()
        for rec in x:
                matches = len(rec.split(':')[1].split(' ')) - 1
		domain = rec.split(' ')[0]
		print domain
		if domain in prime_domains_A:
			pos = []
			for i in interpro_out:
				if len(i.split('\t')) == 13:
					if (i.split('\t')[11] == domain):
						pos = [i.split('\t')[6], i.split('\t')[7]]
			general_out_A.append([domain, matches, all_matches, pos[0], pos[1]])	
	print general_out_A
	general_out_A = np.array(general_out_A)
	general_out_A = general_out_A[np.argsort(general_out_A[:,1])[::-1]]
	for i in general_out_A:
		s = ''
		flag = 0
		for j in i:
			s += str(j)
			flag += 1
			if (flag != len(i)):
				s += ','
		out_file_A.write(s+'\n')
	print ('Addition complex structure information')
        for i in dist_A:
                out = []
                for j in dist_A[i]:
                        domain = j.split(' ')[0]
                        match = j.split('-')[1].split('matches')[0].split(']')[0].split('[')[0]
                        dist = j.split('[')#[1].split('\n')[0].split(',')
                        if (len(dist) == 1):
                                dist = None
                                for k in domains_A:
                                        if (k[0] == domain):
                                                dist = [k[1], k[2], None]
                        else:
                                dist = j.split('[')[1].split('\n')[0].split(',')
                                dist = [dist[0].split(' ')[0], dist[0].split(' ')[1], dist[0].split(' ')[2], dist[0].split(' ')[3]]#dist[1].split(' ')[1]]

                        if (dist != None):
                                if (len(dist) == 3):
                                        out.append([domain, int(match.split(']')[0]), None, None, dist[0], dist[1]])
                                else:
                                        out.append([domain, int(match.split(']')[0]), dist[0], dist[1], dist[2], dist[3]])
                out = np.array(out)
                out = out[np.argsort(out[:,1])[::-1]]
                uni_set = set()
                file_A = open(domain_file_A, 'r').readlines()
                for rec in file_A:
                        x = rec.split(':')[1].split(' ')
                        for j in x:
                                uni_set.update([j])
                seq_len = len(open(prime_A, 'r').readlines()[1])
                for j in out:
                        print (i.split('\n')[0] + ',' + str(j[0]) + ',' + str(len(uni_set)) + ',' + str(j[1]) + ',' + str(j[2]) + ',' + str(j[3]) + ',' + str(j[4]) + ',' + str(j[5]))



	out_file_B = open(fout + '_ligand.csv', 'w')
        out_file_B.write('domain,matches,all_matches,start,end\n')

        print prime_B
        print tmp_prime_domains_B

        file_B = open(domain_file_B, 'r').readlines()
        uni_set = set()
        for rec in file_B:
                x = rec.split(':')[1].split(' ')
                for j in x:
                        uni_set.update([j])
        all_matches = len(uni_set)
        general_out_B = list()
        interpro_out = open('../tmp/interpro_' + prime_B, 'r').readlines()
        print ('General frequency information')
        x = open(domain_file_B, 'r').readlines()
        for rec in x:
                matches = len(rec.split(':')[1].split(' ')) - 1
                domain = rec.split(' ')[0]
                print domain
                if domain in prime_domains_B:
                        pos = []
                        for i in interpro_out:
                                if len(i.split('\t')) == 13:
                                        if (i.split('\t')[11] == domain):
                                                pos = [i.split('\t')[6], i.split('\t')[7]]
                        general_out_B.append([domain, matches, all_matches, pos[0], pos[1]])
        print general_out_B
        general_out_B = np.array(general_out_B)
        general_out_B = general_out_B[np.argsort(general_out_B[:,1])[::-1]]
	for i in general_out_B:
                s = ''
                flag = 0
                for j in i:
                        s += str(j) 
                        flag += 1
                        if (flag != len(i)):
                                s += ','
                out_file_B.write(s+'\n')


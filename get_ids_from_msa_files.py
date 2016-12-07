#from domain_work import *
#from get_ids_from_msa_files import *
#from id_converters import *
#from interaction_search import *
#from script_running import *
#from binder import *
import numpy as np
import pandas as pd
from Bio import SeqIO
from multiprocessing import Pool
from Bio import SeqIO
import ConfigParser
from Bio.PDB import *
import os

##Simple alignment for full sequence and part of this sequence

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

##Gets list of PDB ID's + Prime ID from the header from pdb_msa file

def get_pdb_ids (infile, overlap):

	## Get original ID's from file
        flag = 1
        pdb_ids = set()
        uniprot_ids = set()
	prime = ''
        ## Get pdb ids from pdb_MSA file
        handle = open('pdb_' + infile, 'r')

        for record in SeqIO.parse(handle, "fasta"):
                buf = record.id
                if (seq_test(record.seq, overlap) | flag):
                        if (flag):
                                pdb_ids.update ([buf.split(':')[0].lower() + '_' + buf.split(':')[1].split('|')[0]])
                                flag = 0
				prime = buf.split(':')[0].lower() + '_' + buf.split(':')[1].split('|')[0]
                        else:
                                pdb_ids.update([buf.split('_')[0].lower() + '_' + buf.split('_')[1]])
	return prime, pdb_ids 


##Gets Uniprot Gene ID's from both pdb_ and uni_ msa files

def get_biogrid_chain_ids (infile, overlap):

	## Get original ID's from file
        flag = 1
        pdb_ids = set()
        uniprot_ids = set()

	## Get pdb ids from pdb_ MSA file
        handle = open('pdb_' + infile, 'r')

        for record in SeqIO.parse(handle, "fasta"):
                buf = record.id
		if (seq_test(record.seq, overlap) | flag):
                	if (flag):
				pdb_ids.update ([buf.split(':')[0].lower() + '_' + buf.split(':')[1].split('|')[0]])
                        	flag = 0
                	else:
				pdb_ids.update([buf.split('_')[0].lower() + '_' + buf.split('_')[1]])
	
	flag = 1
	## Get Uniprot ID's from uni_ msa file
	handle = open('uni_' + infile, 'r')

        for record in SeqIO.parse(handle, "fasta"):
                buf = record.id
                if (seq_test(record.seq, overlap) | flag):
                        if (flag):
                                pdb_ids.update ([buf.split(':')[0].lower() + '_' + buf.split(':')[1].split('|')[0]])
                                flag = 0
                        else:
                                uniprot_ids.update([buf.split('|')[2]])
	buf_ids = set()
	for x in pdb_ids:
		os.system ('python get_chain_uni_db.py ' + str(x))
                tmp = open('../tmp/chain_uni_out.txt', 'r')
                for i in tmp:
                        buf_ids.update([i.split('\n')[0]])

	for x in buf_ids:
		os.system ('python get_uni_id.py ' + str(x))
                tmp = open('../tmp/uni_id_out.txt', 'r')
                for i in tmp:
                        uniprot_ids.update([i.split('\n')[0]])
	
	result = (list(uniprot_ids))
        return result


##Gets String ID's from both pdb_ and uni_ msa files

def get_string_ids_sql (infile, overlap):

	## Get uniprot notation

	uniprot_ids = get_biogrid_chain_ids (infile, overlap)
	string_ids = set()
	
	for x in uniprot_ids:
		os.system ('python get_full_uni_db.py ' + str(x))
		tmp = open('../tmp/full_out.txt', 'r')
		for i in tmp:
			string_ids.update([i.split('\n')[0]])
	
	return string_ids



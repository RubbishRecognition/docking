#from domain_work import *
#from get_ids_from_msa_files import *
#from id_converters import *
#from interaction_search import *
#from script_running import *
#from binder import *
import os

def run_biogrid_finder (data):
	
	## Find biogrid interactions using finder
		
        query_string = './finder.o ' + 'B ' + data[1] + ' ' + data[0].split('_')[0]
        #print query_string
	os.system (query_string)


def run_string_finder (data):

	## Find string interactions using finder	

        query_string = './finder.o ' + 'S ' + data[1] + ' ' + data[0]
        #print query_string
	os.system (query_string)


def run_string_finder_sql (gene):
	
	## Find string interactions using SQL finder

	query_string = 'python find_sqlite_db_paired.py ' + gene.split('\n')[0]
	os.system (query_string)


def run_uni_id_finder_sql (gene):

	## Find uniprot ID

	query_string = 'python get_uni_id.py ' + gene
	os.system (query_string)


def run_biogrid_id_sql (gene):

	## Find Biogrid gene ID from Uniprot chain ID

	query_string = 'python get_chain_uni_db.py ' + gene
	os.system (query_string)


def run_hhblits (uniprot20, pdb100, cpu, score, sequence, ident_file, align_file):
	
	## Find similar proteins

	query_string = 'hhblits -i ' + str(sequence)  + ' -d ' + str(uniprot20) + ' -oa3m uni_' + str(ident_file) + ' -cpu ' + str(cpu) + ' -qid ' + str(score) + ' -id 100 ' + '-v 0' + ' -o ' + 'uni_' + str(align_file)  
	os.system (query_string)	
	query_string = 'hhblits -i ' + str(sequence) + ' -d ' + str(pdb100) + ' -oa3m pdb_' + str(ident_file) + ' -cpu ' + str(cpu) + ' -qid ' + str(score) + ' -id 100 ' + '-v 0' + ' -o ' + 'pdb_' + str(align_file) 
	os.system (query_string)
	print query_string


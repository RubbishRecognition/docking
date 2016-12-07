#from domain_work import *
#from get_ids_from_msa_files import *
#from id_converters import *
#from interaction_search import *
#from script_running import *
#from binder import *
import os

##Converts Uniprot ID like P14778 to list of it's domains [IPR000157, ...]
##Using InterPro Domain Database

def get_domain_id (uniprot_id):

        domain_id = set()
        query_string = 'python get_domain_db.py ' + uniprot_id
        os.system (query_string)
        tmp = open('../tmp/domain_out.txt', 'r')
        for i in tmp:
                domain_id.update([i.split('\n')[0]])

        return list (domain_id)


##Coverts Uniprot ID like Q00959 into it's gene ID Grin2a and vice versa
##Using Uniprot database

def get_gene_id (uniprot):

        os.system ('python get_uniprot_gene_id.py ' + str(uniprot))
        tmp = open('../tmp/uniprot_gene_out.txt', 'r')
	res = []
        for i in tmp:
		res.append(i.split('\n')[0])
	if (len(res) == 0):
                return uniprot
        else:
		if (len(res) == 1):
                	return res[0]
		else:
			return res


##Converts PDB chain ID like 1g0y_R into it's Uniprot ID P14778 and vice versa
##Using pdb_chain_uniprot database

def uniprot_pdb_chain_converter (gene):


	os.system ('python get_chain_uni_db.py ' + str(gene))
        tmp = open('../tmp/chain_uni_out.txt', 'r')
        res = []
        for i in tmp:
                res.append(i.split('\n')[0])
        if (len(res) == 0):
                return uniprot
        else:
                if (len(res) == 1):
                        return res[0]
                else:
                        return res


##Converts Uniprot ID like P14778 into it's uniprot.org Entry name IL1R1_HUMAN and vice versa
##Using pdbtosp database

def get_uniprot_id (gene):
	os.system ('python get_uni_id.py ' + str(gene))
        tmp = open('../tmp/uni_id_out.txt', 'r')
	res = []
        for i in tmp:
        	#return i.split('\n')[0]
		res.append(i.split('\n')[0])
	if (len(res) == 0):
		return gene
	else:
		return res[0]


##Converts list of Entry names like [IL1R1_HUMAN, ...] into list of its pdb chain ID's [1g0y_R, ...] 

def uniprot_to_pdb (uniprot_ids):

	pdb_ids = set()
	buf_ids = set()
	
        for x in uniprot_ids:
                os.system ('python get_uni_id.py ' + str(x))
                tmp = open('../tmp/uni_id_out.txt', 'r')
                for i in tmp:
                        buf_ids.update([i.split('\n')[0]])

        for x in buf_ids:
                os.system ('python get_chain_uni_db.py ' + str(x))
                tmp = open('../tmp/chain_uni_out.txt', 'r')
                for i in tmp:
                        pdb_ids.update([i.split('\n')[0]])
        
	result = (list(pdb_ids))	
	return result


##Gets list of domains which are in sequence in file 'gene'

def get_domains_list (gene, interpro):

	## Get domains ids with positions

	domain_id = list()
	os.system (interpro + '/interproscan.sh -i ' + gene + ' -iprlookup -o ../tmp/interpro_' + gene + ' -f TSV')
        tmp = open ('../tmp/interpro_' + gene, 'r').readlines()
	for i in tmp:
		if len(i.split('\t')) == 13:
			domain_id.append([i.split('\t')[11], i.split('\t')[6], i.split('\t')[7]])
        return domain_id



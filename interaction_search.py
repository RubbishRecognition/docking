#from domain_work import *
#from get_ids_from_msa_files import *
from id_converters import *
#from interaction_search import *
from script_running import *
#from binder import *
#from id_converters import get_uniprot_id
from multiprocessing import Pool
from functools import partial

##Gets list of pdb interactions in PDB ID format from 2 lists of proteins in Uniprot Gene format

def find_biogrid_interactions (ids_A, ids_B, cpu, biogrid):
	
	interactions = set()
        pool = Pool(cpu)

	grid_ids_1 = [get_uniprot_id(j) for j in ids_A]
	grid_ids_A = [[get_gene_id (i), biogrid] for i in grid_ids_1]
	grid_ids_1 = [get_uniprot_id(j) for j in ids_B]
        grid_ids_B = [get_gene_id (i) for i in grid_ids_1]
	pool.map(run_biogrid_finder, grid_ids_A)
	grid_ids_1 = [get_uniprot_id(j) for j in ids_A]
        grid_ids_A = [get_gene_id (i) for i in grid_ids_1]
        for gene in grid_ids_A:
                output = open ('../tmp/' + gene.split('_')[0], 'r')
                for entry in output:
                        for gene_B in grid_ids_B:
				if (type(gene_B) == list):
					gene_B = gene_B[0]
                               	if (entry.translate(None, '\n') == gene_B.split('_')[0]):
					res_A = get_gene_id(gene)
					res_B = get_gene_id(gene_B)
					if (type(res_A) == str):
						res_A = [res_A]
					if (type(res_B) == str):
						res_B = [res_B]
					for i in range(len(res_A)):
						res_A[i] = get_uniprot_id(res_A[i])
					for i in range(len(res_B)):
                               	                res_B[i] = get_uniprot_id(res_B[i])
					for k in res_A:
						for l in res_B:
							interactions.update([k + ' ' + l + '\n'])

	return list(interactions)


##Gets list of string interactions in String ID format from 2 lists of proteins in String ID format

def find_string_interactions (ids_A, ids_B, cpu, string):

        interactions = set()
        pool = Pool(cpu)

	#formatted_ids_A = [[i, string] for i in list(ids_A)]

        pool.map(run_string_finder_sql, list(ids_A))

        for gene in ids_A:
                output = open ('../tmp/' + gene.split('\n')[0], 'r')
                for entry in output:
                        for gene_B in ids_B:
                                if (entry.translate(None, '\n') == gene_B):
                                        interactions.update([str(gene + ' ' + gene_B + '\n')])
        return list(interactions)


##Converts list of interactions in String format into list of interactions in Uniprot format 

def find_uniprot_id_sql (interactions, in_format):

        ## String to PDB

        uniprot_interactions = set()

        for x in interactions:
        	interaction = ''
                pat1 = set()
                pat2 = set()

                prot1 = x.split(' ')[0]
                prot2 = x.split(' ')[1]

                prot2 = prot2.split('\n')[0]

                        ## Map first protein                    

                os.system ('python get_full_uni_db.py ' + str(prot1))
                tmp = open('../tmp/full_out.txt', 'r')
                for i in tmp:
                	pat1.update([i.split('\n')[0]])

                        ## Map second protein

                os.system ('python get_full_uni_db.py ' + str(prot2))
                tmp = open('../tmp/full_out.txt', 'r')
                for i in tmp:
                	pat2.update([i.split('\n')[0]])

                for p1 in pat1:
                	for p2 in pat2:
                        	uniprot_interactions.update([str(p1) + ' ' + str(p2) + '\n'])
	return list(uniprot_interactions)


##Converts list of interactions in String format (in_format = 'S') or list of interactions in Uniprot format (in_format = 'U') into list of interactions of PDB format

def find_pdb_id_sql (interactions, in_format):

	if (in_format == 'S'):

        ## String to PDB

                # Read string databases

		uniprot_interactions = set()

                for x in interactions:

                        interaction = ''
			pat1 = set()
                        pat2 = set()

                        prot1 = x.split(' ')[0]
                        prot2 = x.split(' ')[1]

                        prot2 = prot2.split('\n')[0]
	
			## Map first protein			

			os.system ('python get_full_uni_db.py ' + str(prot1))
                	tmp = open('../tmp/full_out.txt', 'r')
                	for i in tmp:
                        	pat1.update([i.split('\n')[0]])

			## Map second protein

			os.system ('python get_full_uni_db.py ' + str(prot2))        
			tmp = open('../tmp/full_out.txt', 'r')
                        for i in tmp:
                                pat2.update([i.split('\n')[0]])
			
			for p1 in pat1:
                                for p2 in pat2:
                                        uniprot_interactions.update([str(p1) + ' ' + str(p2) + '\n'])
		pdb_interactions = set()

                for x in uniprot_interactions:

                        interaction = ''

                        prot1 = x.split(' ')[0]
                        prot2 = x.split(' ')[1]
                        prot2 = prot2.split('\n')[0]

                        pat1 = list(set(uniprot_to_pdb ([prot1])))
                        pat2 = list(set(uniprot_to_pdb ([prot2])))

                        for p1 in pat1:
                                for p2 in pat2:
                                        pdb_interactions.update([str(p1) + ' ' + str(p2) + '\n'])

                return list(pdb_interactions)
	else:
                # Uniprot to pdb                

                pdb_interactions = set()

                for x in interactions:

                        interaction = ''

                        prot1 = x.split(' ')[0]
                        prot2 = x.split(' ')[1]
                        prot2 = prot2.split('\n')[0]

                        pat1 = list(set(uniprot_to_pdb ([prot1])))
                        pat2 = list(set(uniprot_to_pdb ([prot2])))

                        for p1 in pat1:
                                for p2 in pat2:
                                        pdb_interactions.update([str(p1) + ' ' + str(p2) + '\n'])

                return list(pdb_interactions)





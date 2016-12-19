import os
import argparse
import numpy as np
import pandas as pd
from Bio import SeqIO
from multiprocessing import Pool
from Bio import SeqIO
import ConfigParser
from Bio.PDB import *

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

score = int(ConfigSectionMap("hhblits options")['score'])
cpu = int(ConfigSectionMap("hhblits options")['cpu'])
overlap = int (ConfigSectionMap("hhblits options")['overlap'])
#print score + 1
#k = input()

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


def get_domain_id (uniprot_id):

        domain_id = set()

        query_string = 'python get_domain_db.py ' + uniprot_id
        os.system (query_string)
        tmp = open('../tmp/domain_out.txt', 'r')
        for i in tmp:
                domain_id.update([i.split('\n')[0]])

        return list (domain_id)


def run_biogrid_finder (gene):

        ## Find biogrid interactions using finder

        query_string = './finder.o ' + 'B ' + biogrid + ' ' + gene.split('_')[0]
        #print query_string
        os.system (query_string)


def uniprot_pdb_chain_converter (gene):

        ## ID converter unsing pdb_chain_uniprot db

        os.system ('python get_chain_uni_db.py ' + str(gene))
        tmp = open('../tmp/chain_uni_out.txt', 'r')
        #f = tmp.readlines()
        res = []
        #if (len(f) > 0):
        for i in tmp:
                res.append(i.split('\n')[0])
        if (len(res) == 0):
                return uniprot
        else:
                if (len(res) == 1):
                        return res[0]
                else:
                        return res


def get_gene_id (uniprot):

        #for i in range(len(biogrid_ids_A)):
        os.system ('python get_uniprot_gene_id.py ' + str(uniprot))
        tmp = open('../tmp/uniprot_gene_out.txt', 'r')
        #f = tmp.readlines()
        res = []
        #if (len(f) > 0):
        for i in tmp:
                res.append(i.split('\n')[0])
        if (len(res) == 0):
                return uniprot
        else:
                if (len(res) == 1):
                        return res[0]
                else:
                        return res


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


def run_string_finder (gene):

        ## Find string interactions using finder

        query_string = './finder.o ' + 'S ' + string + ' ' + gene
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


def run_hhblits (sequence, ident_file, align_file):

        ## Find similar proteins

        #query_string = 'hhblits -i ' + str(sequence)  + ' -d ' + str(uniprot20) + ' -d ' + str(pdb100) + ' -oa3m ' + str(ident_file) + ' -cpu ' + str(cpu) + ' -qid ' + str(score) + ' -id 100 ' + '-v 0' + ' -o ' + str(align_file)
        #CHANGE 95 TO SCORE!!!!!!!
        query_string = 'hhblits -i ' + str(sequence)  + ' -d ' + str(uniprot20) + ' -oa3m uni_' + str(ident_file) + ' -cpu ' + str(cpu) + ' -qid ' + str(score) + ' -id 100 ' + '-v 0' + ' -o ' + 'uni_' + str(align_file)
        os.system (query_string)
        query_string = 'hhblits -i ' + str(sequence) + ' -d ' + str(pdb100) + ' -oa3m pdb_' + str(ident_file) + ' -cpu ' + str(cpu) + ' -qid ' + str(score) + ' -id 100 ' + '-v 0' + ' -o ' + 'pdb_' + str(align_file)
        os.system (query_string)
        ## Create alignment
        print query_string
        #query_string = 'hhblits -i ' + sequence + ' -d ' + uniprot20 + ' -d ' + pdb100 + ' -o ' + align_file + ' -cpu 5 -score -qid 70 -id 0 -v 0'
        #os.system (query_string)


def pdb_to_uniprot (pdb_ids):

        ## Get uniprot ID's from pdb ID's

        #Get all uniprot ids
        full_pdb = []
        full_uniprot = []
        uniprot_ids = []
        n = 0
        x = open(pdbtosp, 'r')
        for i in x:
                n += 1
                if (n > 24) and (i != '\n') and (i != '-----------------------------------------------------------------------\n'):
                        pattern = i.split(' ')
                        if (pattern[0] != ''):
                                full_pdb.append (pattern[0])
                                k = 5
                                while (len (pattern[k]) < 6):
                                        k += 1
                                if (len(pattern[k].split('_')) == 2):
                                        full_uniprot.append (pattern[k].split('_')[0])
        #Map ID's
        for i in range(len(pdb_ids)):
                for index, value in enumerate(full_pdb):
                        if (value == pdb_ids[i]):
                                uniprot_ids.append(full_uniprot[index])

        return uniprot_ids


def uniprot_to_pdb (uniprot_ids):

        ## Get uniprot ID's from pdb ID's
        pdb_ids = set()
        '''
        #Get all uniprot ids
        full_pdb = []
        full_uniprot = []
        pdb_ids = []
        n = 0
        x = open(pdbtosp, 'r')
        for i in x:
                n += 1
                if (n > 24) and (i != '\n') and (i != '-----------------------------------------------------------------------\n'):
                        pattern = i.split(' ')
                        if (pattern[0] != ''):
                                full_pdb.append (pattern[0])
                                k = 5
                                while (len (pattern[k]) < 6):
                                        k += 1
                                if (len(pattern[k].split('_')) == 2):
                                        full_uniprot.append (pattern[k].split('_')[0])
        #Map ID's
        for i in range(len(uniprot_ids)):
                for index, value in enumerate(full_uniprot):
                        if (value == uniprot_ids[i]):
                                pdb_ids.append(full_pdb[index])

        '''

        buf_ids = set()

        #print 'first'
        #print uniprot_ids

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

        #uniprot_ids = list(set(list(uniprot_ids) + list(pdb_to_uniprot (pdb_ids)))) #restyling
        result = (list(pdb_ids))

        return result


def get_pdb_ids (infile):

        ## Get original ID's from file
        flag = 1
        pdb_ids = set()
        uniprot_ids = set()
        prime = ''
        ## Get pdb ids from pdb_ MSA file
        handle = open('pdb_' + infile, 'r')

        for record in SeqIO.parse(handle, "fasta"):
                buf = record.id
                #print(len(buf))
                if (seq_test(record.seq, overlap) | flag):
                        if (flag):
                                #pdb_ids.append (buf.split(':')[0]) #restyling
                                pdb_ids.update ([buf.split(':')[0].lower() + '_' + buf.split(':')[1].split('|')[0]])
                                flag = 0
                                prime = buf.split(':')[0].lower() + '_' + buf.split(':')[1].split('|')[0]
                        else:
                                #if (len(buf) <= 8):
                                #print pdb_ids
                                #print buf
                                #pdb_ids.append(buf.split('_')[0].upper()) #restyling
                                pdb_ids.update([buf.split('_')[0].lower() + '_' + buf.split('_')[1]])
                                #else:
                                        #print buf
                                        #uniprot_ids.update([buf.split('|')[2]])
        return prime, pdb_ids



def get_biogrid_chain_ids (infile):

        ## The same as get_biogrid_ids but with another db for mapping

        ## Get original ID's from file
        flag = 1
        pdb_ids = set()
        uniprot_ids = set()

        ## Get pdb ids from pdb_ MSA file
        handle = open('pdb_' + infile, 'r')

        for record in SeqIO.parse(handle, "fasta"):
                buf = record.id
                #print(len(buf))
                if (seq_test(record.seq, overlap) | flag):
                        if (flag):
                                #pdb_ids.append (buf.split(':')[0]) #restyling
                                pdb_ids.update ([buf.split(':')[0].lower() + '_' + buf.split(':')[1].split('|')[0]])
                                flag = 0
                        else:
                                #if (len(buf) <= 8):
                                        #pdb_ids.append(buf.split('_')[0].upper()) #restyling
                                pdb_ids.update([buf.split('_')[0].lower() + '_' + buf.split('_')[1]])
                                #else:
                                        #print buf
                                        #uniprot_ids.update([buf.split('|')[2]])

        flag = 1
        ## Get uniprot ids from uni_ MSA file
        handle = open('uni_' + infile, 'r')

        for record in SeqIO.parse(handle, "fasta"):
                buf = record.id
                #print(len(buf))
                if (seq_test(record.seq, overlap) | flag):
                        if (flag):
                                #pdb_ids.append (buf.split(':')[0]) #restyling
                                pdb_ids.update ([buf.split(':')[0].lower() + '_' + buf.split(':')[1].split('|')[0]])
                                flag = 0
                        else:
                                #if (len(buf) <= 8):
                                        #pdb_ids.append(buf.split('_')[0].upper()) #restyling
                                #       pdb_ids.update([buf.split('_')[0].lower() + '_' + buf.split('_')[1]])
                                #else:
                                        #print buf
                                uniprot_ids.update([buf.split('|')[2]])

        #print pdb_ids
        #print uniprot_ids
        #print '##'
        #k = input()
        ## TO DO: 1) pdb to uni_id, 2) uni_id to biogrid_gene (to check: with _HUMAN or without)
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

        #uniprot_ids = list(set(list(uniprot_ids) + list(pdb_to_uniprot (pdb_ids)))) #restyling
        result = (list(uniprot_ids))
        #print result
        #k = input()
        return result



def get_biogrid_ids (infile):

        ## Get original ID's from file
        flag = 1
        pdb_ids = []
        uniprot_ids = []
        handle = open(infile, 'r')

        for record in SeqIO.parse(handle, "fasta"):
                buf = record.id
                #print(len(buf))

                if (flag):
                        pdb_ids.append (buf.split(':')[0])
                        flag = 0
                else:
                        if (len(buf) == 6):
                                pdb_ids.append(buf.split('_')[0].upper())
                        else:
                                uniprot_ids.append(buf.split('|')[2])

        ## ID mapping
        #print 'kek'
        #print infile
        #print uniprot_ids
        #print pdb_ids
        #print 'kek'
        #print pdb_to_uniprot (pdb_ids)
        #print "kek"
        #k = input()
        uniprot_ids = list(set(list(uniprot_ids) + list(pdb_to_uniprot (pdb_ids))))

        return uniprot_ids



def get_string_ids_sql (infile):

        ## Get uniprot notation

        #uniprot_ids = get_biogrid_ids (infile) // previous variant
        uniprot_ids = get_biogrid_chain_ids (infile)

        string_ids = set()

        for x in uniprot_ids:
                os.system ('python get_full_uni_db.py ' + str(x))
                tmp = open('../tmp/full_out.txt', 'r')
                for i in tmp:
                        string_ids.update([i.split('\n')[0]])
        '''
        for x in uniprot_ids:
                os.system ('python get_all_uni_db.py ' + str(x))
                tmp = open('../tmp/all_out.txt', 'r')
                for i in tmp:
                        string_ids.update([i.split('\n')[0]])
        ''' #restyling
        return string_ids



def get_string_ids (infile):

        ## Get uniprot notation

        uniprot_ids = get_biogrid_ids (infile)

        ## Read string databases

        d = dict()

        #data = pd.read_csv('../id_mapping/full_uniprot_2_string.04_2015.tsv', sep = '\t')
        #print full_uniprot
        data = pd.read_csv (full_uniprot, sep = '\t')
        data = np.array(data)
        for i in range(data.shape[0]):
                #print data[i][1]
                buf = d.setdefault(data[i][1].split('|')[1].split('_')[0])
                if (buf == None):
                        d[data[i][1].split('|')[1].split('_')[0]] = [data[i][2]]
                else:
                        buf.append(data[i][2])
                        d[data[i][1].split('|')[1].split('_')[0]] = buf

        #data = pd.DataFrame.from_csv('../id_mapping/all_go_knowledge_full.tsv', sep = '\t')
        data = pd.DataFrame.from_csv (all_knowledge, sep = '\t')
        data = np.array(data)

        for i in range(data.shape[0]):
                buf = d.setdefault(data[i][1])
                if (buf == None):
                        d[data[i][1]] = [data[i][0]]
                else:
                        buf.append(data[i][0])
                        d[data[i][1]] = buf

        ## ID mapping

        string_ids = set()

        for x in uniprot_ids:
                buf = d.setdefault(x)
                if (buf != None):
                        string_ids.update(buf)
                        #string_ids.update([buf])

                buf = d.setdefault(x.split('_')[0])
                if (buf != None):
                        string_ids.update(buf)
                        #string_ids.update([buf])

        return string_ids


def find_biogrid_interactions (ids_A, ids_B):

        ## Find biogrid interactions using run_biogrid_finder in 16 threads

        interactions = set()
        pool = Pool(cpu)

        #pool.map(run_biogrid_finder, list(ids_A))
        #grid_ids_A = [get_gene_id(i) for i in ids_A]
        grid_ids_1 = [get_uniprot_id(j) for j in ids_A]
        grid_ids_A = [get_gene_id (i) for i in grid_ids_1]
        grid_ids_1 = [get_uniprot_id(j) for j in ids_B]
        grid_ids_B = [get_gene_id (i) for i in grid_ids_1]
        #print grid_ids_A
        #print grid_ids_B
        #print 'kek'
        #k = input()
        print grid_ids_A
        pool.map(run_biogrid_finder, grid_ids_A)#[get_uniprot_id (i) for i in [get_gene_id(j) for j in ids_A])
        for gene in grid_ids_A:
                output = open ('../tmp/' + gene.split('_')[0], 'r')
                for entry in output:
                        for gene_B in grid_ids_B:
                                #print gene_B
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


def find_string_interactions (ids_A, ids_B):

        ## Find string interactions using run_string_finder in 'cpu' threads

        interactions = set()
        pool = Pool(cpu)

        pool.map(run_string_finder_sql, list(ids_A))

        for gene in ids_A:
                output = open ('../tmp/' + gene.split('\n')[0], 'r')
                for entry in output:
                        for gene_B in ids_B:
                                if (entry.translate(None, '\n') == gene_B):
                                        interactions.update([str(gene + ' ' + gene_B + '\n')])


        return list(interactions)



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
                        '''
                        os.system ('python get_all_uni_db.py ' + str(prot1))
                        tmp = open('../tmp/all_out.txt', 'r')
                               for i in tmp:
                                pat1.update([i.split('\n')[0]])
                        '''
                        ## Map second protein

                        os.system ('python get_full_uni_db.py ' + str(prot2))
                        tmp = open('../tmp/full_out.txt', 'r')
                        for i in tmp:
                                pat2.update([i.split('\n')[0]])
                        '''
                        os.system ('python get_all_uni_db.py ' + str(prot2))
                        tmp = open('../tmp/all_out.txt', 'r')
                        for i in tmp:
                                pat2.update([i.split('\n')[0]])
                        '''
                        for p1 in pat1:
                                for p2 in pat2:
                                        uniprot_interactions.update([str(p1) + ' ' + str(p2) + '\n'])

                        #print x
                        #print pat1
                        #print pat2
                        #k = input()


                pdb_interactions = set()

                for x in uniprot_interactions:

                        interaction = ''

                        prot1 = x.split(' ')[0]
                        prot2 = x.split(' ')[1]
                        prot2 = prot2.split('\n')[0]

                        pat1 = list(set(uniprot_to_pdb ([prot1])))
                        pat2 = list(set(uniprot_to_pdb ([prot2])))
                        #interaction = str(uniprot_to_pdb (prot1)) + ' ' + str(uniprot_to_pdb (prot2))

                        #print interaction

                        for p1 in pat1:
                                for p2 in pat2:
                                        pdb_interactions.update([str(p1) + ' ' + str(p2) + '\n'])
                        #pdb_interactions.update([interaction])

                #print pdb_interactions

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
                        #interaction = str(uniprot_to_pdb (prot1)) + ' ' + str(uniprot_to_pdb (prot2))

                        #print interaction

                        for p1 in pat1:
                                for p2 in pat2:
                                        pdb_interactions.update([str(p1) + ' ' + str(p2) + '\n'])
                        #pdb_interactions.update([interaction])

                #print pdb_interactions

                return list(pdb_interactions)



def find_pdb_id (interactions, in_format):

        if (in_format == 'S'):

        ## String to PDB

                # Read string databases

                d = dict()

                #data = pd.read_csv('../id_mapping/full_uniprot_2_string.04_2015.tsv', sep = '\t')
                data = pd.read_csv (full_uniprot, sep = '\t')
                data = np.array(data)
                for i in range(data.shape[0]):
                        #buf = d.setdefault(data[i][1].split('|')[1])
                        #if (buf == None):
                        #        d[data[i][1].split('|')[1]] = [data[i][2]]
                        #else:
                        #        buf.append(data[i][2])
                        #        d[data[i][1].split('|')[1]] = buf


                        buf = d.setdefault(data[i][2])
                        if (buf == None):
                                d[data[i][2]] = [data[i][1].split('|')[1].split('_')[0]]
                        else:
                                buf.append(data[i][1].split('|')[1].split('_')[0])
                                d[data[i][2]] = buf

                #data = pd.DataFrame.from_csv('../id_mapping/all_go_knowledge_full.tsv', sep = '\t')
                data = pd.DataFrame.from_csv (all_knowledge, sep = '\t')
                data = np.array(data)
                for i in range(data.shape[0]):
                        #print data[i][0]
                        buf = d.setdefault(data[i][0])
                        if (buf == None):
                                d[data[i][0]] = [data[i][1]]
                        else:
                                buf.append(data[i][1])
                                d[data[i][0]] = buf

                # String to uniprot

                #print "BASE GENERATED"

                uniprot_interactions = set()

                for x in interactions:

                        interaction = ''

                        pat1 = set()
                        pat2 = set()

                        prot1 = x.split(' ')[0]
                        prot2 = x.split(' ')[1]

                        prot2 = prot2.split('\n')[0]

                        buf = d.setdefault(prot1)
                        if (buf != None):
                                interaction = list(set(buf))
                                pat1.update(buf)
                        else:
                                interaction = 'NONE'
                                pat1.update(['NONE'])

                        buf = d.setdefault(prot2)
                        if (buf != None):
                                pat2.update(buf)
                                interaction = str(list(set(interaction))) + ' ' + str(list(set(buf)))
                        else:
                                pat2.update(['NONE'])
                                interaction = str(list(set(interaction))) + ' NONE'

                        #print interaction

                        for p1 in pat1:
                                for p2 in pat2:
                                        uniprot_interactions.update([str(p1) + ' ' + str(p2) + '\n'])

                #print uniprot_interactions

                # Uniprot to pdb

                pdb_interactions = set()

                for x in uniprot_interactions:

                        interaction = ''

                        prot1 = x.split(' ')[0]
                        prot2 = x.split(' ')[1]
                        prot2 = prot2.split('\n')[0]

                        pat1 = list(set(uniprot_to_pdb ([prot1])))
                        pat2 = list(set(uniprot_to_pdb ([prot2])))
                        #interaction = str(uniprot_to_pdb (prot1)) + ' ' + str(uniprot_to_pdb (prot2))

                        #print interaction

                        for p1 in pat1:
                                for p2 in pat2:
                                        pdb_interactions.update([str(p1) + ' ' + str(p2) + '\n'])
                        #pdb_interactions.update([interaction])

                #print pdb_interactions

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
                        #interaction = str(uniprot_to_pdb (prot1)) + ' ' + str(uniprot_to_pdb (prot2))

                        #print interaction

                        for p1 in pat1:
                                for p2 in pat2:
                                        pdb_interactions.update([str(p1) + ' ' + str(p2) + '\n'])
                        #pdb_interactions.update([interaction])

                #print pdb_interactions

                return list(pdb_interactions)



def find_uniprot_id_sql (interactions, in_format):

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

                        '''
                        os.system ('python get_all_uni_db.py ' + str(prot1))
                        tmp = open('../tmp/all_out.txt', 'r')
                        for i in tmp:
                                pat1.update([i.split('\n')[0]])
                        '''
                        ## Map second protein

                        os.system ('python get_full_uni_db.py ' + str(prot2))
                        tmp = open('../tmp/full_out.txt', 'r')
                        for i in tmp:
                                pat2.update([i.split('\n')[0]])

                        '''
                        os.system ('python get_all_uni_db.py ' + str(prot2))
                        tmp = open('../tmp/all_out.txt', 'r')
                        for i in tmp:
                                pat2.update([i.split('\n')[0]])
                        '''
                        for p1 in pat1:
                                for p2 in pat2:
                                        uniprot_interactions.update([str(p1) + ' ' + str(p2) + '\n'])

                        #print x
                        #print pat1
                        #print pat2
                        #k = input()

                #print '##'
                #print interactions
                #print uniprot_interactions
                return list(uniprot_interactions)
                '''
                pdb_interactions = set()

                for x in uniprot_interactions:

                        interaction = ''

                        prot1 = x.split(' ')[0]
                        prot2 = x.split(' ')[1]
                        prot2 = prot2.split('\n')[0]

                        pat1 = list(set(uniprot_to_pdb ([prot1])))
                        pat2 = list(set(uniprot_to_pdb ([prot2])))
                        #interaction = str(uniprot_to_pdb (prot1)) + ' ' + str(uniprot_to_pdb (prot2))

                        #print interaction

                        for p1 in pat1:
                                for p2 in pat2:
                                        pdb_interactions.update([str(p1) + ' ' + str(p2) + '\n'])
                        #pdb_interactions.update([interaction])

                #print pdb_interactions

                return list(pdb_interactions)
                '''


def find_uniprot_id (interactions, in_format):

        if (in_format == 'S'):

        ## String to PDB

                # Read string databases

                d = dict()

                #data = pd.read_csv('../id_mapping/full_uniprot_2_string.04_2015.tsv', sep = '\t')
                data = pd.read_csv (full_uniprot, sep = '\t')
                data = np.array(data)
                for i in range(data.shape[0]):
                        #buf = d.setdefault(data[i][1].split('|')[1])
                        #if (buf == None):
                        #        d[data[i][1].split('|')[1]] = [data[i][2]]
                        #else:
                        #        buf.append(data[i][2])
                        #        d[data[i][1].split('|')[1]] = buf


                        buf = d.setdefault(data[i][2])
                        if (buf == None):
                                d[data[i][2]] = [data[i][1].split('|')[1].split('_')[0]]
                        else:
                                buf.append(data[i][1].split('|')[1].split('_')[0])
                                d[data[i][2]] = buf

                #data = pd.DataFrame.from_csv('../id_mapping/all_go_knowledge_full.tsv', sep = '\t')
                data = pd.DataFrame.from_csv (all_knowledge, sep = '\t')
                data = np.array(data)

                for i in range(data.shape[0]):
                        #print data[i][0]
                        buf = d.setdefault(data[i][0])
                        if (buf == None):
                                d[data[i][0]] = [data[i][1]]
                        else:
                                buf.append(data[i][1])
                                d[data[i][0]] = buf

                # String to uniprot

                #print "BASE GENERATED"

                uniprot_interactions = set()

                for x in interactions:

                        interaction = ''

                        pat1 = set()
                        pat2 = set()

                        prot1 = x.split(' ')[0]
                        prot2 = x.split(' ')[1]

                        prot2 = prot2.split('\n')[0]

                        buf = d.setdefault(prot1)
                        if (buf != None):
                                interaction = list(set(buf))
                                pat1.update(buf)
                        else:
                                interaction = 'NONE'
                                pat1.update(['NONE'])

                        buf = d.setdefault(prot2)
                        if (buf != None):
                                pat2.update(buf)
                                interaction = str(list(set(interaction))) + ' ' + str(list(set(buf)))
                        else:
                                pat2.update(['NONE'])
                                interaction = str(list(set(interaction))) + ' NONE'

                        #print interaction

                        for p1 in pat1:
                                for p2 in pat2:
                                        uniprot_interactions.update([str(p1) + ' ' + str(p2) + '\n'])

                return list(uniprot_interactions)


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

        #print full_uniprot

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
                for j in tmp:
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
                for j in tmp:
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


def get_domains_list (gene):

        ## Get domains ids with positions

        domain_id = list()

        query_string = 'python get_domain_list.py ' + gene
        os.system (query_string)
        tmp = open('../tmp/domain_list_out.txt', 'r')
        for i in tmp:
                domain_id.append(i.split('\n')[0].split(' '))

        return domain_id


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



def distance_calculator (domain_list, structure1, structure2, chain_A, chain_B):

        ## Distance calculator from domains of chain_A to chain_B
        distance_list = []
        #print domain_list
        #print chain_A
        #print chain_B
        #print '##'
        for model in structure1:
                for chain in model:
                        #print chain.id
                        if (chain.id == chain_A):
                                #print chain.id
                                for domain in domain_list:
                                        #print 'kwk'
                                        place = 0
                                        distance = float(0)
                                        num = 0
                                        int_num = 0
                                        for residue in chain:
                                                if ((place >= domain[1]) & (place <= domain[2])):
                                                        #min_dist = 99999999
                                                        for atom in residue:
                                                                min_dist = 99999999
                                                                for model2 in structure2:
                                                                        for chain2 in model2:
                                                                                #print 'idi nahui'
                                                                                if (chain2.id == chain_B):
                                                                                        #print 'heh'
                                                                                        for residue2 in chain2:
                                                                                                for atom2 in residue2:
                                                                                                        #print atom - atom2
                                                                                                        if (abs(atom - atom2) < min_dist):
                                                                                                                min_dist = abs(atom - atom2)
                                                                distance += min_dist
                                                                ## Means that atom exist in interation area
                                                                if (min_dist < 5):
                                                                        int_num += 1
                                                                num += 1
                                                place += 1
                                                #print residue
                                        #distance = distance / num
                                        #print place
                                        #print chain_A
                                        #print "######"
                                        #k = input()
                                        distance = str(int_num) +  ' ' + str(num)# + ' ' + str(1.0 * (domain[2] - domain[1]) / place) + ' ' + str(place) + ' ' + str(domain[1]) + ' ' + str(domain[2])
                                        distance_list.append(domain + [distance])
        #print 'kek'
        #print distance_list
        #k = input()
        return distance_list


def find_interacting_domains (align_A, align_B, complexes, prime_A, prime_B, pdb_ids_A, pdb_ids_B, outfile):

        ## Calculate distances and so on
        #print align_B
        #print align_A
        #print complexes
        #print prime_A
        #print pdb_ids_A
        uni_prime_A = uniprot_pdb_chain_converter (prime_A)
        uni_prime_B = uniprot_pdb_chain_converter (prime_B)
        result = []
        domains_A = get_domains_list (uni_prime_A)
        domains_B = get_domains_list (uni_prime_B)

        complex_list = open(complexes, 'r').readlines()
        complex_positions = list()
        #complex_list = open('/disk1/alekseev/docking/docking/scripts/d', 'r').readlines()
        #print '##'
        #print open('/disk1/alekseev/docking/docking/scripts/d', 'r').readlines()
        #print '##'
        #print complex_list
        #print pdb_ids_A
        #print pdb_ids_B
        #print complexes
        #print '##'                python binder.py ../fasta_data/1g0y.fasta.txt ../fasta_data/1ilr.fasta.txt pdb_1g0y_1ilr uniprot_1g0y_1ilr complex_1g0y_1ilr msa_1g0y hhr_1g0y msa_1ilr hhr_1ilr domains_1g0y domains_1ilr distance_1g0y_1ilr
        for comp in complex_list:
                #print comp
                chain_A = comp.split(':')[0]
                chain_B = comp.split('_')[0] + '_' + comp.split(':')[1].split('\n')[0]

                comp_domains_A = list()
                comp_domains_B = list()
                #print chain_A
                #print chain_B
                #print '###'
                if ((chain_A in pdb_ids_A) & (chain_B in pdb_ids_B) & (chain_A != chain_B)):
                        #print domains_A
                        #print domains_B
                        #print 'kekek'
                        ## Find pdb struct for chain_A
                        handle = open('pdb_' + align_A, 'r')
                        for record in SeqIO.parse(handle, "fasta"):
                                #print record.id
                                if (len(record.id) > 10):
                                        record.id = record.id.split(':')[0].lower() + '_' + record.id.split(':')[1].split('|')[0]
                                #print record.id
                                if (record.id == chain_A):
                                        start, end = find_align_pos (record.seq)
                                        #print [start, end]
                                        for domain in domains_A:
                                                #print domain[0]
                                                #print 'haha'
                                                if ((start <= int(domain[1])) & (end >= int(domain[2]))):
                                                        comp_domains_A.append([domain[0], int(domain[1]) - start, int(domain[2]) - start])

                        ## Find pdb struct for chain_B
                        handle = open('pdb_' + align_B, 'r')
                        for record in SeqIO.parse(handle, "fasta"):
                                #print 'kk'
                                if (len(record.id) > 10):
                                        record.id = record.id.split(':')[0].lower() + '_' + record.id.split(':')[1].split('|')[0]
                                #print record.id
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
                        #print domains_A
                        #print comp_domains_B
                        #print 'first'
                        distance_list_A = distance_calculator (comp_domains_B, structure1, structure2, chain_B.split('_')[1], chain_A.split('_')[1])
                        #print 'second'
                        structure1 = parser.get_structure('X', '../tmp/pdb' + comp.split('_')[0] + '.ent')
                        structure2 = parser.get_structure('X', '../tmp/pdb' + comp.split('_')[0] + '.ent')
                        distance_list_B = distance_calculator (comp_domains_A, structure1, structure2, chain_A.split('_')[1], chain_B.split('_')[1])
                        #print distance_list_A
                        #print distance_list_B
                        #print comp
                        #print prime_B
                        #k = input()
                        dist_dict = dict()
                        for i in distance_list_A:
                                buf = dist_dict.get(i[0])
                                if (buf == None):
                                        dist_dict[i[0]] = [i[3], i[2] - i[1]]
                                else:
                                        #if (float(i[3].split(' ')[0]) < float(buf.split(' ')[0])):
                                        dist_dict[i[0]] += [i[3], i[2] - i[1]]
                        for i in distance_list_B:
                                buf = dist_dict.get(i[0])
                                if (buf == None):
                                        dist_dict[i[0]] = [i[3], i[2] - i[1]]
                                else:
                                        #if (float(i[3].split(' ')[0]) < float(buf.split(' ')[0])):
                                        dist_dict[i[0]] += [i[3], i[2] - i[1]]
                        result.append([comp, dist_dict])

                        output = open(outfile, 'a')
                        output.write(comp + '\n')
                        for i in distance_list_A:
                                output.write(i[0] + ' ' + str(i[1]) + ' ' + str(i[2]) + ' ' + str(i[3]) + '\n')
                        for i in distance_list_B:
                                output.write(i[0] + ' ' + str(i[1]) + ' ' + str(i[2]) + ' ' + str(i[3]) + '\n')
                        #output.write(distance_list_A)
                        #output.write(distance_list_B)
                        output.write('\n\n')
        return domains_A, domains_B, result


def prime_domains (prime_A, prime_B, domain_A, domain_B, domains_co, domains_contra, distance_out):

        uni_prime_A = uniprot_pdb_chain_converter (prime_A)
        uni_prime_B = uniprot_pdb_chain_converter (prime_B)

        infile = open(domain_A, 'r').readlines()
        outfile = open(distance_out, 'w')

        if (len(domains_contra) > len(domains_co)):
                right = domains_contra
        else:
                right = domains_co

        #print right
        output_A = dict()
        output_B = dict()
        outfile.write('Prime: ' + uni_prime_A + '\n' + '\n')

        for rec in right:
                outfile.write(rec[0])
                output_A[rec[0]] = []
                for y in infile:
                        x = str(y)
                        domain = x.split(':')[0]
                        prot = x.split(':')[1].split(' ')

                        dist = rec[1].get(str(domain.split(' ')[0].split('\n')[0]))
                        if (uni_prime_A in prot):
                                outfile.write(domain + ' - ' + str(len(prot)) + ' matches, distance = ' + str(dist) + '\n')
                        outfile.write('\n')
                        output_A[rec[0]] += [domain + ' - ' + str(len(prot)) + ' matches, distance = ' + str(dist) + '\n']



        infile = open(domain_B, 'r').readlines()
        outfile = open(distance_out, 'a')

        outfile.write('Prime: ' + uni_prime_B + '\n' + 'Domains:\n')

        for rec in right:
                outfile.write(rec[0])
                output_B[rec[0]] = []
                for y in infile:
                        x = str(y)
                        domain = x.split(':')[0]
                        prot = x.split(':')[1].split(' ')

                        dist = rec[1].get(str(domain.split(' ')[0].split('\n')[0]))
                        if (uni_prime_B in prot):
                                outfile.write(domain + ' - ' + str(len(prot)) + ' matches, distance = ' + str(dist) + '\n')
                        outfile.write('\n')
                        output_B[rec[0]] += [domain + ' - ' + str(len(prot)) + ' matches, distance = ' + str(dist) + '\n']
        return output_A, output_B


def output_binder (domains_A, domains_B, dist_A, dist_B, prime_A, prime_B, domain_file_A, domain_file_B):

        #out_A = []
        for i in dist_A:
                print ('complex: ' + i)
                out = []
                for j in dist_A[i]:
                        #print j
                        domain = j.split(' ')[0]
                        match = j.split('-')[1].split('matches')[0].split(']')[0].split('[')[0]
                        print ['here', match]
                        dist = j.split('[')#[1].split('\n')[0].split(',')
                        if (len(dist) == 1):
                                dist = None
                                for k in domains_A:
                                        if (k[0] == domain):
                                                dist = [k[1], k[2], None]
                        else:
                                dist = j.split('[')[1].split('\n')[0].split(',')
                                dist = [dist[0].split(' ')[0], dist[0].split(' ')[1], dist[1].split(' ')[1]]
                        #print dist

                        if (dist != None):
                        #        print ('        domain: ' + domain)
                        #        print ('        matches: ' + match)
                        #        print ('        amin_num: ' + dist[2])
                        #        print ('        good atoms: ' + dist[0] + ' of ' + dist[1])
                                out.append([domain, int(match.split(']')[0]), dist[2], dist[0], dist[1]])
                out = np.array(out)
                out = out[np.argsort(out[:,1])]
                print out
                uni_set = set()
                file_A = open(domain_file_A, 'r').readlines()
                for rec in file_A:
                        x = rec.split(':')[1].split(' ')
                        for j in x:
                                uni_set.update([j])
                print ('all matches: ' + str(len(uni_set)))
                seq_len = len(open(prime_A, 'r').readlines()[1])
                print ('all sequence: ' + str(seq_len))



        #out_A = []
        for i in dist_B:
                print ('complex: ' + i)
                out = []
                for j in dist_B[i]:
                        #print j
                        domain = j.split(' ')[0]
                        match = j.split('-')[1].split('matches')[0].split(']')[0].split('[')[0]
                        print ['here', match]
                        dist = j.split('[')#[1].split('\n')[0].split(',')
                        if (len(dist) == 1):
                                dist = None
                                for k in domains_B:
                                        if (k[0] == domain):
                                                dist = [k[1], k[2], None]
                        else:
                                dist = j.split('[')[1].split('\n')[0].split(',')
                                dist = [dist[0].split(' ')[0], dist[0].split(' ')[1], dist[1].split(' ')[1]]
                        #print dist

                        if (dist != None):
                        #       print ('        domain: ' + domain)
                        #       print ('        matches: ' + match)
                        #       print ('        amin_num: ' + dist[2])
                        #       print ('        good atoms: ' + dist[0] + ' of ' + dist[1])
                                out.append([domain, int(match.split(']')[0]), dist[2], dist[0], dist[1]])
                out = np.array(out)
                out = out[np.argsort(out[:,1])]
                print out
                uni_set = set()
                file_B = open(domain_file_B, 'r').readlines()
                for rec in file_B:
                        x = rec.split(':')[1].split(' ')
                        for j in x:
                                uni_set.update([j])
                print ('all matches: ' + str(len(uni_set)))
                seq_len = len(open(prime_B, 'r').readlines()[1])
                print ('all sequence: ' + str(seq_len))


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

        run_hhblits (args.protein_A, args.fout_align_prot_A, args.fout_hhr_prot_A)
        run_hhblits (args.protein_B, args.fout_align_prot_B, args.fout_hhr_prot_B)

        print 'START ID MAPPING'

        prime_A, pdb_ids_A = get_pdb_ids (args.fout_align_prot_A)
        biogrid_ids_A = get_biogrid_chain_ids (args.fout_align_prot_A)
        string_ids_A = get_string_ids_sql (args.fout_align_prot_A)

        prime_B, pdb_ids_B = get_pdb_ids (args.fout_align_prot_B)
        biogrid_ids_B = get_biogrid_chain_ids (args.fout_align_prot_B)
        string_ids_B = get_string_ids_sql (args.fout_align_prot_B)

        ## TO KEEP
        to_keep_a = biogrid_ids_A
        to_keep_b = biogrid_ids_B

        ## Find interactions

        if (len (biogrid_ids_A) > len (biogrid_ids_B)):
                biogrid_interactions = find_biogrid_interactions (biogrid_ids_B, biogrid_ids_A)
                string_interactions = find_string_interactions (string_ids_B, string_ids_A)
        else:
                biogrid_interactions = find_biogrid_interactions (biogrid_ids_A, biogrid_ids_B)
                string_interactions = find_string_interactions (string_ids_A, string_ids_B)

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
        domains_A_co, domains_B_co, domains_co = find_interacting_domains (args.fout_align_prot_A, args.fout_align_prot_B, args.fout_complex, prime_A, prime_B, pdb_ids_A, pdb_ids_B, args.fout_distance)
        domains_A_contra, domains_B_contra, domains_contra = find_interacting_domains (args.fout_align_prot_B, args.fout_align_prot_A, args.fout_complex, prime_B, prime_A, pdb_ids_B, pdb_ids_A, args.fout_distance)
        print domains_A_co, domains_A_contra
        print domains_B_co, domains_B_contra
        dist_A, dist_B = prime_domains (prime_A, prime_B, args.fout_domain_A, args.fout_domain_B, domains_co, domains_contra, args.fout_distance)


        print dist_A, dist_B
        output_binder (domains_A_co, domains_contra, dist_A, dist_B, args.protein_A, args.protein_B, args.fout_domain_A, args.fout_domain_B)

if __name__ == '__main__':
        main()


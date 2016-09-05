import os
import argparse
import numpy as np
import pandas as pd
from Bio import SeqIO
from multiprocessing import Pool
from Bio import SeqIO
import ConfigParser

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
#print score + 1
#k = input()

def run_biogrid_finder (gene):
    
    ## Find biogrid interactions using finder
        
    query_string = './finder.o ' + 'B ' + biogrid + ' ' + gene.split('_')[0]
    #print query_string
    os.system(query_string)


def run_string_finder (gene):

    ## Find string interactions using finder    

    query_string = './finder.o ' + 'S ' + string + ' ' + gene
    #print query_string
    os.system(query_string)

def run_string_finder_sql (gene):
    
    ## Find string interactions using SQL finder

    query_string = 'python find_sqlite_db_paired.py ' + gene.split('\n')[0]
    os.system(query_string)


def run_hhblits (sequence, ident_file, align_file):
    
    ## Find similar proteins

    query_string = 'hhblits -i ' + sequence + ' -d ' + uniprot20 + ' -d ' + pdb100 + ' -oa3m ' + ident_file + ' -cpu ' + str(cpu) + '-qid ' + str(score) + '-id 0 -v 0' + ' -o ' + align_file 
    #print query_string
    os.system (query_string)    

    ## Create alignment

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

    return pdb_ids



def get_biogrid_chain_ids (infile):

        ## The same as get_biogrid_ids but with another db for mapping

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

        ## TO DO: 1) pdb to uni_id, 2) uni_id to biogrid_gene (to check: with _HUMAN or without)        
        buf_ids = []
        for x in pdb_ids:
                os.system ('python get_chain_uni_db.py ' + str(x))
                tmp = open('../tmp/chain_uni_out.txt', 'r')
                for i in tmp.readlines():
                        buf_ids.update([i.split('\n')[0]])

        for x in buf_ids:
                os.system ('python get_uni_id.py ' + str(x))
                tmp = open('../tmp/uni_id_out.txt', 'r')
                for i in tmp.readlines():
                        uniprot_ids.update([i.split('\n')[0]])

        uniprot_ids = list(set(list(uniprot_ids) + list(pdb_to_uniprot (pdb_ids))))

        return uniprot_ids



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

    uniprot_ids = get_biogrid_ids (infile)

    string_ids = set()
    
    for x in uniprot_ids:
        os.system ('python get_full_uni_db.py ' + str(x))
        tmp = open('../tmp/full_out.txt', 'r')
        for i in tmp:
            string_ids.update([i.split('\n')[0]])

    for x in uniprot_ids:
        os.system ('python get_all_uni_db.py ' + str(x))
        tmp = open('../tmp/all_out.txt', 'r')
        for i in tmp:
            string_ids.update([i.split('\n')[0]])

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

        buf = d.setdefault(x.split('_')[0])
        if (buf != None):
            string_ids.update(buf)

    return string_ids


def find_biogrid_interactions (ids_A, ids_B):
    
    ## Find biogrid interactions using run_biogrid_finder in 16 threads

    interactions = set()
    pool = Pool(16)

    pool.map(run_biogrid_finder, list(ids_A))

    for gene in ids_A:
        output = open ('../tmp/' + gene.split('_')[0], 'r')
        for entry in output:
            for gene_B in ids_B:
                if (entry.translate(None, '\n') == gene_B.split('_')[0]):
                    interactions.update([str(gene.split('_')[0] + ' ' + gene_B.split('_')[0] + '\n')])


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
            
            os.system ('python get_all_uni_db.py ' + str(prot1))
            tmp = open('../tmp/all_out.txt', 'r')
            for i in tmp:
                pat1.update([i.split('\n')[0]])

            ## Map second protein

            os.system ('python get_full_uni_db.py ' + str(prot2))
            tmp = open('../tmp/full_out.txt', 'r')
            for i in tmp:
                pat2.update([i.split('\n')[0]])
            
            os.system ('python get_all_uni_db.py ' + str(prot2))
            tmp = open('../tmp/all_out.txt', 'r')
            for i in tmp:
                pat2.update([i.split('\n')[0]])

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
            #    d[data[i][1].split('|')[1]] = [data[i][2]]
            #else:
            #    buf.append(data[i][2])
            #    d[data[i][1].split('|')[1]] = buf


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
        
        #print 'kewk'

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

            os.system ('python get_all_uni_db.py ' + str(prot1))
            tmp = open('../tmp/all_out.txt', 'r')
            for i in tmp:
                pat1.update([i.split('\n')[0]])

            ## Map second protein

            os.system ('python get_full_uni_db.py ' + str(prot2))
            tmp = open('../tmp/full_out.txt', 'r')
            for i in tmp:
                pat2.update([i.split('\n')[0]])

            os.system ('python get_all_uni_db.py ' + str(prot2))
            tmp = open('../tmp/all_out.txt', 'r')
            for i in tmp:
                pat2.update([i.split('\n')[0]])

            for p1 in pat1:
                for p2 in pat2:
                    uniprot_interactions.update([str(p1) + ' ' + str(p2) + '\n'])

            #print x
            #print pat1
            #print pat2
            #k = input()

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
            #    d[data[i][1].split('|')[1]] = [data[i][2]]
            #else:
            #    buf.append(data[i][2])
            #    d[data[i][1].split('|')[1]] = buf


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


def main ():

    #print "START"

    ## Input sequences as a command line args
    
    parser = argparse.ArgumentParser (description='Finding protein interactions')
    parser.add_argument ('protein_A', metavar='prot_A', help='file with sequence A')
    parser.add_argument ('protein_B', metavar='prot_B', help='file with sequence B')
    parser.add_argument ('fout_biogrid', metavar='pdb_out', help='output file for pdb interactions')
    parser.add_argument ('fout_string', metavar='uniprot_out', help='output file for uniprot interactions')
    parser.add_argument ('fout_complex', metavar='complex_out', help='output file for complexes')
    parser.add_argument ('fout_align_prot_A', metavar='align_out', help='output file for multiple sequence alignment for sequence A')
    parser.add_argument ('fout_hhr_prot_A', metavar='hhr_out', help='output file for hhalign result for sequence A')
    parser.add_argument ('fout_align_prot_B', metavar='align_out', help='output file for multiple sequence alignment for sequence B')
    parser.add_argument ('fout_hhr_prot_B', metavar='hhr_out', help='output file for hhalign result for sequence B')
    

    args = parser.parse_args ()

    #print args.protein_A
    #print args.protein_B

    ## Parse config file

    #config = open (args.conf, 'r')
    #for x in config.readlines():
    #x = config.readlines()
    
    #pdbtosp = x[0].split('|')[1]
    #biogrid = x[1].split('|')[1]
    #uniprot = x[2].split('|')[1]
    #string  = x[3].split('|')[1]
    #uniprot20 = x[4].split('|')[1]
    #pdb100 = x[5].split('|')[1]
    #full_uniprot = x[6].split('|')[1]
    #all_knowledge = x[7].split('|')[1]
    
    #config_parser (args.conf)
    #print full_uniprot
    #k = input()

    print "START HHBLITS"

    ## Run hhblits for the given sequences

    #run_hhblits (args.protein_A, '../tmp/' + args.protein_A.split('/')[2].split('.')[0]+ '_identies.txt', '../tmp/' + args.protein_A.split('/')[2].split('.')[0] + '_multi_align.txt')
    #run_hhblits (args.protein_B, '../tmp/' + args.protein_B.split('/')[2].split('.')[0] + '_identies.txt', '../tmp/' + args.protein_B.split('/')[2].split('.')[0] + '_multi_align.txt')
    
    run_hhblits (args.protein_A, args.fout_align_prot_A, args.fout_hhr_prot_A)
    run_hhblits (args.protein_B, args.fout_align_prot_B, args.fout_hhr_prot_B)

    print 'START ID MAPPING'

    ## ID mapping
    #print '../tmp/' + args.protein_A.split('/')[2].split('.')[0] + '_identies.txt'
    #biogrid_ids_A = get_biogrid_ids ('../tmp/' + args.protein_A.split('/')[2].split('.')[0] + '_identies.txt')
    #string_ids_A = get_string_ids ('../tmp/' + args.protein_A.split('/')[2].split('.')[0] + '_identies.txt')
    #print 'OTHER'
    #print '../tmp/' + args.protein_B.split('/')[2].split('.')[0] + '_identies.txt'
    #biogrid_ids_B = get_biogrid_ids ('../tmp/' + args.protein_B.split('/')[2].split('.')[0] + '_identies.txt')
    #string_ids_B = get_string_ids ('../tmp/' + args.protein_B.split('/')[2].split('.')[0] + '_identies.txt')
    #print full_uniprot
    biogrid_ids_A = get_biogrid_ids (args.fout_align_prot_A)
    string_ids_A = get_string_ids_sql (args.fout_align_prot_A)

    biogrid_ids_B = get_biogrid_ids (args.fout_align_prot_B)
    string_ids_B = get_string_ids_sql (args.fout_align_prot_B)

    print "FINDING INTERACTIONS"
    #k = input()
    #print biogrid_ids_A
    #print biogrid_ids_B

    #print string_ids_A
    #print string_ids_B

    #k = input()
    ## Find interactions

    if (len (biogrid_ids_A) > len (biogrid_ids_B)):
        biogrid_interactions = find_biogrid_interactions (biogrid_ids_B, biogrid_ids_A)
        string_interactions = find_string_interactions (string_ids_B, string_ids_A)
    else:
        biogrid_interactions = find_biogrid_interactions (biogrid_ids_A, biogrid_ids_B)
        string_interactions = find_string_interactions (string_ids_A, string_ids_B)

    #print biogrid_interactions
    #print string_interactions    

    ## Get pdb and uniprot ID's for interacting proteins

    pdb_interactions = find_pdb_id_sql (biogrid_interactions, 'B') + find_pdb_id_sql (string_interactions, 'S')
    uniprot_interactions = biogrid_interactions + find_uniprot_id_sql (string_interactions, 'S')

    #print pdb_interactions
    #print uniprot_interactions

    ## Find pdb complex ID's

    ## TO DO
    #pdb_complexes = find_pdb_complex_id (pdb_interactions)

    ## Run hhpred
    ## TO DO

    ## Write down the results

    output = open (args.fout_biogrid, 'w')
    #output.write (str(pdb_interactions))
    for x in pdb_interactions:
        output.write(str(x))

    output = open (args.fout_string, 'w')
    #output.write (uniprot_interactions)

    for x in uniprot_interactions:
        output.write(str(x))

    #Finding complexes
    complexes = []
    for x in pdb_interactions:
        if (x.split(' ')[0] == x.split(' ')[1].split('\n')[0]):
            complexes.append(x.split(' ')[0] + '\n')

    output = open (args.fout_complex, 'w')
    #output.write (complexes)
    for x in complexes:
        output.write(str(x))


if __name__ == '__main__':
    main()


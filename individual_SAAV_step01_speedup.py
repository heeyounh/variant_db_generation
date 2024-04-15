# -*- coding: utf-8 -*-
"""
Created on Tue Jan 15 09:40:11 2019

@author: Heeyoun
"""
import os
import pickle
#from numba import jit



#@autojit
def search(path):
    res = []
    for root, dirs, files in os.walk(path):
        rootpath = os.path.join(os.path.abspath(path), root)
        
        for file in files:
            filepath = os.path.join(rootpath, file)
            
            res.append(filepath)
    return res

#@jit(nopython=True, cache=True)
def vcf_file_cook (file_):
        vcf_file = open(file_)
        vcf_file_lines = vcf_file.readlines()
        VCF_DIC = []
        result_file = open(file_[:-4] + '.txt', 'w')
        dmp_file = open(file_[:-4] + '.dmp', 'wb')
        for line in vcf_file_lines:
            if line[0] == '#':
                pass
            else:
                lines = line[:-1].split('\t')
                title = ';'.join(lines[:5])
                Chr = lines[0]
                try:
                    ENST_list = CHR_ENST_dic[Chr]
                    for ENST in ENST_list:
                        GTF_ = GTF_fasta_seq[ENST]
                        Exon_list = GTF_[6].split(';')
                        for EXON_ in Exon_list:
                            EXON = EXON_.split(':')
                            #print  int(EXON[0]), lines[1], int(EXON[1])
                            if int(EXON[0]) <= int(lines [1]) <=int(EXON[1]):
                                VCF_DIC.append(lines[:5] + GTF_)
                                print(file_ + '\t', end = '', file= result_file)
                                for line in lines[:5]:
                                    print( line + '\t', end = '', file= result_file)
                                for GTF in GTF_:
                                    print( GTF + '\t', end = '', file = result_file)
                                print( '', file= result_file)
                except:
                    pass
        pickle.dump(VCF_DIC, dmp_file)
        result_file.close()
        dmp_file.close()
        vcf_file.close()
    

nowdir = os.getcwd()

file_list = search(nowdir)
GTF_file = open('GTF_transcript_sequence.dmp', 'rb')
GTF_fasta_seq = pickle.load(GTF_file)
GTF_file.close()

GTF_ENST_file = open('GTF_CHR_ENST.dmp', 'rb')
CHR_ENST = pickle.load(GTF_ENST_file)
GTF_ENST_file.close()
CHR_ENST_dic = {}
for chr_ in CHR_ENST.keys():
    CHR_ENST_list = CHR_ENST[chr_]
    ENST_list = []
    for list_ in CHR_ENST_list:
        list__ = list_[0].split(';')
        ENST_list.append(list__[2])
    CHR_ENST_dic [chr_] = ENST_list

print('dmp loading is over.')



for file_ in file_list:
    if ".vcf" == file_[-4:]:
        print (file_)
        vcf_file_cook(file_)
        
                #VCF_DICprint title
        
            
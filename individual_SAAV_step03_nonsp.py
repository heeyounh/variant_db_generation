# -*- coding: utf-8 -*-
"""
Created on Mon Jan 21 19:30:32 2019

@author: Heeyoun
"""
import os
import pickle
import re

def search(path):
    res = []
    for root, dirs, files in os.walk(path):
        rootpath = os.path.join(os.path.abspath(path), root)
        
        for file in files:
            filepath = os.path.join(rootpath, file)
            
            res.append(filepath)
    return res


nowdir = os.getcwd()

total_gene_dic = {}
gene_file_dic = {}
file_list = search(nowdir)
all_files = []
Fasta_file = open('SAAV3.fasta', 'w')
for file_ in file_list:
    if "_SAAV.dmp2" in file_:
        print (file_ + ' is loaded.')
        dmp_file = open(file_, 'rb')
        dmp_file_dic = pickle.load(dmp_file)

        file_names = file_[:-5].split('\\')
        for key in dmp_file_dic:
            data_list = dmp_file_dic[key]
            keys_ = key.split('|')
            keys = keys_[-1].split(':')
            print('>TCONS_' + key + '_ counts:'+ str(len(data_list)), file = Fasta_file)
            print(keys[1], file = Fasta_file)
#             title_list = data_list[0]
#             SAAV_list = data_list[2]
#             for SAAV in SAAV_list:
#                 try:
#                     total_gene_dic [title_list[1]] = total_gene_dic [title_list[1]] + [SAAV]
#                     gene_file_dic [SAAV] = gene_file_dic [SAAV] + [file_names[-2] + '|' + file_names[-1]]
#                 except:
#                     total_gene_dic [title_list[1]] = [SAAV]
#                     gene_file_dic [SAAV] = [file_names[-2] + '|' + file_names[-1]]
#                 all_files.append(file_names[-2] + '|' + file_names[-1])    
       
# genes = list(total_gene_dic.keys())
# genes.sort()             
# all_files = list(set(all_files))
# all_files.sort()
                    

# result_File = open('Fasta_txt2.txt', 'w')
# result_File2 = open('Fasta_txt2_1.txt', 'w')
# print ( '\t\t\t\t', end = '', file = result_File2)
# for file__ in all_files:
#     print ( file__ + '\t', end = '', file = result_File2)
# print ( '', file = result_File2)


# for gene in genes:
#     pep_list = total_gene_dic [gene]
#     pep_list = list(set(pep_list))

#     #file_list2 = list(set(file_list2))
#     print(">TCONS_" + gene + '|' + str(len(pep_list)) + '_SAAV peptides', file = Fasta_file)
#     print(gene  + '\t' + str(len(pep_list)) + '\t' + ';'.join(pep_list), file = result_File)
#     for pep_ in pep_list:
#         try:
#             file_list2 = gene_file_dic [pep_]
#         except:
#             file_list2 = []
#         print(gene + '\t' + str(len(pep_list)) + '\t' + pep_ + '\t' + str(len(file_list2)) + '\t', end = '', file = result_File2)
#         for file__ in all_files:
#             if file__ in file_list2:
#                 print('1\t', end = '', file = result_File2)
#             else:
#                 print('0\t', end = '', file = result_File2)
#             print('', file = result_File2)
#     fasta = ''
#     for pep_ in pep_list:
#         fasta = fasta + pep_
    
#     for i in range(0, len(fasta), 60):
#         print(fasta[i:i + 60], file = Fasta_file)

#genes = gene_file_dic.keys()
#genes.sort()
#
#for gene in genes:
#    file_list = gene_file_dic [gene]
#    file_list = list(set(file_list))
    #print>>result_File, gene  + '\t' + str(len(pep_list)) + '\t' + ';'.join(pep_list)  + '\t' + str(len(file_list)) + '\t' + ';'.join(file_list)
    

Fasta_file.close()
# result_File.close()
# result_File2.close()

        
                
        
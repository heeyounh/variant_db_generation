# -*- coding: utf-8 -*-
"""
Created on Tue Dec 18 15:10:14 2018

@author: Heeyoun
"""
import pickle

GENCODE_fasta_file = open('GRCm38.p4.genome.fa')
GENCODE_file_lines = GENCODE_fasta_file.readlines()
GENCODE_fasta_file.close()

print ('fasta file is opened')

# GENCODE_fasta_file2 = open('GRCh37.p13.genome.fasta', 'w')
#test = [16916394,16916496]
#test_chr = 1

title_index_dic = {}
for i in range(len(GENCODE_file_lines)):
    if GENCODE_file_lines[i][0] == ">":
        lines = GENCODE_file_lines[i][:-1].split(' ')
        title = lines[0][1:]
        title_index_dic [title] = i



# GENCODE_fasta_dic = {}
# titles = []
# for line in GENCODE_file_lines:
#     if '>' in line[0]:
#         title_ = line[:-1].split(' ')
#         title = title_[0][1:]
#         titles.append(title)
#         # print (line[:-1], file =GENCODE_fasta_file2)
#     #elif '' in line[0]:
#     #    pass
#     else:
#         try:
#             GENCODE_fasta_dic [title] = GENCODE_fasta_dic [title] + [line[:-1]]
#         except:
#             GENCODE_fasta_dic [title] = [line[:-1]]
#         # try:
        #     for i in range(0, len(line[:-1]), 60):
        #         print(line[i:i + 59], file= GENCODE_fasta_file2)
                
        # except:
        #     print (line[i:], file =GENCODE_fasta_file2)
            
# GENCODE_fasta_file2.close()

# print ('fasta file is loaded.')
#GENCODE_fasta_file3 = open('test.fasta', 'w')
#chrY = GENCODE_fasta_dic ['>chrY Y']
#
#try:
#    for i in range(0, len(chrY[:-1]), 60):
#        print>>GENCODE_fasta_file3, chrY[i:i + 59]
#        
#except:
#    print>>GENCODE_fasta_file3, chrY[i:]
#
#print>>GENCODE_fasta_file3, chrY[13324079-1], 'T -> C', chrY[13324079-2:13324079]
#print>>GENCODE_fasta_file3, chrY[13324118-1], 'G -> A', chrY[13324118-2:13324118]
#print>>GENCODE_fasta_file3, chrY[13468285-1], 'T -> G', chrY[13468285-2:13468285]
#
#GENCODE_fasta_file3.close()
GTF_file = open('gencode.vM10.annotation.gtf')
GTF_file_line = GTF_file.readlines()
GTF_file.close()
CHR_GTF_dic = {}
transcript_strand_dic = {}
transcript_exon_dic = {}
CHR_list = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chrX', 'chrY', 'chrM']

for line in GTF_file_line:
    if line[0] == '#':
        pass
    else:
        lines = line.split('\t')
        CHR = lines[0]
        
        transcript_list = []
        if lines[2] == 'transcript':
            gene_id = ''
            transcript_id = ''
            gene_type = ''
            gene_name = ''
            strand = lines[6]
            infos = lines[-1].split(';')
            for i in range(0, len(infos)):
                info = infos[i].split('"')
                if 'gene_id' in infos[i]:
                    gene_id = info[1]
                elif 'transcript_id' in infos[i]:
                    transcript_id = info[1]
                elif 'gene_type' in infos[i]:
                    gene_type = info[1]
                elif 'gene_name' in infos[i]:
                    gene_name = info[1]
            transcript_list.append(';'.join([gene_name,gene_id,transcript_id,gene_type]))
            transcript_strand_dic [transcript_id] = strand 
            try:
                CHR_GTF_dic [CHR] = CHR_GTF_dic [CHR]  + [transcript_list] 
            except:
                CHR_GTF_dic [CHR] = [transcript_list]
        elif lines[2] == 'exon':
            gene_id = ''
            transcript_id = ''
            gene_type = ''
            gene_name = ''
            start = lines[3]
            end = lines[4]
            #strand = lines[6]
            
            infos = lines[-1].split(';')
            for i in range(0, len(infos)):
                info = infos[i].split('"')
                if 'gene_id' in infos[i]:
                    gene_id = info[1]
                elif 'transcript_id' in infos[i]:
                    transcript_id = info[1]
                elif 'gene_type' in infos[i]:
                    gene_type = info[1]
                elif 'gene_name' in infos[i]:
                    gene_name = info[1]
            try:
                transcript_exon_dic [transcript_id] = transcript_exon_dic [transcript_id] + [':'.join([start,end])]
            except:
                transcript_exon_dic [transcript_id] = [':'.join([start,end])]                
transcript_dic = {}                            
GTF_seq_file = open('GTF_transcript_sequence.txt', 'w')            
for CHR_ in CHR_list:
    #Fasta_seq = GENCODE_fasta_dic [CHR_]
    
    CHR_index = title_index_dic [CHR_] + 1
    
    for CHR_GTF in CHR_GTF_dic [CHR_]:
        for GTF_ in CHR_GTF:
            GTF_lines = GTF_.split(';')
            transcript_ = GTF_lines[2]
            transcript_exon_list = transcript_exon_dic [transcript_]
            fastas = ''
            for exons in transcript_exon_list:
                exon = exons.split(':')
                
                exon1 = int(exon[0])
                exon2 = int(exon[1])
                
                exon_line_index1 = exon1 // 60
                exon_line_index2 = exon2 // 60
                
                exon_site_index1 = exon1 % 60 - 1
                exon_site_index2 = exon2 % 60 - 1
                
                
                exon_seq = ''
                for i in range(exon_line_index1, exon_line_index2 + 1):
                    if i == exon_line_index1:
                        exon_seq = exon_seq + GENCODE_file_lines[CHR_index + exon_line_index1][exon_site_index1:-1]
                    
                    elif i == exon_line_index2:
                        exon_seq = exon_seq + GENCODE_file_lines[CHR_index + exon_line_index2][:exon_site_index2 + 1]
                    else:
                        
                        exon_seq = exon_seq + GENCODE_file_lines[CHR_index + i][:-1]
                
                
                
                
                
                #fasta = Fasta_seq [int(exon[0])-1:int(exon[1])]
                fastas = fastas + exon_seq
            print(CHR_ + '\t' + GTF_lines[0] + '\t' + GTF_lines[1] + '\t' + GTF_lines[2] + '\t' + GTF_lines[3] + '\t' + transcript_strand_dic[transcript_] + '\t' + ';'.join(transcript_exon_list) + '\t' + str(len(transcript_exon_list)) + '\t' + str(len(fastas)) + '\t' + fastas, file =  GTF_seq_file)
            transcript_dic [GTF_lines[2]] = [CHR_, GTF_lines[0], GTF_lines[1], GTF_lines[2], GTF_lines[3], transcript_strand_dic[transcript_], ';'.join(transcript_exon_list), str(len(transcript_exon_list)), str(len(fastas)), fastas]
GTF_seq_file.close()
print ('dumping...')
GTF_seq_dmp_file = open('GTF_transcript_sequence.dmp', 'wb')
pickle.dump(transcript_dic, GTF_seq_dmp_file)
CHR_ENST_file = open('GTF_CHR_ENST.dmp', 'wb')
pickle.dump(CHR_GTF_dic, CHR_ENST_file)
CHR_ENST_file.close()
GTF_seq_dmp_file.close()

        
    
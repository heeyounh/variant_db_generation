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

def REVERSE(seq):
    Seq = seq[::-1]
    Seq_1 = Seq.replace('A', 't')
    Seq_2 = Seq_1.replace ('G', 'c')
    Seq_3 = Seq_2.replace ('C', 'g')
    Seq_4 = Seq_3.replace ('T', 'a')
    return Seq_4.upper()

def Translator (ori_seq, start, end):
    

    seqfinal = ''
    Sequence = ori_seq[start-1 : end]
    CodonList = []
    for x in range(0, len(Sequence), 3):
        CodonList.append(Sequence[x:x+3])
        
    seqfinal = ''
    for codon in CodonList:
        if codon in CodonDict:
            seqfinal = seqfinal + CodonDict[codon]
            
        else:
            break
    seqlist = seqfinal.split("X")
    print(VCF_[:10], Sequence, seqlist[0], file = syno_file)
    
    return seqlist[0]
        
    
    
    
    
    
def Translator2(Sequence, SNV_Sequence):
    seqfinal2 = ''
    #translation_frames = []
    seq1 = []
    seqfinal = ''
    for frame in range(0,3):
        #--- Translate sequence into a list of codons
        CodonList = [ ]
        for x in range(frame, len(Sequence), 3):
            CodonList.append(Sequence[x:x+3])       
        #--- For each codon, translate it into amino acid, and create a list
        ProteinSeq = ''
        for codon in CodonList:
            if codon in CodonDict:
                ProteinSeq = ProteinSeq + CodonDict[codon]
                
            else:
                break

        
        #seqtemp = ''.join(ProteinSeq)
        seqlist = ProteinSeq.split("X")
        #print(frame, seqlist)
        for i in seqlist:
            #print i
            k = i.find('M')
            if k > -1:
                seq0 = i[k:]
                #print seq0
                if len(seq0) >= 30:
                    seq1.append(seq0)
    print(VCF_[:10], seq1, file = syno_file)
    length = []
    #print(seq1)
    for l in seq1:
        #print l, len(l)
        length.append(len(l))
    index = length.index(max(length))
    seqfinal = seq1[index]
    # for i in range(len(seq1)):
        
    #     if len(seq1[i]) == max(length):
    #         seqfinal = seq1[i]
    #         index = i
    if seq1 == []:
        pass
    else:
        seq2 = []
        #seqfinal2 = ''
        for frame2 in range(0,3):
            #--- Translate sequence into a list of codons
            CodonList2 = [ ]
            for x2 in range(frame2, len(SNV_Sequence), 3):
                CodonList2.append(SNV_Sequence[x2:x2+3])       
            #--- For each codon, translate it into amino acid, and create a list
            ProteinSeq2 = ''
            for codon2 in CodonList2:
                if codon2 in CodonDict:
                    ProteinSeq2 = ProteinSeq2 + CodonDict[codon2]
                    
                else:
                    break
    
            
            #seqtemp = ''.join(ProteinSeq)
            seqlist2 = ProteinSeq2.split("X")
            for i2 in seqlist2:
                #print i
                k2 = i2.find('M')
                if k2 > -1:
                    seq02 = i2[k2:]
                    #print seq0
                    if len(seq02) >= 30:
                        seq2.append(seq02)
        print(VCF_[:10], seq2, file = syno_file)
        length2 = []
        #print(seq1)
        for l2 in seq2:
            #print l, len(l)
            length2.append(len(l2))
        index2 = length2.index(max(length2))

    # length2 = []
    # for l in seq2:
    #     #print l, len(l)
    #     length2.append(len(l))
    
    # for m in seq1:
    #     #print m
    #     if len(m) == max(length):
    #         seqfinal = m
    #print(seq1, seq2)
        seqfinal2 = seq2 [index2]
        # try:    
        #     seqfinal2 = seq2 [index]
        # except:
        #     pass
        #     #print (seqfinal, index, seq2)
    
    #temp = ((frame+1), seqfinal)
    #translation_frames.append(temp)
    return seqfinal, seqfinal2

def T_digestion(protein_seq):
    K_dig = re.compile('K')
    R_dig = re.compile('R')
    PK_dig = re.compile('K P')
    RK_dig = re.compile('R P')
    seq_1 = K_dig.sub('K ', protein_seq)
    seq_2 = R_dig.sub('R ', seq_1)
    seq_3 = PK_dig.sub('KP', seq_2)
    seq_4 = RK_dig.sub('RP', seq_3)
    peplist = seq_4.split(' ')
    return peplist

def slicing(Ori_pro_seq, target_list, MIN, MAX):
    min_target = min(target_list)
    max_target = max(target_list)
    min_seq = ''
    max_seq = ''
    
    if min_target - MIN < 0:
        min_seq = Ori_pro_seq[:min_target]
        pass
    else:
        min_seq = Ori_pro_seq[min_target - MIN: min_target]
    
    if max_target + MAX > len(Ori_pro_seq):
        max_seq = Ori_pro_seq[max_target:]
    else:
        max_seq = Ori_pro_seq[max_target:max_target + MAX]
        
    return min_seq + max_seq
        
    
    
    


CodonDict={'ATT':'I',   'ATC':'I',  'ATA':'I',  'CTT':'L',  'CTC':'L',  
'CTA':'L',  'CTG':'L',  'TTA':'L',  'TTG':'L',  'GTT':'V',  'GTC':'V',  
'GTA':'V',  'GTG':'V',  'TTT':'F',  'TTC':'F',  'ATG':'M',  'TGT':'C',  
'TGC':'C',  'GCT':'A',  'GCC':'A',  'GCA':'A',  'GCG':'A',  'GGT':'G',  
'GGC':'G',  'GGA':'G',  'GGG':'G',  'CCT':'P',  'CCC':'P',  'CCA':'P',  
'CCG':'P',  'ACT':'T',  'ACC':'T',  'ACA':'T',  'ACG':'T',  'TCT':'S',  
'TCC':'S',  'TCA':'S',  'TCG':'S',  'AGT':'S',  'AGC':'S',  'TAT':'Y',  
'TAC':'Y',  'TGG':'W',  'CAA':'Q',  'CAG':'Q',  'AAT':'N',  'AAC':'N',  
'CAT':'H',  'CAC':'H',  'GAA':'E',  'GAG':'E',  'GAT':'D',  'GAC':'D',  
'AAA':'K',  'AAG':'K',  'CGT':'R',  'CGC':'R',  'CGA':'R',  'CGG':'R',  
'AGA':'R',  'AGG':'R',  'TAA':'X',  'TAG':'X',  'TGA':'X'}



nowdir = os.getcwd()

Gencode_file = open('gencode.vM10.pc_transcripts.fa')
file_lines = Gencode_file.readlines()
gencode_dic = {}
for line in file_lines:
    if line[0] == '>':
        lines = line[:-1].split('|')
        for i in range(len(lines)):
            if 'CDS:' in lines[i]:
                lines_cds_ = lines[i].split(':')
                lines_cds = lines_cds_[1].split('-')
                gencode_dic [lines[0][1:]] = [int(lines_cds[0]), int(lines_cds[1])]
                
    

file_list = search(nowdir)
for file_ in file_list:
    if 'GTF_CHR_ENST.dmp' in file_ or 'GTF_transcript_sequence.dmp' in file_:
        pass
    else:
        if file_[-4:] == '.dmp':
            dmp_file = open(file_, 'rb')
            log_file = open(file_[:-4] + '_log_file.txt', 'w')
            dmp2_file = open(file_[:-4] + '_SAAV.dmp2', 'wb')
            syno_file = open(file_[:-4] + 'syno.txt', 'w')
            #vcf_file_lines = vcf_file.readlines()
            VCF_DIC = pickle.load(dmp_file)
            print (file_ + " loading is over.")
            #VCF_keys = list(VCF_DIC.keys())
            #VCF_keys.sort()
            m = 0
            VCF_DIC2 = {}
            for VCF_ in VCF_DIC:
                key1 = ''
                #VCF_ = VCF_DIC [VCF_key]
                VCF_list = VCF_[5:]
                VCF_titles = VCF_[:10]
                #print(VCF_titles)
                site = int(VCF_titles[1])
                ref = VCF_titles[3]
                alt = VCF_titles[4]
                if VCF_list[4] == 'protein_coding':
                    ori_seq = VCF_list[9]
                    exon_list = VCF_list[6].split(';')
                    snv_seq = ''
                    if VCF_list[5] == '-':
                        length = []
                        lengths = []
                        site2 = 0
                        site3 = 0
                        for exon_list_ in exon_list:
                            exon_list___ = exon_list_.split(':')
                            #length.append(int(exon_list___[1]) - int(exon_list___[0]))
    
                            
                            if int(exon_list___[1]) >= site >= int(exon_list___[0]):
                                site3 = site2 + site - int(exon_list___[0])
                                lengths.append(int(exon_list___[1]) - int(exon_list___[0]) - (len(ref) - len(alt)))
                                length.append(int(exon_list___[1]) - int(exon_list___[0]))
                            else:
                                site2 = site2 + int(exon_list___[1]) - int(exon_list___[0]) + 1
                                length.append(int(exon_list___[1]) - int(exon_list___[0]))
                                lengths.append(int(exon_list___[1]) - int(exon_list___[0]))
                                
                        if ori_seq[site3:site3 + len(ref)] == ref:
                            snv_seq = ori_seq[:site3] + alt + ori_seq [site3 + len(ref):]
                            print(VCF_[:10], ori_seq[site3:site3 + len(ref)], ref, alt, site3, ori_seq, snv_seq, file = syno_file)
                        else:
                            print( '|'.join(VCF_[:10]) + '\tref site in error!\t' + ori_seq[site3:site3 + len(ref)] + ref, file = syno_file)
                            key1 = 'pass'
                        n = 0
                        ns = 0
                        reverse_seq_list = []
                        reverse_snv_seq_list = []
                        for i in range(0, len(length)):
                            reverse_seq_list.append(ori_seq[n:n + length[i] + 1])
                            reverse_snv_seq_list.append(snv_seq[ns:ns + lengths[i] + 1])
                            n = n + length[i] + 1
                            ns = ns + lengths[i] + 1
                        print(VCF_[:10], reverse_seq_list, reverse_snv_seq_list, file = syno_file)
                        reverse_seq_2 = []
                        reverse_snv_seq_2 = []
                        for i in range(0, len(reverse_seq_list)):
                            reverse_seq_2.append(REVERSE (reverse_seq_list[i]))
                            reverse_snv_seq_2.append(REVERSE (reverse_snv_seq_list[i]))
                        print(VCF_[:10], reverse_seq_2, reverse_snv_seq_2, file = syno_file)
                        reverse_seq_final = ''
                        reverse_snv_seq_final = ''
                        for i in range(0, len(reverse_seq_2)):
                            reverse_seq_final = reverse_seq_final + reverse_seq_2[i]
                            reverse_snv_seq_final = reverse_snv_seq_final + reverse_snv_seq_2[i]
                            
                        ori_seq = reverse_seq_final
                        snv_seq = reverse_snv_seq_final
                        print(VCF_[:10], ori_seq, snv_seq, file = syno_file)
                    else:
                        length = []
                        lenghs = []
                        site2 = 0
                        site3 = 0
                        for exon_list_ in exon_list:
                            exon_list___ = exon_list_.split(':')
                            
    
                            
                            if int(exon_list___[1]) >= site >= int(exon_list___[0]):
                                site3 = site2 + site - int(exon_list___[0])
                                lengths.append(int(exon_list___[1]) - int(exon_list___[0]) - (len(ref) - len(alt)))
                                length.append(int(exon_list___[1]) - int(exon_list___[0]))
                            else:
                                site2 = site2 + int(exon_list___[1]) - int(exon_list___[0]) + 1
                                length.append(int(exon_list___[1]) - int(exon_list___[0]))
                                lengths.append(int(exon_list___[1]) - int(exon_list___[0]))
                                
                        if ori_seq[site3:site3 + len(ref)] == ref:
                            snv_seq = ori_seq[:site3] + alt + ori_seq [site3 + len(ref):]
                            print(VCF_[:10], ori_seq, snv_seq, file = syno_file)
                        else:
                            key1 = 'pass'
                            print('|'.join(VCF_[:10]) + '\tref site in error!\t' + ori_seq[site3:site3 + len(ref)] + ref, file = syno_file)
                        
                    #print>>log_file, 'DNA_seq', VCF_
                    #if ori_seq == snv_seq:
                        #print>>log_file, 'Same sequence in DNA'
                    #else:
                        
                        #print>>log_file, ori_seq
                        #print>>log_file, snv_seq
                    if key1 == 'pass':
                        pass
                    else:
                        mut_info = ''
                        mut_site = 0
                        # try:
                        if len(ori_seq) <= len(snv_seq):
                            for i in range(0, len(ori_seq)):
                                if ori_seq[i] == snv_seq[i]:
                                    pass
                                else:
                                    mut_info = ori_seq[i: i + len(ref)] + ':' + snv_seq[i: i + len(alt)] + ':' + str(i + 1)
                                    mut_site = i + 1
                                    break
                        else:
                            for i in range(0, len(snv_seq)):
                                if ori_seq[i] == snv_seq[i]:
                                    pass
                                else:
                                    mut_info = ori_seq[i: i + len(ref)] + ':' + snv_seq[i: i + len(alt)] + ':' + str(i + 1)
                                    mut_site = i + 1
                                    break
                            
                        # except:
                        #     pass
                        try:    
                            CDS_list = gencode_dic [VCF_titles[8]]     
                            if mut_site >= CDS_list[0] and mut_site <= CDS_list[1]:
                                #Ori_pro_seq, SNV_pro_seq = Translator (ori_seq, snv_seq, mut_site, CDS_list[0], CDS_list[1])
                                Ori_pro_seq = Translator (ori_seq, CDS_list[0], CDS_list[1])
                                SNV_pro_seq = Translator (snv_seq, CDS_list[0], CDS_list[1])
                                #VCF__ = VCF_.split(';')
                                
                                mut_prot_info = ''
                                target_index = []
                                try:
                                    
                                    for i in range(0, len(SNV_pro_seq)):
                                        if Ori_pro_seq[i] == SNV_pro_seq[i]:
                                            pass
                                        else:
                                            mut_prot_info = Ori_pro_seq[i] + ':' + SNV_pro_seq[i] + ':' + str(i + 1)
                                            target_index.append(i)
                                            break
                                    
                                except:
                                    pass
                                print(CDS_list, mut_site, mut_prot_info, Ori_pro_seq, SNV_pro_seq, file = syno_file)
                                if target_index == []:
                                    if Ori_pro_seq == SNV_pro_seq:
                                        print (file_[:-4] + '\t' + '|'.join(VCF_[:10]) + '\t' + 'Same sequence in Protein', file = syno_file)
                                    elif Ori_pro_seq == '' or SNV_pro_seq == '':
                                        print (file_[:-4] + '\t' + '|'.join(VCF_[:10]) + '\t' + 'only short protein sequences', file = syno_file)
                                    
                                    else:
                                        if SNV_pro_seq in Ori_pro_seq:
                                            i = len(SNV_pro_seq) - 1
                                            target_index.append(i)
                                            mut_prot_info = Ori_pro_seq[i + 1] + ':*:' + str(i + 1)
                                        elif Ori_pro_seq in SNV_pro_seq:
                                            i = len(Ori_pro_seq) - 1
                                            target_index.append(i)
                                            mut_prot_info =  '*:' + SNV_pro_seq[i + 1] + ':' + str(i + 1)
                                        pass
                                
                                
            
                                    
                                    
                                    
                                        
                                    
                                else:
                                    m = m + 1
                                    
                                    #print>>log_file, Ori_pro_seq
                                    #print>>log_file, SNV_pro_seq
                                    
                                    #ref_peps = T_digestion (Ori_pro_seq)
                                    #var_peps = T_digestion (SNV_pro_seq)
                                    
                                    ref_peps = slicing(Ori_pro_seq, target_index, 30, 30)
                                    var_peps = slicing(SNV_pro_seq, target_index, 30, 30)
                                    
                                    #print (Ori_pro_seq, SNV_pro_seq, mut_prot_info, target_index, var_peps, ref_peps)
                                    # if len(target_index) > 1:
                                    #     print (VCF_titles, ori_seq, snv_seq, mut_info)
                                    # v_pep_list = []
                                    # r_pep_list = []
                                    # for v_pep in var_peps:
                                        
                                    #     if v_pep in ref_peps:
                                    #         pass
                                    #     else:
                                    #         if 5 < len(v_pep) < 50:
                                    #             v_pep_list.append(v_pep)
                                    # for r_pep in ref_peps:
                                    #     if r_pep in var_peps:
                                    #         pass
                                    #     else:
                                    #         if 5 < len(r_pep) < 50:
                                    #             r_pep_list.append(r_pep)
                                    if var_peps == '':
                                        for i in range(0, len(VCF_titles)):
                                            print(VCF_titles[i] + '\t', end = '', file = log_file)
                                        for i in range(0, len(VCF_list[:-1])):
                                            print(VCF_list[i] + '\t', end = '', file = log_file)
                                        # print(';'.join(r_pep_list) + '\t' + 'no_detectable_SAAVs\t', file = log_file)
                                        print(ref_peps + '\t' + 'no_detectable_SAAVs\t' + mut_info + '\t'+ mut_prot_info + '\t' + ori_seq + '\t' + snv_seq +'\t' + Ori_pro_seq + '\t' + SNV_pro_seq, file = log_file)
                                        pass
                                    else:
                                        for i in range(0, len(VCF_titles)):
                                            print( VCF_titles[i] + '\t', end = '', file = log_file)
                                        for i in range(0, len(VCF_list[:-1])):
                                            print(VCF_list[i] + '\t', end = '', file = log_file)
                                        # print(';'.join(r_pep_list) + '\t' + ';'.join(v_pep_list) + '\t' + mut_info + '\t'+ mut_prot_info + '\t' + ori_seq + '\t' + snv_seq +'\t' + Ori_pro_seq + '\t' + SNV_pro_seq, file = log_file)
                                        print(ref_peps + '\t' + var_peps + '\t' + mut_info + '\t'+ mut_prot_info + '\t' + ori_seq + '\t' + snv_seq +'\t' + Ori_pro_seq + '\t' + SNV_pro_seq, file = log_file)
                                        #VCF_DIC2 [':'.join(VCF_titles)] = [VCF_list, r_pep_list, v_pep_list, mut_info, mut_prot_info]
                                        try:
                                            VCF_DIC2 [VCF_titles[6] + '|' + mut_prot_info + '|' + ref_peps + ':' + var_peps] = VCF_DIC2 [VCF_titles[6] + '|' + mut_prot_info + '|' + ref_peps + ':' + var_peps] + [[VCF_titles, ref_peps, var_peps, mut_info]]
                                        except:
                                            VCF_DIC2 [VCF_titles[6] + '|' + mut_prot_info + '|' + ref_peps + ':' + var_peps] = [[VCF_titles, ref_peps, var_peps, mut_info]]
                            else:
                                print (file_[:-4] + '\t' + '|'.join(VCF_[:10]) + '\t' + str(mut_site) + '\t' + str(CDS_list[0]) + '\t' + str(CDS_list[1]) + 'no protein sequences', file = syno_file)
                        except:
                            print (file_[:-4] + '\t' + '|'.join(VCF_[:10]) + '\t' + 'not protein sequences', file = syno_file)
                    #print>> log_file, ';'.join(r_pep_list) + '\t',
            print('total: ' + str(len(VCF_DIC)) + ', SAAVs: ' + str(n), file = log_file)
            pickle.dump(VCF_DIC2, dmp2_file)
            dmp2_file.close()
            log_file.close()                        
                    
            syno_file.close()      
#!/usr/bin/env python

#####################################################################
#   This file is the main part of SHAPELoop.
#   SHAPELoop is an RNA secondary structure prediction tool developed 
#   based on conserved SHAPE patterns of various loop motifs.
#   
#   Shapepattern can be used for multiple purposes:
#           1. Predict RNA secondary structures
#           2. Evaluate munally inspected or software predicted structures
#           3. Identify kissing loops in RNA
#           4. Detect structure switches of RNA
#   
#   Author: Jingyi Cao
#
#   Copyright (C) 2020  Jingyi Cao
#
#   SHAPELoop is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
#   Contact: caojingyi1258@outlook.com
#####################################################################

import pandas as pd
import numpy as np
from scipy.stats import wilcoxon 
from collections import Counter
from shape_rule_summary import *
from loop_motifs import *
from loop_score_for_candidates import *
from loop_classification import *
from select_candidate import *
from sys import argv
import os

#------
#get input
#------
ex_name = argv[1]
ex_seq = argv[2]
ex_shape = argv[3]
ens_num = argv[4]
bp = int(argv[5])
out_dir = argv[6]
subopt = int(argv[7])
tmp_path = os.path.split(os.path.realpath(__file__))[0]
shpatterns = tmp_path.replace("scripts","shape_pattern/") + "SHAPE_pattern.txt"
tool_fold = tmp_path.replace("scripts","rna_tools/fold") 
tool_pfs = tmp_path.replace("scripts","rna_tools/sample_pfs")

ex_ct = out_dir + '/guidance_structure/' + ex_name + '.guidance.ct'
ex_dot = out_dir + '/guidance_structure/' + ex_name + '.guidance.dot'
ex_pfs = out_dir + '/candidate_ensemble/' + ex_name + '.pfs'
ex_pfs_ct = out_dir + '/candidate_ensemble/' + ex_name + '.ensemble.ct'
ex_pfs_dot = out_dir + '//candidate_ensemble/' + ex_name + '.ensemble.dot' 


#------
#get guidance structure with SHAPE restraints
#------
print("---------------------------------------------------")
print("###Generating guidance structure with SHAPE restraints...")
print("---------------------------------------------------")
os.system("sh " + tool_fold + " "  + ex_seq + " " + ex_shape + " " + ex_ct + " " + ex_dot)


#------
#sample candidate structures with SHAPE renstraints
#------
print("---------------------------------------------------")
print("###Generating candidate ensemble with SHAPE restraints...")
print("---------------------------------------------------")
os.system("sh " + tool_pfs + " " + ex_seq + " " + ex_shape + " " + ex_pfs + " " + ex_pfs_ct + " " + ex_pfs_dot + " " + ens_num)
os.system("rm " + ex_pfs)


#------
#load seq, structures, and conserved SHAPE patterns
#-----
e_seq = open(ex_seq).readlines()[1].strip()
e_dot = open(ex_dot).readlines()[2].strip()
e_shape = pd.read_csv(ex_shape,header=None,sep='\t').iloc[:,1].tolist()
e_shape = [np.nan if x<=-999 else x for x in e_shape]
pfs_dot = pd.read_csv(ex_pfs_dot+".tmp",header=None).iloc[:,0].tolist()
#conserved SHAPE patterns
rule_dict = {}
for each_rule in open(shpatterns).readlines():
    tmp = each_rule.strip().split('\t')
    if tmp[0] in rule_dict.keys():
        rule_dict[tmp[0]].append((tmp[1],tmp[2]))
    else:
        rule_dict[tmp[0]]=[(tmp[1],tmp[2])]


#------
#get penalties for the guidance structure
#------
print("---------------------------------------------------")
print("###Get penalties for the guidance structure...\t\t\t\tdone.")
print("---------------------------------------------------")
score_mfe = get_score_for_RNA_candidate(e_dot,e_seq,e_shape,rule_dict,bp)
with open(out_dir+'/penalty/'+ex_name+'.guidance.penalty','w') as D:
    for k,v in score_mfe.items():
        D.write('\t'.join(k.split('|'))+'\t'+str(v)+'\n')


#------
#classify loops in guidance structure
#------
print("---------------------------------------------------")
print("###Classify loops in guidance structure...\t\t\t\tdone.")
print("---------------------------------------------------")
poor_loop,middle_loop,good_loop,poor_po,middle_po = get_good_middle_poor_loops(score_mfe)
print("#poor loop:")
for pl in poor_loop:
    print (pl)
print("\n#middle loop:")
for ml in middle_loop:
    print (ml)
print("\n#good loop:")
for gl in good_loop:
    print (gl)


#------
#get penalties for candidate structures
#------
print("---------------------------------------------------")
print("###Get penalties for candidate structures...\t\t\t\tdone.")
print("---------------------------------------------------")
can_score_tmp = []
for n,x in enumerate(pfs_dot):
    can_score_tmp.append(get_score_for_RNA_candidate(x,e_seq,e_shape,rule_dict))
with open(out_dir+'/penalty/'+ex_name+'.candidates.penalty','w') as D:
    for n,x in enumerate(can_score_tmp):
        for k,v in x.items():
            D.write(str(n+1)+'\t'+'\t'.join(k.split('|'))+'\t'+str(v)+'\n')


#------
#select candidate
#------
print("---------------------------------------------------")
print("###Select candidate...\t\t\t\tdone.")
print("---------------------------------------------------")

#evaluate candidates
if len(poor_loop)>0:
    score_for_can = second_socre_for_candidates(\
    can_score_tmp,poor_loop,poor_po,good_loop,middle_po)
    while len(poor_loop)>0 and len(score_for_can.keys())==0:
        poor_loop,poor_po,middle_loop,middle_po = \
        mod_middle_poor_loops(score_mfe,poor_loop,middle_loop)
        score_for_can = second_socre_for_candidates(\
        can_score_tmp,poor_loop,poor_po,good_loop,middle_po)
else:
    #if there is no poor loops
    score_for_can = poor_loop_0(can_score_tmp,good_loop,middle_po)


#select predicted structure based on evaluation scores
selected_id,selected_can = get_selected_candidate(score_for_can,e_dot,pfs_dot,\
np.sum([score_mfe[k] for k in poor_loop]))

print("#loop penalties for the selected candidate:")
if len(selected_id) ==1:
    print("guidance:")
    for k,v in score_mfe.items():
        print(k+"\t"+str(v))
    print(e_dot)
else:
    for sel_id in get_centroid(selected_id,pfs_dot):
        print(str(sel_id+1)+":")
        for k,v in can_score_tmp[sel_id].items():
            print(k+"\t"+str(v))
        print(pfs_dot[sel_id])


#------
#write output
#------
if len(selected_id) ==1:
    with open(out_dir+'/'+ex_name+'.SHAPELoop.dot','w') as D:
        D.write('> guidance '+ex_name+'\n'+e_seq+'\n'+e_dot+'\n')
else:
    with open(out_dir+'/'+ex_name+'.SHAPELoop.dot','w') as D:
        for sel_id in get_centroid(selected_id,pfs_dot):
            D.write('> '+str(sel_id+1)+" "+ex_name+'\n'+e_seq+'\n'+pfs_dot[sel_id]+'\n')
    if subopt == 1:
        with open(out_dir+'/'+ex_name+'.SHAPELoop.sup.dot','w') as D:
            for sel_id in selected_id:
                D.write('> '+str(sel_id+1)+" "+ex_name+'\n'+e_seq+'\n'+pfs_dot[sel_id]+'\n')



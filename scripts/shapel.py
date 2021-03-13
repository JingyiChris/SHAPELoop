#!/usr/bin/env python
#################################################
#    This file is the main part of SHAPELoop    #
#################################################

import pandas as pd
import numpy as np
from scipy.stats import wilcoxon 
from scipy import stats
from collections import Counter
from loop_motifs import *
from get_penalty import *
from evaluate_loops import *
from select_candidate import *
from sys import argv
import sys
import os

#------ get input ------
ex_name = argv[1]
ex_seq = argv[2]
ex_shape = argv[3]
ens_num = argv[4]
bp = int(argv[5])
out_dir = argv[6]
subopt = int(argv[7])
mode = argv[8]
mcfold = argv[9]

tmp_path = os.path.split(os.path.realpath(__file__))[0]
tool_fold = tmp_path.replace("scripts","rna_tools/fold") 
tool_pfs = tmp_path.replace("scripts","rna_tools/sample_pfs")

ex_ct = out_dir + '/guidance_structure/' + ex_name + '.guidance.ct'
ex_dot = out_dir + '/guidance_structure/' + ex_name + '.guidance.dot'
ex_pfs = out_dir + '/candidate_ensemble/' + ex_name + '.pfs'
ex_pfs_ct = out_dir + '/candidate_ensemble/' + ex_name + '.ensemble.ct'
ex_pfs_dot = out_dir + '/candidate_ensemble/' + ex_name + '.ensemble.dot' 


#------ get guidance structure with SHAPE restraints ------
print("-----------------------------------------------------------")
print(" Generating guidance structure with SHAPE restraints...")
print("----------------------------------------------------------")
os.system("sh " + tool_fold + " "  + ex_seq + " " + ex_shape + " " + ex_ct + " " + ex_dot)

#------ sample candidate structures with SHAPE renstraints ------
print("-----------------------------------------------------------")
print(" Generating candidate ensemble with SHAPE restraints...")
print("-----------------------------------------------------------")
os.system("sh " + tool_pfs + " " + ex_seq + " " + ex_shape + " " + ex_pfs + " " + ex_pfs_ct + " " + ex_pfs_dot + " " + ens_num)
os.system("rm " + ex_pfs)


#------ load seq, structures, and characteristic SHAPE patterns -----
e_dot = open(ex_dot).readlines()[2].strip()
e_seq = open(ex_seq).readlines()[1].strip()
e_shape = pd.read_csv(ex_shape,header=None,sep='\t').iloc[:,1].tolist()
e_shape = [np.nan if x<=-999 else x for x in e_shape]

if mcfold == 'false':
    pfs_dot = pd.read_csv(ex_pfs_dot+".tmp",header=None).iloc[:,0].tolist()
else:
    pfs_dot = pd.read_csv(mcfold,header=None).iloc[2:,0].tolist()

#------ get updated characteristic SHAPE patterns ------
if mode == 'default':
    shpatterns = open(tmp_path.replace("scripts","shape_pattern/") + "SHAPE_pattern.txt").readlines()
else:
    try:
        shpatterns = open(out_dir + "/SHAPE_pattern.txt").readlines()
    except IOError:
        print '\x1B[31mError!: No characteristic SHAPE patterns are identified\x1B[0m'
        raise SystemExit(1)
rule_dict = {}
for each_rule in shpatterns:
    tmp = each_rule.strip().split('\t')
    if tmp[0] in rule_dict.keys():
        rule_dict[tmp[0]].append((tmp[1],tmp[2]))
    else:
        rule_dict[tmp[0]]=[(tmp[1],tmp[2])]

#------ get penalties for the guidance structure ------
print("-----------------------------------------------------------")
print(" Get penalties for the guidance structure...")
print("-----------------------------------------------------------")
all_loops = get_dshape_for_loops(pfs_dot+[e_dot],e_seq,e_shape,rule_dict)
all_penalty_dict = get_score_for_loops(all_loops)

score = {}
can_hl = get_h(e_dot,e_seq,e_shape,'hairpin')
can_itl = get_it(e_dot,e_seq,e_shape,'internal')
can_bul = get_bu(e_dot,e_seq,e_shape,'bulge')
can_myall = merge_dicts([can_hl,can_itl,can_bul])
for k in can_myall.keys():
    if k in all_penalty_dict.keys():
        score[k] = all_penalty_dict[k]

with open(out_dir+'/penalty/'+ex_name+'.guidance.penalty','w') as D:
    for k,v in score.items():
        D.write('\t'.join(k.split('|'))+'\t'+str(v)+'\n')
print("Done.")


#------ classify loops in the guidance structure ------
print("-----------------------------------------------------------")
print(" Evaluate loops in guidance structure...")
print("-----------------------------------------------------------")
poor_loop,good_loop,middle_loop = get_good_poor_middle_loops(score)
print("Done.")
print("#poor loop:")
for pl in poor_loop:
    print (pl)
print("\n#good loop:")
for gl in good_loop:
    print (gl)
print("\n#middle loop:")
for ml in middle_loop:
    print (ml)


#------ get penalties for candidate structures ------
print("-----------------------------------------------------------")
print(" Get penalties for candidate structures...")
print("-----------------------------------------------------------")
can_score_tmp = []
for n,x in enumerate(pfs_dot):
    can_hl = get_h(x,e_seq,e_shape,'hairpin')
    can_itl = get_it(x,e_seq,e_shape,'internal')
    can_bul = get_bu(x,e_seq,e_shape,'bulge')
    can_myall = merge_dicts([can_hl,can_itl,can_bul])
    tmp_dic = {}
    for k in can_myall.keys():
        if k in all_penalty_dict.keys():
            tmp_dic[k] = all_penalty_dict[k]
    can_score_tmp.append(tmp_dic)

with open(out_dir+'/penalty/'+ex_name+'.candidates.penalty','w') as D:
    for n,x in enumerate(can_score_tmp):
        for k,v in x.items():
            D.write(str(n+1)+'\t'+'\t'.join(k.split('|'))+'\t'+str(v)+'\n')
print("Done.")

#------ select candidate ------
print("-----------------------------------------------------------")
print(" Select candidate...")
print("-----------------------------------------------------------")
if len(poor_loop)>0:
    score_for_can = second_socre_for_candidates(can_score_tmp,poor_loop,good_loop,middle_loop)
else:
    score_for_can = poor_loop_0(can_score_tmp,good_loop,middle_loop)
ctrl_score = np.sum(score.values())
selected_id = select_can(score_for_can,can_score_tmp,e_dot,ctrl_score)

if -1 in selected_id:
    print("guidance:")
    for k,v in score.items():
        print(k+"\t"+str(v))
    print(e_dot)
    with open(out_dir+'/'+ex_name+'.SHAPELoop.dot','w') as D:
        D.write('> guidance '+ex_name+'\n'+e_seq+'\n'+e_dot+'\n')
else:
    for x in get_centroid(selected_id,pfs_dot):
        print(str(x+1)+":")
        for k,v in can_score_tmp[x].items():
            print(k+"\t"+str(v))
        print(pfs_dot[x])
    with open(out_dir+'/'+ex_name+'.SHAPELoop.dot','w') as D:
        for x in get_centroid(selected_id,pfs_dot):
            D.write('> '+str(x+1)+" "+ex_name+'\n'+e_seq+'\n'+pfs_dot[x]+'\n')
    if subopt == 1:
        with open(out_dir+'/'+ex_name+'.SHAPELoop.sup.dot','w') as D:
            for x in get_uniq(selected_id,pfs_dot):
                D.write('> '+str(x+1)+" "+ex_name+'\n'+e_seq+'\n'+pfs_dot[x]+'\n')
print("done.")

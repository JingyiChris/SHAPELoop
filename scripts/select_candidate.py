#!/usr/bin/env python3
##############################################################
#    Select candidates according to their loop penalties     #
##############################################################
from evaluate_loops import *
import numpy as np

#----------------------------------------------------------------------------------
def get_centroid(strucid, pfs_dot):
    dis = []
    for x in strucid:
        tmp_sc = []
        for j in strucid:
            if x!=j:
                tmp_sc.append(np.sum(\
                [1 for n in range(len(pfs_dot[0])) if pfs_dot[x][n]!=pfs_dot[j][n]]))
        dis.append(np.sum(tmp_sc))
    choo_id = [strucid[n] for n,x in enumerate(dis) if x==np.min(dis)]
    return(get_uniq(choo_id,pfs_dot))

#----------------------------------------------------------------------------------
def get_uniq(strucid,pfs_dot):
    tmp_dot = [pfs_dot[strucid[0]]]
    f_id = [strucid[0]]
    for i in strucid[1:]:
        if pfs_dot[i] not in tmp_dot:
            f_id.append(i)
            tmp_dot.append(pfs_dot[i])
    return(f_id)

#----------------------------------------------------------------------------------
def get_selected_candidate(rna_score,poor_loop_score):
    f_selected = [-1]
    if len(rna_score.keys()) >0:
        can_score = [(k,np.sum(x)) for k,x in rna_score.items()]
        hscore = sorted([x[1] for x in can_score],reverse=True)[0]
        if hscore >= poor_loop_score:
            f_selected = [x[0] for x in can_score if x[1] == hscore]
    return(f_selected)

#----------------------------------------------------------------------------------
def select_can(rna_score,can_score_tmp,e_dot,poor_loop_score):

    if len(rna_score.keys())>0 and len(e_dot)<=250:
        sel_id = get_selected_candidate(rna_score,poor_loop_score)
    
    else:
        can_score_tmp = [np.sum(list(v.values())) for v in can_score_tmp]
        h_score = sorted(set(can_score_tmp),reverse=True)[:2]
        sel_id = [n for n,x in enumerate(can_score_tmp) if x in h_score]
    return(sel_id)

#################################################################################################################
#   Evaluate candidates according to their loop penalties and loop positions.                                   #
#                                                                                                               #
#   Desired RNA candidates are supposed to satisfy:                                                             #
#       1) There are no more loops beyond the regions of "poor", "middle" and "good" loops;                     #
#       2) all "good" loops are kept;                                                                           #
#       3) There are no loops with penalty higher than 1 in the "middle" loop regions;                          #
#       4) The sum of penalties for loops within the "poor" loop regions is the lowest among all candidates.    #
#                                                                                                               #
#   Input: loop_scores_for_candidates_dict, poor_loop, poor_loop_position,                                      #
#          good_loop, middle_loop_position, mfe_dot, pfs_dot                                                    #
#                                                                                                               #
#   Parameter: cutoff N used for defining poor/middle/good loops                                                #
#   Output: selected structure (dot_bracket)                                                                    #
#################################################################################################################
from loop_classification import *
import numpy as np

#=============================================================================================
def second_socre_for_candidates(can_scores,poor_loop,poor_po,good_loop,middle_po,cutoff=-1):
    '''
    Score RNA candidates according to their loop scores and loop positions.
    '''
    f_can_score = {}
    for n,can_score in enumerate(can_scores):
        poor_loops_count = np.sum([0 if x not in can_score.keys() else 1 for x in poor_loop ])
        good_loops_count = np.sum([0 if x in can_score.keys() else 1 for x in good_loop ])
        if poor_loops_count == 0 and good_loops_count ==0 :
            sum_score = []
            for k,v in can_score.items():
                if k not in good_loop:
                    can_po = get_positions(k,accurate=True)
                    if loop_contain_or_not(can_po,poor_po)==1:
                        sum_score.append(v)
                    elif loop_contain_or_not(can_po,middle_po)==1 :
                        if v=="NULL" or v>=cutoff:
                            pass
                        else:
                            sum_score.append(np.nan)
                    else:
                        sum_score.append(np.nan)
            f_can_score[n] = sum_score
            
    return(f_can_score)  

#=============================================================================================
def poor_loop_0(can_scores,good_loop,middle_po,cutoff=-1):
    '''If there is no poor loop'''
    f_can_score = {}
    for n,can_score in enumerate(can_scores):

        good_loops_count = np.sum([0 if x in can_score.keys() else 1 for x in good_loop ])
    
        if good_loops_count ==0:
            sum_score = []
            for k,v in can_score.items():
                if k not in good_loop:
                    can_po = get_positions(k,accurate=True)
                    if loop_contain_or_not(can_po,middle_po)==1:
                        if v=="NULL" or v>=cutoff:
                            pass
                        else:
                            sum_score.append(np.nan)
                    else:
                        sum_score.append(np.nan)
            f_can_score[n] = sum_score
            
    return(f_can_score) 
#============================================================================================
def get_centroid(strucid, pfs_dot):
    euclid_dis = []
    for x in strucid:
        tmp_sc = []
        for j in strucid:
            if x!=j:
                tmp_sc.append(np.sum(\
                [1 for n in range(len(pfs_dot[0])) if pfs_dot[x][n]!=pfs_dot[j][n]]))
        euclid_dis.append(np.sum(tmp_sc))
    choo_id = [strucid[n] for n,x in enumerate(euclid_dis) if x==np.min(euclid_dis)]
    #choo_struc = [pfs_dot[x] for x in choo_id]

    return(choo_id)

#===========================================================================================
def get_selected_candidate(rna_score,e_dot,pfs_dot,poor_loop_score):
    """
    Select the candidate with the highest score
    """
    f_selected = ['ctrl']
    if len(rna_score.keys()) >0:
        highest_score=-999
        candidates = [(k,np.sum(x)) for k,x in rna_score.items()
                      if np.nan not in x and 'NULL' not in x ]
        
        candidates_null = [(k,[i for i in x if i!= 'NULL'])
                           for k,x in rna_score.items()
                           if np.nan not in x and 'NULL' in x]
        if len(candidates)>0:
            hscore = sorted([x[1] for x in candidates],reverse=True)[0]
            if hscore >= poor_loop_score:
                f_selected = [x[0] for x in candidates if x[1] == hscore] 

        elif len([x for x in candidates_null if len(x[1])>0])>0:
            hscore = sorted([np.sum(x[1]) for x in candidates_null],reverse=True)[0]
            if hscore >= poor_loop_score:
                f_selected = [x[0] for x in candidates_null if np.sum(x[1]) == hscore]

    if 'ctrl' not in f_selected:
        se_can = [pfs_dot[x] for x in f_selected]
    else:
        se_can = e_dot

    return(f_selected,se_can)

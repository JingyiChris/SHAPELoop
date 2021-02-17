#!/usr/bin/env python3
######################################################################################
#   Evaluation for loops in guidance structures based on their penalties             #
######################################################################################
import numpy as np

#------------------------------------------------------------------------------------
def get_good_poor_loops(mfe_score):

    poor_loop = []
    good_loop = []
    
    for k,v in mfe_score.items():
        if v<-1:
            poor_loop.append(k)
        elif v>=-1:
            good_loop.append(k)

    return(poor_loop,good_loop)

#------------------------------------------------------------------------------------
def poor_loop_0(can_scores,good_loop):
    '''There is no poor loop'''
    f_can_score = {}
    for n,can_score in enumerate(can_scores):
        good_loops_count = np.sum([0 if x in can_score.keys() else 1 for x in good_loop ])
        if good_loops_count ==0 :
            sum_score = []
            for k,v in can_score.items():
                if k not in good_loop:
                    sum_score.append(v)
            f_can_score[n] = sum_score
    return(f_can_score)

#------------------------------------------------------------------------------------
def second_socre_for_candidates(can_scores,poor_loop,good_loop):
    f_can_score = {}
    for n,can_score in enumerate(can_scores):
        poor_loops_count = np.sum([0 if x not in can_score.keys() else 1 for x in poor_loop ])
        good_loops_count = np.sum([0 if x in can_score.keys() else 1 for x in good_loop ])
        if good_loops_count ==0 and poor_loops_count == 0:
            sum_score = []
            for k,v in can_score.items():
                if k not in good_loop:
                    sum_score.append(v)
            f_can_score[n] = sum_score
    return(f_can_score)

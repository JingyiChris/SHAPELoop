#!/usr/bin/env python3
###############################################
#    Get penalties for loops in candidates    #
###############################################
from scipy import stats
from loop_motifs import *

#-----------------------------------------------------------------------------------
def dshape_true_pdf(dshape_true):
    fit = (0.753, 0.593, 0.100, 0.260)
    return stats.norminvgauss.pdf(dshape_true,*fit)

def dshape_neg_pdf(dshape_neg):
    fit = (-0.283, 0.825, 0.003, 0.246)
    return stats.johnsonsu.pdf(dshape_neg,*fit)

def get_prob_for_dshape(dsh, prior_t=0.549,prior_f=0.451):
    prob_t = dshape_true_pdf(dsh)
    prob_f = dshape_neg_pdf(dsh)
    prob_dsh = 1-prob_t*prior_t/(prob_t*prior_t+prob_f*prior_f)
    return prob_dsh

#----------------------------------------------------------------------------------
def get_dshape_for_loops(dot_list,seq_tmp,shape_tmp,summary_dic):
    '''get loop positions for all candidates'''
    loop_pos = {}
    for d in dot_list:
        h_loop = get_h(d,seq_tmp,shape_tmp,'hairpin')
        it_loop = get_it(d,seq_tmp,shape_tmp,'internal')
        bu_loop = get_bu(d,seq_tmp,shape_tmp,'bulge')
        myall_tmp = merge_dicts([h_loop,it_loop,bu_loop])
        for k,v in myall_tmp.items():
            k2 = k.split("|")[0]+"|"+k.split("|")[1]+'|'
            if k2 in summary_dic.keys() and k not in loop_pos.keys():
                loop_pos[k] = []
                for v2 in summary_dic[k2]:
                    pos1_tmp = int(v2[0].split("_")[0])
                    pos2_tmp = int(v2[0].split("_")[1])
                    loop_pos[k].append(v[1][pos1_tmp]-v[1][pos2_tmp])
    return loop_pos

#----------------------------------------------------------------------------------
def get_score_for_loops(loop_pos):
    '''get panelties for loops in candidates'''
    loop_score = {}
    for k,v in loop_pos.items():
        loop_score[k] = 0
        for x in v:  
            p = get_prob_for_dshape(x)
            if p >0.5:
                loop_score[k] += np.log((1-p)/p)
        loop_score[k] = 0 if loop_score[k] > -0.1 else loop_score[k]
    return(loop_score)

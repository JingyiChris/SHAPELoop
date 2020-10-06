#####################################################################################
#   Get loop penalties for RNA candidates based on conserved SHAPE patterns         #
#   Input: dot_bracket structure, sequence, SHAPE reactivity, SHAPE_pattern_dict    #
#   Parameter: number of flanking base pairs (addbp=2)                              #
#   Ouput: penalty_dic, {loop_name|loop_length|loop_start_position: penalty, ...}   #
#####################################################################################

from loop_motifs import *
from shape_rule_summary import merge_dicts

#==================================================================================
def get_score_for_loops(candic,summary_dic):
    
    '''
	Get penalties for all loop motifs in a candidate
	'''
    
    s_each = {}
    for k,v in candic.items():
        k_tmp = k.split('|')[0]+'|'+k.split('|')[1]+'|'
        if k_tmp in summary_dic.keys():
            s_each[k] = 0
            for x in summary_dic[k_tmp]:
                po_1 = int(x[0].split("_")[0])
                po_2 = int(x[0].split("_")[1])
                if v[1][po_1] <= v[1][po_2] or v[1][po_1] < -500 or v[1][po_2] < -500:
                    s_each[k] +=-1        
        ##if there are no SHAPE rules summaried for this loop type: 'NULL'
        else:
            s_each[k] = "NULL"
    return s_each

#==================================================================================
def get_score_for_RNA_candidate(e_dot,e_seq,e_shape,summary_dics,addbp=2):
    ''' 
	Get loop score for all RNA candidates
	'''
    hl = get_h(e_dot,e_seq,e_shape,'hairpin',addbp)
    itl = get_it(e_dot,e_seq,e_shape,'internal',addbp)
    bul = get_bu(e_dot,e_seq,e_shape,'bulge',addbp)
    myall = merge_dicts([hl,itl,bul])
    score = get_score_for_loops(myall,summary_dics)
    return(score)

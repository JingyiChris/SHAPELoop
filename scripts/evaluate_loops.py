#!/usr/bin/env python3
############################################################################
#    Evaluate loops in the guidance structures based on their penalties    #
############################################################################
import numpy as np

#------------------------------------------------------------------------------------
def get_good_poor_middle_loops(mfe_score):

    poor_loop = []
    good_loop = []
    middle_loop = []
    
    for k,v in mfe_score.items():
        if v<-1:
            poor_loop.append(k)
        elif v==0:
            good_loop.append(k)
        else:
            middle_loop.append(k)
    return(poor_loop,good_loop,middle_loop)

#------------------------------------------------------------------------------------
def get_positions(name,accurate=True,h_extend=4,it_extend=4,bu_extend=4):

    '''get the positions of loops'''

    if "internal" in name:
        loop_5_len = int(name.split("|")[1].split("_")[0])
        loop_3_len = int(name.split("|")[1].split("_")[1])
        first_position = int(name.split("|")[2].split("%")[0])
        last_position = int(name.split("|")[2].split("%")[1])
        if accurate == True:
            node = [first_position, first_position+loop_5_len, \
            last_position-loop_3_len, last_position]
        else:
            node = [first_position-it_extend, first_position+loop_5_len+it_extend, \
            last_position-loop_3_len-it_extend, last_position+it_extend]

    elif "hairpin" in name:
        start = int(name.split("|")[2])
        length = int(name.split("|")[1])
        if accurate == True:
            node = [start,start+length]
        else:
            node = [start-h_extend,start+length+h_extend]


    elif "bulge" in name:
        first_position = int(name.split("|")[2].split("%")[0])
        last_position = int(name.split("|")[2].split("%")[1])
        length = int(name.split("|")[1])
        if accurate == True:
            node = [first_position, first_position+length, \
            last_position-1, last_position]
        else:
            node = [first_position-bu_extend, first_position+length+bu_extend, \
            last_position-bu_extend, last_position+bu_extend]
    return node
    
#------------------------------------------------------------------------------------
def extended_positions(position_list):

    '''get extended positions for loops'''

    son_set = []
    for x in position_list:
        #hairpin
        if len(x)==2:
            son_set += list(range(x[0],x[1]))
        #internal/bulge
        elif len(x)==4:
            son_set += list(range(x[0],x[1]))+list(range(x[2],x[3]))
    return(list(set(son_set)))

#------------------------------------------------------------------------------------
def loop_contain_or_not(loop1_positions,loop2_positions):

    ''' determine if the postion of loop1 is in the position of loop2'''

    original_s = 0
    son_set = []
    if len(loop1_positions)==2:
        son_set += list(range(loop1_positions[0],loop1_positions[1]))
    elif len(loop1_positions)==4:
        son_set += list(range(loop1_positions[0],loop1_positions[1]))+\
        list(range(loop1_positions[2],loop1_positions[3]))
    ##issubset?
    if set(son_set).issubset(set(loop2_positions)):
        original_s = 1
    return(original_s)

#------------------------------------------------------------------------------------
def poor_loop_0(can_scores,good_loop,middle_loop):

    '''if there is no poor loop'''

    f_can_score = {}
    for n,can_score in enumerate(can_scores):
        good_loops_count = np.sum([0 if x in can_score.keys() else 1 for x in good_loop ])
        pos_m = []
        if len(middle_loop) >0:
            pos_m = extended_positions([get_positions(x,accurate=False) for x in middle_loop])
        if good_loops_count ==0 :
            sum_score = []
            for k,v in can_score.items():
                if k not in good_loop:
                    pos_tmp = extended_positions(get_positions(k))
                    if v <-1 and loop_contain_or_not(pos_tmp,pos_m)==1:
                        pass
                    else:
                        sum_score.append(v)
            f_can_score[n] = sum_score
    return(f_can_score)

#------------------------------------------------------------------------------------
def second_socre_for_candidates(can_scores,poor_loop,good_loop,middle_loop):
    f_can_score = {}
    for n,can_score in enumerate(can_scores):
        poor_loops_count = np.sum([0 if x not in can_score.keys() else 1 for x in poor_loop ])
        good_loops_count = np.sum([0 if x in can_score.keys() else 1 for x in good_loop ])
        pos_m = []
        if len(middle_loop) >0:
            pos_m = extended_positions([get_positions(x,accurate=False) for x in middle_loop]) 
        if good_loops_count ==0 and poor_loops_count == 0:
            sum_score = []
            for k,v in can_score.items():
                if k not in good_loop:
                    pos_tmp = extended_positions([get_positions(k)])
                    if v <-1 and loop_contain_or_not(pos_tmp,pos_m)==1:
                        pass
                    else:
                        sum_score.append(v)
            f_can_score[n] = sum_score
    return(f_can_score)

######################################################################################
#   Classification of loop motifs in the guidance structure based on their penalties #
#	                                                                                 #
#	Poor loop: loops with penalty > 1                                                #
#   Middle loops: loops with "NULL" or 1 penalty                                     #
#   Good loops: Loops with penalty < 1                                               #
#                                                                                    #
######################################################################################
import numpy as np

#====================================================================================
def get_positions(name,accurate=True,h_extend=4,it_extend=4,bu_extend=4):
    
    '''get the positions of each loop in its RNA'''

    if "internal" in name:
        loop_5_len = int(name.split("|")[1].split("_")[0])
        loop_3_len = int(name.split("|")[1].split("_")[1])
        first_position = int(name.split("|")[2].split("%")[0])
        last_position = int(name.split("|")[2].split("%")[1])
        if accurate == True:
            node = [first_position, first_position+loop_5_len,
            last_position-loop_3_len, last_position]
        else:
            node = [first_position-it_extend, first_position+loop_5_len+it_extend,
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
            node = [first_position, first_position+length,
            last_position-1, last_position]
        else:
            node = [first_position-bu_extend, first_position+length+bu_extend,
            last_position-bu_extend, last_position+bu_extend]
    return node
#=====================================================================================
def combine_positions(position_list):
    '''get combined positions for loops'''
    son_set = []
    for x in position_list:
        #hairpin
        if len(x)==2:
            son_set += list(range(x[0],x[1]))
        #internal/bulge
        elif len(x)==4:
            son_set += list(range(x[0],x[1]))+list(range(x[2],x[3]))
    return(list(set(son_set)))

#=====================================================================================
def loop_contain_or_not(loop1_positions,loop2_positions):
    '''
    Determine if the postion of loop1 is compeletly 
    contained in the position of loop2
    '''
    original_s = 0
    son_set = []
    if len(loop1_positions)==2:
        son_set += list(range(loop1_positions[0],loop1_positions[1]))
    elif len(loop1_positions)==4:
        son_set += list(range(loop1_positions[0],loop1_positions[1]))+ \
        list(range(loop1_positions[2],loop1_positions[3]))
    ##issubset?
    if set(son_set).issubset(set(loop2_positions)):
        original_s = 1
    return(original_s)

#=====================================================================================
def get_good_middle_poor_loops(mfe_score,cutoff=-1):
    '''
    loops with penalty < -1 are defiend as poor_loops
    loops with "NULL" and -1 penalty are defined as middle loops
    loops with penalty = 0 are defined as good loops
    good loops should be retained
    poor loops should be removed

    '''
    poor_loop = []
    middle_loop = []
    good_loop = []
    poor_po = []
    middle_po = []
    
    if len([x for x in mfe_score.values() if x!="NULL" and x <= cutoff])>0:
        for k,v in mfe_score.items():
            if v!="NULL" and v < cutoff:
                poor_loop.append(k)
                poor_po.append(get_positions(k,accurate=False))
    
    
        #get the positions of poor loops
        poor_po = combine_positions(poor_po)
    
        for k,v in mfe_score.items():
            if k not in poor_loop:
                pos = get_positions(k,accurate=False)
                
                if v=="NULL" or v == cutoff:
                    middle_loop.append(k)
                    middle_po.append(pos)
                elif v > cutoff:
                    good_loop.append(k)
            
        #get the positions of middle loops
        middle_po = combine_positions(middle_po)
    
    return(poor_loop,middle_loop,good_loop,poor_po,middle_po)

#=====================================================================================
def mod_middle_poor_loops(mfe_score,poor_loop,middle_loop):
    ''' 
    If there are too many poor loops in MFE structure, and the variety of pfs sampled candidates,
    however, is low, then the cutoff for defining poor loops needs to be stricter,
    so poor loops with the highest penalty score will be moved to the group of middle loops
    '''
    highest_score = sorted([mfe_score[x] for x in poor_loop],reverse=True)[0]
   
    poor_loop_update = [x for x in poor_loop if mfe_score[x]!=highest_score]
    poor_po = [get_positions(x,accurate=False) for x in poor_loop_update]
    poor_po = combine_positions(poor_po)
    
    middle_loop_update = middle_loop+[x for x in poor_loop if mfe_score[x]==highest_score]
    middle_po = [get_positions(x,accurate=True) for x in middle_loop_update]
    middle_po = combine_positions(middle_po)

    return(poor_loop_update,poor_po,middle_loop_update,middle_po)


#######################################################################
#   Extract hairpin, internal and bulge loops from given structures   # 
#                                                                     #
#   Input: dot_bracket structure, sequence, shape reactivity          #
#                                                                     #
#   Parameters: number of flanking base pairs (addbp=2)               #
#                                                                     #
#   Output: loop dictionary:                                          #
#           {loop_name|loop_length|loop_position:                     #
#           (loop_seq,loop_shape_list)                                #
#           ... }                                                     #
#######################################################################
import pandas as pd
import numpy as np

#===============================================================================
def get_pairs(structure):
    '''
    get base pair relationships from structures: (res1,res2), unpaired: -1
    '''
    pairs = []
    S = []
    for i, s in enumerate(structure):
        if s == '(':
            S.append(i)
        elif s == ')':
            if len(S) == 0:
                raise ValueError('extra ) found in structure at position {}'.format(i))
            j = S.pop()
            pairs.append((j, i))
            pairs.append((i,j))
        elif s == '.':
            pairs.append((i,-1))
    return sorted(pairs)
#===============================================================================
def get_h(dot,seq,shape,name,addbp=2):
    '''
    get hairpin loops from structures
    '''
    h_loop = {}
    pairs = get_pairs(dot)
    for num,se in enumerate(seq):
        if num>=addbp and num<len(seq) and pairs[num][1]==-1:
            if np.sum([1 for x in range(1,1+addbp) if pairs[num-x][1]>-1])==addbp:
                loop_num = 1
                
                while (num+loop_num)<len(seq) and pairs[num+loop_num][1]==-1:
                    loop_num+=1
                if np.sum([1 for x in range(1,1+addbp) if pairs[num-x][1]==num+loop_num-1+x])==addbp:

                    h_loop[name+"|"+str(loop_num)+"|"+str(num)] = ( 
                    ''.join([seq[x] for x in range(num,num+loop_num)]),
                    [shape[x] for x in range(num,num+loop_num)])
    return h_loop
#===============================================================================
def get_bu(dot,seq,shape,name,addbp=2):
    '''
    get bulge loops from structures
    '''
    bu_loop = {}
    pairs = get_pairs(dot)
    
    for num,se in enumerate(seq):
        if num>=addbp and num<len(seq)-addbp and pairs[num][1]==-1:
            if np.sum([1 for x in range(1,1+addbp) if pairs[num-x][1]>-1])==addbp:
            
                loop_num=1
                while (num+loop_num)<len(seq)-2 and pairs[num+loop_num][1]==-1:
                    loop_num+=1
                if np.sum([1 for x in range(1,1+addbp) if pairs[num-x][1] - pairs[num+loop_num-1+x][1]==x*2-1])==addbp:
           
                    if loop_num==1:
                        bu_loop[name+"|"+str(loop_num)+"|"+\
                        str(num)+"%"+str(pairs[num-1][1])] = ( 
                        ''.join([seq[x] for x in range(num-1,num+loop_num+1)]),
                        [shape[x] for x in range(num-1,num+loop_num+1)])
                    else: 
                        bu_loop[name+"|"+str(loop_num)+"|"+\
                        str(num)+"%"+str(pairs[num-1][1])] = ( 
                        ''.join([seq[x] for x in range(num,num+loop_num)]),
                        [shape[x] for x in range(num,num+loop_num)])
                
    return bu_loop
#===============================================================================
def get_it(dot,seq,shape,name,addbp=2):
    '''
    get internal loops from structures
    '''
	
    it_loop = {}
    pairs = get_pairs(dot)
   
    for num,se in enumerate(seq):
        if num>=addbp and num<len(seq)-addbp and pairs[num][1]==-1:
            if np.sum([1 for x in range(1,1+addbp) if pairs[num-x][1]>num-x])==addbp:
                loop_num1 = 1 
                while (num+loop_num1)<len(seq)-addbp and pairs[num+loop_num1][1]==-1:
                    loop_num1+=1
                
                ##keybase stands for the last residue at the 3'-strand of the internal loop
                keybase = pairs[num-1][1]-1
                if pairs[keybase][1]==-1 and keybase!=num and keybase!=(num+loop_num1-1):
                    loop_num2=1
                    while pairs[keybase-loop_num2][1]==-1:
                        loop_num2+=1
                   
                    if pairs[num+loop_num1][1]==keybase-loop_num2 \
                    and pairs[num+loop_num1+1][1]==keybase-loop_num2-1:                
                        it_loop[name+"|"+str(loop_num1)+"_"+str(loop_num2)+"|"+\
                        str(num)+"%"+str(keybase)] = ( 
                        ''.join([seq[x] for x in list(range(num,num+loop_num1))+ 
                        list(range(keybase-loop_num2+1,keybase+1))]),
                        [shape[x] for x in list(range(num,num+loop_num1))+
                        list(range(keybase-loop_num2+1,keybase+1))])
                    
    return it_loop



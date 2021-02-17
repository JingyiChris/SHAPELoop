#!/usr/bin/env python3

#######################################################################
#   Extract hairpin, internal and bulge loops from given structures   # 
#                                                                     #
#   Input: dot_bracket structure, sequence, shape reactivity          #
#                                                                     #
#   Parameters: number of flanking base pairs (addbp=2)               #
#                                                                     #
#######################################################################
import pandas as pd
import numpy as np
from scipy.stats import wilcoxon
from collections import Counter

# --------------------------------------------------------------------------------
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

# --------------------------------------------------------------------------------
def get_h(dot,seq,shape,name,addbp=2):
    #extract hairpin loops from RNA
    h_loop = {}
    pairs = get_pairs(dot)
    
    for num,se in enumerate(seq):
        if num>=addbp and num<len(seq) and pairs[num][1]==-1:

            if np.sum([1 for x in range(1,1+addbp) if pairs[num-x][1]>-1])==addbp:
 
                loop_num = 1
                while (num+loop_num)<len(seq) and pairs[num+loop_num][1]==-1:
                    loop_num+=1
                if np.sum([1 for x in range(1,1+addbp) if pairs[num-x][1]==num+loop_num-1+x])==addbp:

                    h_loop[name+"|"+str(loop_num)+"|"+str(num)] = ( \
                    ''.join([seq[x] for x in range(num-addbp,num+loop_num+addbp)]),\
                    [shape[x] if shape[x]>-500 else np.nan for x in range(num,num+loop_num)])
    return h_loop

# --------------------------------------------------------------------------------
def get_it(dot,seq,shape,name,addbp=2):
    ##extract internal loops from RNA, from 5' to 3' 
    it_loop = {}
    pairs = get_pairs(dot)
   
    for num,se in enumerate(seq):
        if num>=addbp and num<len(seq)-addbp and pairs[num][1]==-1:
            if np.sum([1 for x in range(1,1+addbp) if pairs[num-x][1]>num-x])==addbp:
                loop_num1 = 1 
                while (num+loop_num1)<len(seq)-addbp and pairs[num+loop_num1][1]==-1:
                    loop_num1+=1
           
                keybase = pairs[num-1][1]-1
                if pairs[keybase][1]==-1 and keybase!=num and keybase!=(num+loop_num1-1):
                    loop_num2=1
                    while pairs[keybase-loop_num2][1]==-1:
                        loop_num2+=1
                   
                    ###make sure internal loop
                    if pairs[num+loop_num1][1]==keybase-loop_num2 \
                    and pairs[num+loop_num1+1][1]==keybase-loop_num2-1:                
                        it_loop[name+"|"+str(loop_num1)+"_"+str(loop_num2)+"|"+\
                        str(num)+"%"+str(keybase)] = ( \
                        ''.join([seq[x] for x in list(range(num,num+loop_num1))+ \
                        list(range(keybase-loop_num2+1,keybase+1))]),\
                        [shape[x] if shape[x]>-500 else np.nan for x in list(range(num,num+loop_num1))+\
                        list(range(keybase-loop_num2+1,keybase+1))])
                    
    return it_loop

# --------------------------------------------------------------------------------
def get_bu(dot,seq,shape,name,addbp=2):
    #extract bulge loops from RNA 
    bu_loop = {}
    pairs = get_pairs(dot)
    
    
    for num,se in enumerate(seq):
        if num>=addbp and num<len(seq)-addbp and pairs[num][1]==-1:
            if np.sum([1 for x in range(1,1+addbp) if pairs[num-x][1]>-1])==addbp:
            
                loop_num=1
                while (num+loop_num)<len(seq)-2 and pairs[num+loop_num][1]==-1:
                    loop_num+=1
                if np.sum([1 for x in range(1,1+addbp) if pairs[num-x][1] - pairs[num+loop_num-1+x][1]==x*2-1])==addbp:
           
                    if loop_num!=1: 
                        bu_loop[name+"|"+str(loop_num)+"|"+\
                        str(num)+"%"+str(pairs[num-1][1])] = ( \
                        ''.join([seq[x] for x in range(num,num+loop_num)]),\
                        [shape[x] if shape[x]>-500 else np.nan for x in range(num,num+loop_num)])
                
    return bu_loop
# --------------------------------------------------------------------------------
def merge_dicts(diclist):
    result = {}
    for dictionary in diclist:
        result.update(dictionary)
    return result

# --------------------------------------------------------------------------------
def mask_dic(dic):
    for key in list(dic):
        if not dic.get(key):
            del dic[key]
    return dic
# --------------------------------------------------------------------------------
def get_statistical_summary(loopdic,looptype,pvalue=0.05,topN=2):
    f_dic = {}
    f_dic_tmp = {}
    f_dic_topN = {}

    if looptype != 'internal':
        for k,v in loopdic.items():
            l_name = looptype+"|"+k.split("|")[1]+"|"
            if l_name in f_dic.keys():
                f_dic[l_name].append(v[1])
            else:
                f_dic[l_name]=[v[1]]
    
    elif looptype == 'internal':
        for k,v in loopdic.items():
            l_name = looptype+"|"+k.split("|")[1]+"|"
            strand5 = k.split("|")[1].split('_')[0]
            strand3 = k.split("|")[1].split('_')[1]
            
            if l_name in f_dic.keys():
                f_dic[l_name].append(v[1])
            else:
                f_dic[l_name]=[v[1]]
            
            if strand5 != strand3:
                l_name_rev = looptype+"|"+strand3+"_"+strand5+"|"
                if l_name_rev in f_dic.keys():
                    f_dic[l_name_rev].append(v[1][-int(strand3):]+v[1][:int(strand5)])
                else:
                    f_dic[l_name_rev]=[v[1][-int(strand3):]+v[1][:int(strand5)]]
    
    for k,v in f_dic.items():
        comp_array = np.array(v).T
        f_dic_tmp[k] = []
        for n1 in range(comp_array.shape[0]):
            for n2 in range(comp_array.shape[0]):
                
                comp_tmp = np.delete(comp_array[[n1,n2]],\
                    np.where(np.isnan(comp_array[[n1,n2]]))[1],axis=1)
                if n2>n1 and (comp_tmp[0]==comp_tmp[1]).all() == False:
                    p_g = wilcoxon(comp_tmp[0], comp_tmp[1],alternative='greater')[1]
                    p_l = wilcoxon(comp_tmp[0], comp_tmp[1],alternative='less')[1]
                    if p_g < pvalue:
                        sym = str(n1)+"_"+str(n2)
                        f_dic_tmp[k].append((sym,p_g))
                    elif p_l <pvalue:
                        sym = str(n2)+"_"+str(n1)
                        f_dic_tmp[k].append((sym,p_l))
        
        if f_dic_tmp[k]:
            rank_tmp = sorted(Counter([x[1] for x in f_dic_tmp[k]]).keys())[:topN]
            f_dic_topN[k] = [x for x in f_dic_tmp[k] if x[1] in rank_tmp]

    return(mask_dic(f_dic_topN))

##################################################################################################
#   Summarize conserved SHAPE patterns for loop motifs                                           #
#                                                                                                #
#   example:                                                                                     #
#        hairpin|4|:  1>0 (SHAPE reactivity of residue 2 is higher than that of residue 1),      #
#                     2>0 (SHAPE reactivity of residue 3 is higher than that of residue 1)       #
#                                                                                                #
#   Input: loop_dict, loop_type                                                                  #
#                                                                                                #
#   Parameters:                                                                                  #
#       merge_internal: merge asymmetric internal loops, such as 5'-xxx|3'-xx and 5'-xx|3'-xxx   #
#       pvalue: P-value cutoff for wilcoxon test (default is 0.05),                              #
#       topN: keep SHAPE patterns with the lowest N p-value (default is 2)                       #
#                                                                                                #
#   Output: rule_dict, {loop_type|loop_length|: [(shape_rule, P-value), ...]...}                 #
##################################################################################################

import pandas as pd
import numpy as np
from scipy.stats import wilcoxon 
from collections import Counter

#====================================================================================
def merge_dicts(diclist):
    result = {}
    for dictionary in diclist:
        result.update(dictionary)
    return result
#====================================================================================
def mask_dic(dic):
    for key in list(dic):
        if not dic.get(key):
            del dic[key]
    return dic
#====================================================================================
def get_statistical_summary(loopdic,looptype,pvalue=0.05,topN=2,merge_internal=True):
    f_dic = {}
    f_dic_top2 = {}
    pelist = list(set([x.split("|")[1] for x in loopdic.keys()]))
    
    for pe_tmp in pelist:
        pe = '|'+pe_tmp+'|'
        f_dic[looptype+pe] = []
        compare_list = []
        se = [x  for x in loopdic.keys() if pe in x ]
        
        ##Get array for wilcoxon test
        if looptype != 'internal' or merge_internal == False:
            for x in range(len(loopdic[se[0]][1])):
                compare_list.append([loopdic[p][1][x] for p in se])
        
        else:
           
            ##Merge asymmetric internal loops: such as 5'-GAA|3'-CG and 5'-GC|3'-AAG
            tmp1 = pe_tmp.split("_")[0]
            tmp2 = pe_tmp.split("_")[1]
            update_dic = {}
            for k,v in loopdic.items():
                if pe_tmp in k:
                    update_dic[k] = v 
                if "|"+tmp2+"_"+tmp1+"|" in k and tmp1!=tmp2:
                    update_dic[k] = (v[0],v[1][-int(tmp1):]+v[1][:int(tmp2)])
            
            for x in range(int(tmp1)+int(tmp2)):
                compare_list.append([v[1][x] for k,v in update_dic.items()])
            
        
        #Compare SHAPE reactivities between each two residue with wilcoxon test
        compare_list = np.array(compare_list)
        for n1,j in enumerate(compare_list):
            for n2,k in enumerate(compare_list):
                if n1<n2:
                    test_2 = compare_list[[n1,n2]]
                    test_2_tmp = np.delete(test_2,np.where(test_2 < -500)[1],axis=1)
                    
                    if (test_2_tmp[0]==test_2_tmp[1]).all() == False :
                        p_value_g = wilcoxon(test_2_tmp[0],test_2_tmp[1],alternative='greater')[1]
                        p_value_l = wilcoxon(test_2_tmp[0],test_2_tmp[1],alternative='less')[1]
                        
                        if p_value_g <pvalue:
                            sym = str(n1)+"_"+str(n2)
                            f_dic[looptype+pe].append((sym,p_value_g))
                        if p_value_l <pvalue:
                            sym = str(n2)+"_"+str(n1)
                            f_dic[looptype+pe].append((sym,p_value_l))
            
        #filter SHAPE rules for each loop type: keep SHAPE rules with the two smallest P-value
        if f_dic[looptype+pe]:
            f_dic_top2[looptype+pe] = []
            p_list = [xp[1] for xp in f_dic[looptype+pe]]
            rank_tmp = sorted(Counter(p_list).keys())[:topN]
            for n_p,xp in enumerate(p_list):
                if xp in rank_tmp:
                    f_dic_top2[looptype+pe].append((f_dic[looptype+pe][n_p][0],xp))
    return(mask_dic(f_dic_top2))

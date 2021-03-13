#!/usr/bin/env python3
#################################################################
#    Summarize characteristic SHAPE patterns for loop motifs    #
#################################################################
import pandas as pd
import numpy as np
from scipy.stats import wilcoxon
from collections import Counter
import warnings
warnings.filterwarnings('ignore')
from loop_motifs import *
from sys import argv
import os

#--------------------------- get input ----------------------------------
input_data = open(argv[1]).readlines()
mode_tmp = argv[2]
out_dir = argv[3]
tmp_path = os.path.split(os.path.realpath(__file__))[0]

if mode_tmp == 'combined':
    h_ori = pd.read_csv(tmp_path.replace("scripts","shape_pattern/") +\
        "hairpins.txt",header=None,sep='\t')
    it_ori = pd.read_csv(tmp_path.replace("scripts","shape_pattern/") + \
        "internals.txt",header=None,sep='\t')
    bu_ori = pd.read_csv(tmp_path.replace("scripts","shape_pattern/") + \
        "bulges.txt",header=None,sep='\t')

#--------------------------- get SHAPE patterns -------------------------
h_loop = []
it_loop = []
bu_loop = []
for tmp in range(int(len(input_data)/4)):
    r_name = input_data[tmp*4].strip()
    true_seq = input_data[tmp*4+1].strip()
    true_dot = input_data[tmp*4+2].strip()
    true_shape = [float(x) for x in input_data[tmp*4+3].strip().split(',')]
    
    h_loop.append(get_h(true_dot,true_seq,true_shape,r_name))
    it_loop.append(get_it(true_dot,true_seq,true_shape,r_name))
    bu_loop.append(get_bu(true_dot,true_seq,true_shape,r_name))
    
if mode_tmp == 'combined':
    h_loop.append(dict(zip(h_ori.iloc[:,0],\
        [(0,[float(x1) for x1 in x.split(",")]) for x in h_ori.iloc[:,1]])))
    it_loop.append(dict(zip(it_ori.iloc[:,0],\
        [(0,[float(x1) for x1 in x.split(",")]) for x in it_ori.iloc[:,1]])))
    bu_loop.append(dict(zip(bu_ori.iloc[:,0],\
        [(0,[float(x1) for x1 in x.split(",")]) for x in bu_ori.iloc[:,1]])))

hairpin_sum = get_statistical_summary(merge_dicts(h_loop),'hairpin')
internal_sum = get_statistical_summary(merge_dicts(it_loop),'internal')
bulge_sum = get_statistical_summary(merge_dicts(bu_loop),'bulge')
rule_dic = merge_dicts([hairpin_sum,internal_sum,bulge_sum])
if rule_dic:
    print('--- SHAPE pattern identification done! ---')
    with open(out_dir+'/SHAPE_pattern.txt','w') as D:
        for k,v in rule_dic.items():
            for x in v:
                D.write(k+'\t'+x[0]+'\t'+str(x[1])+'\n')
else:
    print('--- No SHAPE patterns identified! ---')

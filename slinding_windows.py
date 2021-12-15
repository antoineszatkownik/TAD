# -*- coding: utf-8 -*-
"""
Created on Wed Dec 15 16:57:16 2021

@author: maxim
"""

import numpy as np
import random
import pandas as pd 

with open("dataforstudent\\Arrowhead\\GSE63525_GM12878_primary+replicate_Arrowhead_domainlist.txt", "r") as f:
    mat = np.array(pd.read_csv("Chrom_1.csv", sep = " "))
    chrom = 1
    for l in f.readlines()[1:]:
        tmp = l.split("\t")
        r = random.randint(10,40) #radom offset from arrowhead windows
        coord1 = int(np.ceil(((int(tmp[1])//25000)*25000)/25000))-r #remove r for label 1
        coord2 = int(np.ceil(((int(tmp[2])//25000)*25000)/25000))-r #remove r for label 1
        if chrom != int(tmp[0]):
            chrom = int(tmp[0])
            mat = np.array(pd.read_csv("Chrom_{}.csv".format(chrom), sep = " "))
        print(coord1,coord2)
        try:
            np.savetxt("data\\0\\chrom_{}_{}.csv".format(chrom,coord1), mat[coord1-17:coord1+16,coord1-17:coord1+16])
            np.savetxt("data\\0\\chrom_{}_{}.csv".format(chrom,coord2), mat[coord2-17:coord2+16,coord2-17:coord2+16])
        except IndexError:
            continue
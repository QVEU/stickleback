#!/usr/bin/env python
'''
stickleback.0.0.py

Patrick T. Dolan
Unit Chief, Quantitative Virology and Evolution Unit
6/1/22
USAGE: stickleback.0.0.py pathto/input.sam query templateFasta [] []   ( [arg] = optional argument)
Designed to map inserted sequences in a template sequence (plasmid or transcript).
'''

##### Imports #####
import pandas as pd
import Levenshtein
from plotnine import *
import re
import numpy as np
import sys
from multiprocessing import Pool

##### Functions #####
def Initialize(args):
    '''
    Function: Initialize()
        Reads in and filters data and initializes a bunch of parameters for the mapping.

    Arguments:
        "seq_query" is a *TUPLE* of the template sequence and the query. Need to be sent as a package for multi-threading with Pool.
    Note: replaces position if alternative minimum sites are found. Producing a 3' bias? Need to consider improvements here.
    '''
    #Input SAM
    print("\n\n==============================================")
    print("============     stickleback     =============")
    print("==============================================\n\n")

    #Query String
    query=args[2]
    print("Query Length: {}".format(len(query)))

    #Template File
    template=args[3]
    with open(template,'r') as fasta:
        fa=fasta.readlines()
    templateSeq="".join([i.strip() for i in fa][1:])
    print("Template Length: {}".format(len(templateSeq)))

    Dtemp=distCompute((templateSeq,query))[0] #compute the distance to the template (without insert) to determine sensitivity.
    print("Template-Query Distance: {}".format(Dtemp))

    #Choose Read Size
    if len(args)>2:
        minL=int(args[4]) # read length minimum
        maxL=int(args[5]) # read length maximum
    else:
        minL=len(query)
        maxL=len(templateSeq)

    print("Mapping Reads of Size {} to {} with a cutoff of {}".format(minL, maxL, int(Dtemp)-1))

    #input SAM
    inputData=sys.argv[1]
    print("Loading SAM: {}".format(inputData))
    pdSam=loadSAM(inputData,minL,maxL)

    print("===============================================\n")

    return(query,templateSeq, pdSam, Dtemp)

def loadSAM(inputData, minL, maxL):
    '''
    Function: loadSAM()
        Reads in SAM file.
        Since this is the slow step, I broke this out to work on.

    Arguments:
        "inputData" is a string of the path to the input SAM file.
    Note: quite slow, how can we speed up? pd.read seems to have issues with the headers on the file.
    '''
    with open(inputData,"r") as IF:
        samInput=[lin.split("\t")[0:11] for lin in IF if lin.startswith("@") is False]#Filter out headers with @
    
    names=["read","FLAG","template","pos","mapq","cigar","Rnext","Pnext","Tlen","seq","Qscore"]#SAM columns
    size=len(samInput)
    print(size)
    pdSam=pd.DataFrame(samInput,columns=names)
    pdSam['length']=pdSam.seq.apply(lambda s: len(s))
    pdSam=pdSam[(pdSam.length>minL)&(pdSam.length<maxL)&(pdSam.template!="*")]
    return(pdSam)

def distCompute(seq_query):
    '''
    Function: distCompute
        Computes the Levenshtein distance of the query sequence to each window of the target sequence, returning the position, minimum distance, context, and
    Arguments:
        "seq_query" is a *TUPLE* of the template sequence and the query. Need to be sent as a package for multi-threading with Pool.
    Note: replaces position if alternative minimum sites are found. Producing a 3' bias? Need to consider improvements here.
    '''
    #Unpack Values
    seq=seq_query[0]
    query=seq_query[1]

    #Prepare Query
    quL=len(query)
    minD=quL #set minimum distance to length of query

    #Initial Values
    minPos=0
    matchseq="X"
    contextSeq=""
    context=""
    matchlist=[]
    for i in range(len(seq)-quL):#for each position in
        matchlist=[]
        D=Levenshtein.distance(seq[i:i+quL].upper(),query.upper())
        match=seq[i:i+quL]
        context=seq[max(0,i-quL):i]+"|"+seq[i+quL:min(len(seq),i+(quL)*2)]
        '''print statements
        print(str(i)+"\nD:",str(D))
        print("minD:"+str(minD))
        print("minPos:"+str(minPos))
        '''
        minD=min(D,minD)
        if minD==D:
            minPos=i
            matchseq=match
            contextSeq=context
            matchlist=matchlist.append(minPos)
    return(minD,minPos,matchseq,contextSeq,matchlist)

def poolDists(query,pdSam,nthreads=16):
    with Pool(nthreads) as p:
        print("1. Computing minimum distance hit position for {} reads.".format(len(pdSam)))
        quL=len(query)
        print("Query: {}".format(query))
        #Map the query sequence onto each read in the SAM
        pdSam['minD'], pdSam['minPos'], pdSam['matchseq'], pdSam['context'], pdSam['matchlist']  = zip(*p.map(distCompute, [(c,query) for c in pdSam.seq]))
        print("2. Adjusting minimum distance hit position on target...")
        pdSam['minD_v'], pdSam['minPos_v'], pdSam['matchseq_v'] , pdSam['context_v'],pdSam['matchlist_v'] = zip(*p.map(distCompute,[(templateSeq,c.replace("\|","")) for c in pdSam.context]))
        pdSam['insPos_v']=pdSam["minPos_v"]+int(quL/2)
    return(pdSam)

if __name__=="__main__":
    query, templateSeq, pdSam, Dtemp = Initialize(sys.argv)
    pdSam = poolDists(query,pdSam)
    print("Mapped hits in {} reads.".format(len(pdSam[pdSam.minD<(Dtemp-1)])))
    pdSam.to_csv("SAM_table_annotated.csv")
    print("Wrote {}.".format("stickleback_SAM_table.csv"))
  

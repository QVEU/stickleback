#!/usr/bin/env python
'''
stickleback.py
  ><```º>
Patrick T. Dolan
Unit Chief, Quantitative Virology and Evolution Unit
6/4/22
USAGE: python stickleback.py pathto/input.sam query templateFasta [] []   ( [arg] = optional argument)
Designed to map inserted sequences in a template sequence (plasmid or transcript).
'''

##### Imports #####
import pandas as pd
import Levenshtein
import time
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
    print("\n\n----------------=============------------------")
    print("--==--==--==--==   ><```º>   ==--==--==--==--==")
    print("==--==--==--==-- stickleback --==--==--==--==--")
    print("----------------=============------------------\n\n")


    #Check for Args
    if len(args)<3:
        print("USAGE: python stickleback.py pathto/input.sam query templateFasta [minLength] [maxLength]   ( [arg] = optional argument)")
        exit()

    #Query
    query=args[2]
    print("Query Length: {}".format(len(query)))

    #Template File
    template=args[3]
    with open(template,'r') as fasta:
        templateSeq="".join([i.strip() for i in fasta][1:])
    print("Template Length: {}".format(len(templateSeq)))

    Dtemp=blockDist((templateSeq,query))[0] #compute the distance to the template (without insert) to determine sensitivity.
    print("Template-Query Distance: {}".format(Dtemp))

    #Choose Read Size
    if len(args)>4:
        minL=int(args[4]) # read length minimum
        maxL=int(args[5]) # read length maximum
    else:
        minL=len(query)
        maxL=len(templateSeq)

    cutoff=(int(Dtemp)-2)
    print("Mapping Reads of Size {} to {} with a cutoff of {}".format(minL, maxL, cutoff))

    #input SAM
    inputData=sys.argv[1]
    print("Loading SAM: {}".format(inputData))
    outfile=inputData.replace(".sam",'_stickleback.csv')
    print("Outfile: {}".format(outfile))
    pdSam=loadSAM(inputData,minL,maxL)
    return(str(query),str(templateSeq), pdSam, int(Dtemp), int(cutoff), outfile)

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
    print("Total Candidate Reads: {}".format(size))
    pdSam=pd.DataFrame(samInput,columns=names)
    pdSam['length']=pdSam.seq.apply(lambda s: len(s))
    pdSam=pdSam[(pdSam.length>minL)&(pdSam.length<maxL)&(pdSam.template!="*")]
    return(pdSam)

def blockDist(seq_query):
    '''
    ~ new in v0.1 ~
    Function: blockDist()
        Uses numpy vectorize to speed up the computation (formerly for-loop)

    Arguments:
        `seq_query` is a tuple of read and query string, passed this way for 'pool'

    '''

    #Unpack Values
    seq=seq_query[0]
    query=seq_query[1]
    #print(seq)
    #Prepare Query
    quL=len(query)
    minD=quL #set minimum distance to length of query

    #Make blocks of read sequence
    blocks=np.array([seq[i:(i+quL)] for i in range(len(seq)-quL)])

    squarer = lambda B: Levenshtein.distance(B,query)
    vfunc = np.vectorize(squarer,otypes=[str])
    bDist=vfunc(blocks)

    minD=min(bDist)

    #print(minD)
    minPos=np.argmin(bDist)
    matchlist=minPos
    meanPos=minPos.mean()
    minPos=int(meanPos)

    matchseq=blocks[minPos]
    contextSeq="|".join([blocks[minPos-1],blocks[min(len(blocks),minPos)]])
    return(minD,minPos,matchseq,contextSeq,matchlist)

def poolBlocks(query,pdSam,nthreads=8):#uses multithreading to compute the matches
    with Pool(nthreads) as p:

        print("1. Computing minimum distance hit position for {} reads.".format(len(pdSam)))
        quL=len(query)
        print("Query: {}".format(query))
        pdSam['minD'],   pdSam['minPos'],   pdSam['matchseq'],    pdSam['context'],   pdSam['matchlist']   = zip(*p.map(blockDist, [(c,query) for c in pdSam.seq]))

        print("2. Mapping minimum distance site on template sequence...")
        pdSam['minD_v'], pdSam['minPos_v'], pdSam['matchseq_v'] , pdSam['context_v'], pdSam['matchlist_v'] = zip(*p.map(blockDist, [(templateSeq,c.replace("\|","")) for c in pdSam.context]))
        print("done.")
        pdSam['insPos_v']=pdSam["minPos_v"]+int(quL/2)

    return(pdSam)

if __name__=="__main__":
    t0=(time.time())
    query, templateSeq, pdSam, Dtemp, cutoff, outfileName = Initialize(sys.argv)
    pdSam = poolBlocks(query,pdSam)
    print("Mapped hits in {} reads.".format(np.sum([int(i)<=int(cutoff) for i in pdSam.minD])))
    pdSam.to_csv(outfileName)
    t1=(time.time())
    print("Done in {} minutes.".format((t1-t0)/60))
    print("Wrote {}.".format(outfileName))

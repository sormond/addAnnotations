#!/usr/bin/env python3


# Function of script: takes info out of original vcf for variants given by varAft and adds to varAft annotation file. Identifies matching variants in vcf using chromosome and position info. Position info is also +1 or -1, to correct for complex variants where positional info changes by 1 when using varAft.
# Adds refseq ID

""" Code written by Shannon Ormond 2018. Contact email: s.ormond@massey.ac.nz """


# loads varAft annotation file and adds extra columns to be annotated

def loadf(file) :

    df = pd.read_csv(file, sep = "\t", comment='#')

    df['Refseq'] = np.nan

    return df



# loads vcf file

def loadvcf(file) :

    df = pd.read_csv(file, sep = "\t", comment='#', header=None)

    return df



# makes list containing varAft annotation dataframe chromosome and position 'chunks' e.g. [[1, 540649], [5, 3284579]]

def makelist(df) :

    txtlist = []

    for line in range(0, len(df)) :

        e = df.iloc[line, 0:2].tolist()

        e.append(line)

        f = []

        g = []

        f.append(df.iloc[line, 0])

        g.append(df.iloc[line, 0])

        f.append(df.iloc[line, 1] + 1)

        g.append(df.iloc[line, 1] - 1)

        f.append(line)

        g.append(line)

        txtlist.append(e)

        txtlist.append(f)

        txtlist.append(g)

    return(txtlist)



# compares chromsome and positions from varAft dataframe with same thing from vcf file. Where they match, INFO field from vcf is added to INFO column of varAft dataframe. Also, Exac freq and Refseq ID is extracted from INFO field and added to Refseq or Exac column

# writes out varAft dataframe with added info to a csv

def compare(df, fdf, list) :

    newdf = fdf

    for line in range(0, len(df)) :

        if line % 50000 == 0 :

                print(str(line) + " lines out of " + str(len(df)) + " completed.")

        chr = df.iloc[line, 0]

        chr = chr.replace("chr", "")

        pos = df.iloc[line, 1]

        vcflist = [chr, pos]

        for block in range(0, len(list)) :

            v = list[block]

            if v[0:2] == vcflist :

                d = v[2] #takes row of VarAft dataframe to edit

                INFO = df.iloc[line, 7]

                refseq = re.search('NM_(.+?)\|', INFO)

                if refseq :

                    refseq = refseq.group(1)

                    refseq = "NM_" + refseq

                    newdf.loc[d, 'Refseq'] = refseq

                exac = re.search('exac_AF=(.+?);', INFO)

                if exac :

                    exac = exac.group(1)

                    newdf.loc[d, 'ExacFreq'] = exac

                newdf.loc[d, 'INFO'] = INFO

    newdf.to_csv(args.datafile + ".csv", sep=',', mode='w')



# required package stuff

import pandas as pd

import sys

import argparse

import re

import numpy as np



# command line stuff

parser = argparse.ArgumentParser()

parser.add_argument('-1', '--vcf1', help = "Required argument. CSV/text file to add annotations to.")

parser.add_argument('-2', '--vcf2', help = "Required argument. VCF file to take annotations from.")

parser.add_argument('-d', '--datafile', help = "Required argument. Stem name for output csv file.")

args = parser.parse_args()



# calling function stuff

df1 = loadf(args.vcf1)

df2 = loadvcf(args.vcf2)

list = makelist(df1)

compare(df2, df1, list)
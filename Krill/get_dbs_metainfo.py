#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  5 11:26:42 2021

@author: saulo
"""
import pandas as pd
import pathlib, os


def get_csvs(path,pattern):
    return list(pathlib.Path(path).glob('**/{}'.format(pattern)))

def get(path,ext):
    ref_rename = get_csvs(path,'fastaFilesRenamed.tsv')
    rename = {}
    
    for ref in ref_rename:
        db = os.path.basename(os.path.dirname(ref))
        df = pd.read_csv(ref,dtype='object',sep='\t')
        df = df.apply(lambda col: col.map(lambda x: x.rstrip(f'.{ext}')))
        rename[db] = pd.Series(df.OriginalName.values, index=df.NewName.values).to_dict()
        
    files = ['BGCs_with_Hits.tsv','Complete_BGCs_with_Hits.tsv','DNABases_and_ORFs_count.tsv']
    
    for f in files:
        tables = get_csvs(path,f)
        table_clean = pd.DataFrame()
        countORFsAndDNA = pd.DataFrame()
        
        for t in tables:
            db = os.path.basename(os.path.dirname(t))
            df = pd.read_csv(t,dtype='object',sep='\t')
            df['Database'] = db

            if not 'DNABases' in str(t):
                df['OriginalName'] = df['contig'].str.split('_').str[0].replace(rename[db],regex=True)
                table_clean = table_clean.append(df)
            
            if 'DNABases' in str(t):
                df[['NT','ORFS']] = df[['NT','ORFS']].astype(int)
                df = df.groupby('Database').sum()
                bgcs_count = len(list(pathlib.Path(os.path.join(path,db)).rglob('*region*.gbk')))
                df['BGCs/MegaBases'] = bgcs_count/df['NT']*1000000
                df['BGCs/MegaORFs'] = bgcs_count/df['ORFS']*1000000
                countORFsAndDNA = countORFsAndDNA.append(df)
                
        if not 'DNABases' in str(f):
            table_clean = table_clean[["Database","OriginalName","contig","cluster_number","product","completeness","SMILES","Start","End","Size","strand","genes","regulatory_genes",'KnownResistenceHit','CoreHit','BGCs_Hits','BGCs_Hits_Mean_Similarities(%)']].sort_values(by='Database')
            table_clean.to_csv(os.path.join(path,'DBs_'+f),index=False,sep='\t')
        if 'DNABases' in str(f):
            countORFsAndDNA.to_csv(os.path.join(path,'DBs_normalized_info.tsv'),sep='\t')

if __name__ == '__main__':
    get(os.getcwd(),'.fasta')
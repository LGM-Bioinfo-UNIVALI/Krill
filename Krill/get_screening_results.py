#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 23 19:24:23 2021

@author: saulo
"""
import numpy as np
import pandas as pd
import sys, os, ast
from statistics import mean


def df2xlsx(path, sheet_name, df):
    df = df.copy().reset_index()
    writer = pd.ExcelWriter(path, engine='xlsxwriter')

    df.to_excel(writer, sheet_name=sheet_name, startrow=1, header=False, index=False)
    worksheet = writer.sheets[sheet_name]

    (max_row, max_col) = df.shape
    column_settings = [{'header': column} for column in df.columns]
    worksheet.add_table(0, 0, max_row, max_col - 1, {'columns': column_settings})
    worksheet.set_column(0, max_col - 1, 15)

    writer._save()


def run(path):
    if os.path.exists(os.path.join(path,'ARTS','ARTS_Extractor','KnownHits.tsv')):
        KnownHits = pd.read_csv(os.path.join(path,'ARTS','ARTS_Extractor','KnownHits.tsv'),sep='\t',converters={'Sample':str})
    CoreHits = pd.read_csv(os.path.join(path,'ARTS','ARTS_Extractor','CoreHits.tsv'),sep='\t',converters={'Sample':str})
    DupHits = pd.read_csv(os.path.join(path,'ARTS','ARTS_Extractor','DupHits.tsv'),sep='\t',converters={'Sample':str})
    
    try:
        clusters_blast = pd.read_csv(os.path.join(path,'AntiSMASH','AntiSMASH_Extractor','clusters_blast.tsv'),sep='\t')
    except pd.errors.EmptyDataError:
        clusters_blast = pd.DataFrame()
    clusters = pd.read_csv(os.path.join(path,'AntiSMASH','AntiSMASH_Extractor','clusters.tsv'),sep='\t',converters={'file_name':str})

    hit_selected = clusters[['file_name','contig','cluster_number','product','contig_edge','SMILES','Start','End','strand','genes']].set_index('contig')
    hit_selected['KnownResistenceHit'] = ''
    hit_selected['BGCs_Hits'] = ''
    hit_selected['BGCs_Hits_Mean_Similarities(%)'] = ''
    for i in clusters.index:
        bgc = clusters.loc[i, 'contig']
        cluster_number = clusters.loc[i,'cluster_number']
        start = int(clusters.loc[i,'Start'])
        end = int(clusters.loc[i,'End'])
        
        if os.path.exists(os.path.join(path,'ARTS','ARTS_Extractor','KnownHits.tsv')):
            tmp = KnownHits[(KnownHits['Contig'] == bgc) & (KnownHits['HitStart'].astype(int) >= start) & (KnownHits['HitEnd'].astype(int) <= end)]
            if not tmp.empty:
                row_list = []
                for row in tmp.itertuples():
                    row_list.append(['KnownResistenceHit',row.Model,row.Description,row.evalue,row.HitStart,row.HitEnd,row.HitStrand])
                hit_selected.loc[((hit_selected.index == bgc) & (hit_selected.Start >= start) & (hit_selected.End <= end)),'KnownResistenceHit'] = str(row_list)
            else:
                hit_selected.loc[((hit_selected.index == bgc) & (hit_selected.Start >= start) & (hit_selected.End <= end)),'KnownResistenceHit'] = 'N/A'
        
        if not clusters_blast.empty:
            tmp_blast = clusters_blast[(clusters_blast['contig'] == bgc) & (clusters_blast['cluster'] == cluster_number)]
        else:
            tmp_blast = pd.DataFrame()
            
        if not tmp_blast.empty:
            row_list = []
            mean_similarities = []
            for row in tmp_blast.itertuples():
                row_list.append(['BGCs_Hits',row.accession,row.cluster_type,row.description,row.similarity])
                mean_similarities.append(row.similarity)
            hit_selected.loc[((hit_selected.index == bgc) & (hit_selected.cluster_number == cluster_number)),'BGCs_Hits'] = str(row_list)
            hit_selected.loc[((hit_selected.index == bgc) & (hit_selected.cluster_number == cluster_number)),'BGCs_Hits_Mean_Similarities(%)'] = round(mean(mean_similarities),2)
        else:
            hit_selected.loc[((hit_selected.index == bgc) & (hit_selected.cluster_number == cluster_number)),'BGCs_Hits'] = 'N/A'
            hit_selected.loc[((hit_selected.index == bgc) & (hit_selected.cluster_number == cluster_number)),'BGCs_Hits_Mean_Similarities(%)'] = 'N/A'

    hit_selected['CoreHit'] = ''
    for i in clusters.index:
        bgc = clusters.loc[i, 'contig']
        start = int(clusters.loc[i,'Start'])
        end = int(clusters.loc[i,'End'])
        tmp = CoreHits[(CoreHits['Contig'] == bgc) & (CoreHits['HitStart'].astype(int) >= start) & (CoreHits['HitEnd'].astype(int) <= end)]
        if not tmp.empty:
            row_list = []
            for row in tmp.itertuples():
                row_list.append(['CoreHit',row.Core_gene,row.Description,row.Function,row.HitStart,row.HitEnd,row.HitStrand])
            hit_selected.loc[((hit_selected.index == bgc) & (hit_selected.Start >= start) & (hit_selected.End <= end)),'CoreHit'] = str(row_list)
        else:
            hit_selected.loc[bgc,'CoreHit'] = 'N/A'
    
    hit_selected = hit_selected.replace('N/A',np.NaN).replace('',np.NaN)
    hit_selected.dropna(subset=['KnownResistenceHit','CoreHit'],how='all',inplace=True)
    
    
    hit_selected.rename(columns={'contig_edge':'completeness','file_name':'sample'},inplace=True)
    completeness = {"['True']":"Fragmented","['False']":"Complete"}
    hit_selected['completeness'] = hit_selected['completeness'].replace(completeness)
    hit_selected['Size'] = hit_selected['End']-hit_selected['Start']

    hit_selected['regulatory_genes'] = hit_selected['genes'].str.contains('regulatory',case=False)
    hit_selected = hit_selected[['sample','cluster_number', 'product', 'completeness', 'SMILES', 'Start', 'End', 'Size', 'strand', 'genes', 'regulatory_genes','KnownResistenceHit','CoreHit','BGCs_Hits','BGCs_Hits_Mean_Similarities(%)']]

    hit_selected.to_csv(os.path.join(path,'BGCs_with_Hits.tsv'),sep='\t')
    
    # hit_selected.to_csv(os.path.join(path,'BGCs_with_Hits.tsv'),sep='\t')
    df2xlsx(os.path.join(path,'BGCs_with_Hits.xlsx'), 'BGCs with Hits', hit_selected)

    hit_selected = hit_selected[hit_selected['completeness'] == "Complete"]
    hit_selected.to_csv(os.path.join(path,'Complete_BGCs_with_Hits.tsv'),sep='\t')
    # hit_selected.to_csv(os.path.join(path,'Complete_BGCs_with_Hits.tsv'),sep='\t')
    df2xlsx(os.path.join(path,'Complete_BGCs_with_Hits.xlsx'), 'Complete BGCs with Hits', hit_selected)

if __name__ == '__main__':
    run('/media/bioinfo/6tb_hdd/03_ELLEN/krill_runs/SRR13515398/')
# -*- coding: utf-8 -*-
#!/usr/bin/env python

import requests
import json
import pandas as pd
import sys
import re,os

def upperGenes(genes):
    # The app uses uppercase gene symbols. So it is crucial to perform upperGenes() step.
    return [gene.upper() for gene in genes]

# cosine distance search example
def query(genes,dg_value):
    data = {"genes":genes,"vals":dg_value}
    data['genes'] = upperGenes(data['genes'])
    config = {"aggravate":False,"searchMethod":"CD","share":False,"combination":False,"db-version":"latest"}
    metadata = [{"key":"Tag","value":"CD python example"},{"key":"Cell","value":"VCAP"}]
    payload = {"data":data,"config":config,"meta":metadata}
    headers = {'content-type':'application/json'}
    r = requests.post(url,data=json.dumps(payload),headers=headers)
    resCD2= r.json()  
    return resCD2

# extract
def compound_score_df(genes,dg_value):
    dg_value = query(genes,dg_value)
    top_meta = dg_value["topMeta"]
    score,pubchem_ID,sig_id = [],[],[]
    for element in top_meta:
        sig_id.append(element["sig_id"])
        score.append(element["score"])
        try:
            pubchem_ID.append(element["pubchem_id"])
        except:
            pubchem_ID.append("NaN")
    
    compound_score = { 
                        #'sig_id':sig_id,
                        'pubchem_ID':pubchem_ID,
                        'score':score
                        }
    compound_score = pd.DataFrame(compound_score)
    return compound_score


if __name__ == "__main__":
    url = 'http://amp.pharm.mssm.edu/L1000CDS2/query'
    input_df = sys.argv[1]
    output_df = sys.argv[2]
    df = pd.read_csv(input_df)
    #df = pd.read_csv("/home/tangchen/public/png/deGene.csv")
    # genes=["DDIT4","HIG2","FLT1","ADM","SLC2A3","ZNF331"]
    # dg_value=[9.97,10.16,7.66,17.80,20.29,15.22]
    genes = list(df.iloc[:,0])
    dg_value = list(df.iloc[:,1])
    compound_score = compound_score_df(genes,dg_value) 
    compound_score.to_csv(output_df)
    print(input_df,"is done^_^")                 


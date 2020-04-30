# -*- coding: utf-8 -*-


from Bio import Entrez
from bs4 import BeautifulSoup as BS
import re
import sys
sys.path.append("/home/aimzez/work/omics/pyscripts")
import pandas as pd
import numpy as np
import ncbilib

#%%
# USER VARIABLES
Entrez.email = "amsesk@umich.edu"
BIOPROJECT = "PRJNA510147"
LIBRARY_SOURCE_TARGET = "GENOMIC"


#%%
def print_sratoolkit_commands(df, outdir = "."):
    for t in df.itertuples():
        print (f"prefetch -O {outdir} {t.sra_accession}")
        

#%% Objectified so far
import ncbilib
runs = ncbilib.NcbiSraExperimentPackageSet.from_bioproject(BIOPROJECT)
runs = runs.filter("LIBRARY_SOURCE == GENOMIC")
runs.fetch_metadata()
#%%
runs.show_prefetch()
runs.show_rename_dump_full()
df = runs.to_pandas()
print(df.sort_values("BioSample"))
#%%
for line in open("/home/aimzez/work/pursuit_paper/jgi_search.txt").readlines():
    handle = Entrez.esearch(db="biosample", term=line.strip())
    records = Entrez.read(handle)
    handle.close()
    samples = BS(Entrez.efetch("biosample", id=records["IdList"]).read(), "xml")

    for i in samples("BioSample"):
        biosample = i.attrs["accession"]
        handle = Entrez.esearch(db="sra", term=biosample)
        sraids = Entrez.read(handle)["IdList"]
        handle.close()
        
        handle = Entrez.efetch(db="sra", id=sraids)
        sra = BS(handle.read(), "xml")
        handle.close()
        sra_gen = [x for x in sra("EXPERIMENT_PACKAGE") if x.LIBRARY_SOURCE.string == "GENOMIC"]
        if len(sra_gen) != 0:
            print(f"{'-'*40}\n{line.strip()}\n{'-'*40}")
            print (f"mkdir {line.replace(' ','_')}")
            print (f"cd {line.replace(' ','_')}")
            for exp in sra_gen: 
                out = [
                    exp.EXPERIMENT.attrs["accession"], 
                    exp.EXPERIMENT.PLATFORM.string, 
                    i.Owner.find("Name").string, 
                    exp.LIBRARY_LAYOUT
                    ]
                print(f"prefetch -O . {out[0]}")
                
            print("------------------------------------")
            

    
#print(samples.prettify())
#%%









#%%
##############################################
################# OLD CODE ###################
##############################################
#%% ENTREZ HANDLE
bs_to_strain = {acc: run.BioSample for acc,run in runs}
handle = Entrez.efetch(db="biosample", id = list(bs_to_strain.values()), maxret=100)
sample_xml = BS(handle.read(), "xml")
handle.close()
#%%
for acc,bs in bs_to_strain.items():
    s = sample_xml.select(f'BioSample[accession={bs}]')
    s = s[0]
    try:    
        strain = s.select('Attribute[attribute_name="strain"]')[0].string
    except IndexError:
        strain = s.select('Attribute[attribute_name="Strain"]')[0].string
    try:
        host = s.select('Attribute[display_name="host"]')[0].string
    except IndexError:
        host = np.nan
    bs_to_strain[acc] = (bs, strain, host)
    
'''
for s in (sample_xml.find_all("BioSample")):
    acc = s.attrs["accession"]
    try:    
        strain = s.select('Attribute[attribute_name="strain"]')[0].string
    except IndexError:
        strain = s.select('Attribute[attribute_name="Strain"]')[0].string
    try:
        host = s.select('Attribute[display_name="host"]')[0].string
    except IndexError:
        host = np.nan
    
    bs_to_strain[acc] = strain
'''
    
#%%
#genomic_runs = [x for x in runs if x.LIBRARY_SOURCE.string == "GENOMIC"]
rundict = {
    "sra_accession": None,
    "biosample": None,
    "strain": None
    }
runs_ldict = []
for run in runs:
    experiment = run.find("EXPERIMENT")
    rd = dict(rundict)
    rd["sra_accession"] = experiment.attrs["accession"]
    
    assert len(run.find_all("EXPERIMENT")) == 1, f"More than one SRA experiment associated with {rd['sra_accession']}"
    
    #print(run.LIBRARY_SOURCE.string)
    if run.LIBRARY_SOURCE.string == LIBRARY_SOURCE_TARGET:
        biosample = get_biosample(run.select("EXTERNAL_ID"))

        if biosample is not None:
            rd["biosample"] = biosample
            runs_ldict.append(rd)
        
        else:
            
            try:
                strain = run.select('EXTERNAL_ID[label]')[0].string.split(" ")[2]
                rd["strain"] = strain
                runs_ldict.append(rd)
                print(f"[WARNING]: No BioSample associated with SRA record, inferring `{strain}` from label: {rd['sra_accession']}")
                
            except IndexError:
                print(f"[SKIPPING] No BioSample associated with SRA record, and CANNOT infer from label: {rd['sra_accession']}")
            #'''
    
    else:
        pass
    
df = pd.DataFrame(runs_ldict)

#%%
df2 = pd.DataFrame(ldict)
df2 = pd.DataFrame(ldict).drop_duplicates()

#%%
df = df.merge(df2, on = "biosample", how="left")
df["strain"] = df.strain_x.combine(df.strain_y, func = lambda x,y: x if pd.isna(y) else y )
df = df.drop(["strain_x", "strain_y"], axis=1)
df = df.sort_values("strain").reset_index()

#%%
print_sratoolkit_commands(df)

#%%
'''
new = list(zip(acc, samples, strains))
for samp in new:
    print(f"prefetch -O . {samp[0]}")

#%% ENTREZ HANDLE
handle = Entrez.esearch(db="bioproject", term="PRJNA193498", retmax=1000)
records = Entrez.read(handle)
bp_ids = records["IdList"]
handle.close()
#%%

samples = [x.LocusTagPrefix.attrs["biosample_id"] for x in proj_records if x.LocusTagPrefix is not None and "biosample_id" in x.LocusTagPrefix.attrs]
print(len(samples))
#proj_LTP = [(x.LocusTagPrefix.attrs["biosample_id"], x.LocusTagPrefix.string) for x in proj_records]
#handle.close()
#%%
for i,run in enumerate(records.find_all("RUN_SET")):
    out = [run.IDENTIFIERS.string,  run.Member.attrs["organism"], run.Member.attrs["sample_title"][2:]]
    print (f"mv {out[0]}_1.fastq {out[1].replace(' ','_')}_{out[2]}_1_{i+1}.fastq")
    print (f"mv {out[0]}_2.fastq {out[1].replace(' ','_')}_{out[2]}_2_{i+1}.fastq")
'''
#%%

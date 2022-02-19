#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 22 15:03:34 2020

@author: aimzez
"""
#%%
from Bio import Entrez
from bs4 import BeautifulSoup as BS
import re
import sys
sys.path.append("/home/aimzez/dev/omics/pyscripts")
import pandas as pd
import numpy as np
import ncbilib

Entrez.email = "amsesk@umich.edu"
search_term = '"Mollicutes"[Organism]'

#%%
handle = Entrez.esearch(db="genome", term=search_term, retmax=500)
records = Entrez.read(handle)
handle.close()

for i in records["IdList"]:
    sample = BS(Entrez.esummary(db = "genome", id=i, report="full"), "xml")
    try:
        acc = sample.find(Name = "Assembly_Accession").string
        sp = sample.find(Name="Organism_Name").string
        
        ass_handle = Entrez.esearch(db="assembly", term=acc, retmax=10)
        ass_records = Entrez.read(ass_handle)
        ass_handle.close()
        assert len(ass_records["IdList"]) == 1, "error"
        ass = BS(Entrez.esummary(db = "assembly", id=ass_records["IdList"][0], report="full"), "xml")
        ass_ftppath = ass.FtpPath_GenBank.string
        
        biosample = BS(Entrez.efetch(db = "biosample", id=ass.BioSampleId.string, report="full"), "xml")
        biosample_
        print(f"{sp}\t{i}\t{acc}\t{ass_ftppath}")
    except:
        print(f"{sp}\t{i}\tNo Accession")
    
    
        
#%%
        
    439 GvMRE_I1
    354 GvMRE_I2
    141 GvMRE_Ic1
     80 GvMRE_Ic2
      1 GvMRE_Ic2S461
      1 GvMRE_Ic2S471
     79 GvMRE_Ic3
     94 GvMRE_Ic4
      1 GvMRE_Ic4S513
      1 GvMRE_Ic4S517
     49 GvMRE_Ic5
      1 GvMRE_Ic5S550
     29 GvMRE_Ic6
    287 GvMRE_II

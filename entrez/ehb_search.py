from Bio import Entrez
from bs4 import BeautifulSoup as BS
import re
import sys
sys.path.append("/home/aimzez/dev/omics/entrez")
import pandas as pd
import numpy as np
import ncbilib

search_terms = ['"uncultured Kocuria"[Organism]']

Entrez.email = "amsesk@umich.edu"

#%%
handle = Entrez.esearch(db="nucleotide", term=search_terms[0], retmax=10000)
records = Entrez.read(handle)
handle.close()

for i in records["IdList"]:
    print(i)

    '''
    sample = BS(Entrez.esummary(db="genome", id=i, report="full"), "xml")
    try:
        acc = sample.find(Name="Assembly_Accession").string
        sp = sample.find(Name="Organism_Name").string

        ass_handle = Entrez.esearch(db="assembly", term=acc, retmax=10)
        ass_records = Entrez.read(ass_handle)
        ass_handle.close()
        assert len(ass_records["IdList"]) == 1, "error"
        ass = BS(Entrez.esummary(db="assembly", id=ass_records[
                 "IdList"][0], report="full"), "xml")
        ass_ftppath = ass.FtpPath_GenBank.string

        biosample = BS(Entrez.efetch(
            db="biosample", id=ass.BioSampleId.string, report="full"), "xml")
        biosample_
        print(f"{sp}\t{i}\t{acc}\t{ass_ftppath}")
    except:
        print(f"{sp}\t{i}\tNo Accession")
    '''

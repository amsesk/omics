from Bio import Entrez
from bs4 import BeautifulSoup as BS
import re
import sys

sys.path.append("/home/aimzez/dev/omics/entrez")
import pandas as pd
import numpy as np
import ncbilib

fungal_tax = [
    "Fungi", "Mucoromycota", "Ascomycota", "Basidiomycota", "Zoopagomycota",
    "Chytridiomycota", "Cryptomycota", "Blastocladiomycota",
    "Neocallimastigomycota", "Mortierellomycotina", "Mucoromycotina",
    "Glomeromycotina", "Mortierella", "Rhizopus", None
]
bacterial_tax = ["Burkholderiaceae", "Mollicutes", "Bacteria"]

modifiers = [None, "endosymbiotic", "endohyphal"]

Entrez.email = "amsesk@umich.edu"


def get_host(record):
    gbquals = record.find_all("GBQualifier")
    for g in gbquals:
        names = g.find("GBQualifier_name").contents
        assert len(names) == 1, "Not one name for GBQUalifier"
        if names[0] == "host":
            return g.find("GBQualifier_value").contents[0]


#%%
for b in bacterial_tax:
    for f in fungal_tax:
        for m in modifiers:
            s = f'"{b}"[ORGN] AND "environmental samples"[ORGN]'

            if f is not None:
                s = s + f' AND "{f}"[ALL]'

            if m is not None:
                s = s + f' AND "{m}"[ALL]'

            print(f"SEARCH TERM: {s}")
            handle = Entrez.esearch(db="nucleotide",
                                    term=s,
                                    retmax=1000000,
                                    retmode="xml")
            try:
                records = Entrez.read(handle)
            except ValueError:
                print(f"SEARCH TERM FAILED: {s}")

            handle.close()

            count = 0
            for i in records["IdList"]:
                is_rdna = False
                handle = Entrez.efetch(db="nucleotide", id=i, retmode="xml")
                record = BS(handle.read(), "xml")
                handle.close()

                features = record.find_all("GBFeature")
                for ft in features:
                    if ft.find("GBFeature_key").contents[0] == "rRNA":
                        is_rdna = True
                        break

                if is_rdna:
                    host = get_host(record)

                    if host is None:
                        continue
                    else:
                        accession = record.find(
                            "GBSeq_accession-version").contents[0]
                        organism = record.find("GBSeq_organism").contents[0]
                        taxonomy = record.find("GBSeq_taxonomy").contents[0]

                        print('\t'.join([accession, host, organism, taxonomy]))

        # print(record.prettify())

        # print(record.find("GBSeq_locus").contents[0])
        # sys.exit()
        # print(BS(record, "xml").prettify())

        # print(b, f, count)
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

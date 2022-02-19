#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 11:23:26 2020

@author: aimzez
"""
import re
import ast
import numpy as np
import pandas as pd
from Bio import Entrez
from bs4 import BeautifulSoup as BS


class NcbiNucleotideRecords (object):

    def __init__(self):
        pass

    def from_keyword(soup):


class NcbiSraExperimentPackageSet (object):

    def __init__(self, index):
        self.index = index
        self.META = ["sample_name", "geo_loc_name",
                     "isolation_source", "strain", "host"]

    def __iter__(self):
        yield from self.index.items()

    def from_soup(soup, delim="EXPERIMENT_PACKAGE"):
        index = dict()
        for run in soup.find_all(delim):

            acc = run.find("EXPERIMENT").attrs["accession"]

            if acc not in index:
                index[acc] = NcbiSraExperimentPackage(acc, run)
            else:
                raise KeyError(f"Experiment {exp} duplicated.")

        return NcbiSraExperimentPackageSet(index)

    def from_bioproject(bioproject, delim="EXPERIMENT_PACKAGE"):
        handle = Entrez.esearch(db="sra", term=bioproject, retmax=1000)

        records = Entrez.read(handle)
        sra_ids = records["IdList"]
        handle.close()

        handle = Entrez.efetch(db="sra", id=sra_ids)
        sra_records = BS(handle.read(), "xml")
        handle.close()

        sra_accessions = [record.attrs["accession"]
                          for record in sra_records("EXPERIMENT")]

        handle = Entrez.efetch(db="sra", id=sra_accessions)
        sra_records = BS(handle.read(), "xml")
        handle.close()

        return NcbiSraExperimentPackageSet.from_soup(sra_records)

    def filter(self, expr):
        s = re.search("(.*)(==|!=)(.*)", expr)
        if s is None:
            raise ValueError("Invalid filtering expression.")

        else:
            attr, operand, value = (x.strip() for x in s.groups())

        new_index = dict()
        for acc, run in self:
            if not hasattr(run, attr):
                raise ValueError("Invalid attribute.")

            if operand == "==":
                if getattr(run, attr) == value:
                    new_index[acc] = run
            elif operand == "!=":
                if getattr(run, attr) != value:
                    new_index[acc] = run

        return NcbiSraExperimentPackageSet(new_index)

    def fetch_metadata(self):
        biosamples = [x.BioSample for _, x in self]
        handle = Entrez.efetch(db="biosample", id=biosamples, maxret=100)
        sample_xml = BS(handle.read(), "xml")
        handle.close()

        for _, run in self:
            s = sample_xml.select(f'BioSample[accession={run.BioSample}]')
            s = s[0]
            meta_dict = dict()
            for m in self.META:
                try:
                    meta_dict[m] = s.select(f'Attribute[harmonized_name="{m}"]')[0].string
                except IndexError:
                    meta_dict[m] = np.nan

            try:
                meta_dict["taxonomy_name"] = s.Description.Organism.attrs[
                    "taxonomy_name"]
            except:
                meta_dict["taxonomy_name"] = np.nan

            try:
                meta_dict["title"] = s.Description.Title.string
            except:
                meta_dict["title"] = np.nan

            run.metadata = meta_dict

        return None

    def to_ldict(self):
        ldict = []
        for _, run in self:
            run_dict = {
                "SRX": run.srx,
                "SRR": run.srr,
                "BioSample": run.BioSample
            }
            run_dict.update(run.metadata)

            ldict.append(run_dict)

        return ldict

    def to_pandas(self):
        return pd.DataFrame(self.to_ldict())[[
            "SRX",
            "SRR",
            "BioSample",
            "taxonomy_name",
            "strain",
            "title",
            "host",
            "isolation_source",
            "geo_loc_name",
            "sample_name"
        ]]

    def show_prefetch(self):
        for srx, _ in self:
            print(f"prefetch -O . {srx}")

    def show_rename_dump_strain(self):
        for _, run in self:
            print(f"mv {run.srr}_1.fastq {run.metadata['strain']}_1.fastq")
            print(f"mv {run.srr}_2.fastq {run.metadata['strain']}_2.fastq")

    def show_rename_dump_full(self):
        for _, run in self:
            print(f"mv {run.srr}_1.fastq {run.metadata['taxonomy_name'].replace(' ','_')}_{run.metadata['strain']}_1.fastq")
            print(f"mv {run.srr}_2.fastq {run.metadata['taxonomy_name'].replace(' ','_')}_{run.metadata['strain']}_2.fastq")


class NcbiSraExperimentPackage (object):

    def __init__(self, srx, soup):
        self.srx = srx
        self.soup = soup
        self.LIBRARY_SOURCE = self.soup.LIBRARY_SOURCE.string
        self.BioSample = self.get_biosample(soup.select("EXTERNAL_ID"))
        self.srr = self.get_srr()
        self.metadata = None

    def get_biosample(self, tags, attr="namespace"):
        if tags is None:
            return None

        for t in tags:
            namespace_lower = t.attrs["namespace"].lower()
            if namespace_lower == "biosample":
                return t.string

        print(f"[WARNING]: No BioSample associated with SRA record: {self.srx}")
        return None

    def get_srr(self):
        return self.soup.find("RUN").attrs["accession"]

if __name__ == "__main__":
    pass

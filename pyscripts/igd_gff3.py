#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import pandas as pd
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
#%%
if sys.argv[1] in ["-h", "--help", "-help"]:
    print("USAGE: python igd_gff3.py [path/to/gff3] [path/to/output]")
    sys.exit()
#%%
gff3 = pd.read_csv(sys.argv[1], sep="\t", comment="#", header=None)

gff3_sub = gff3.iloc[:,[0,3,4,6]]
gff3_sub.columns = ["contig", "start", "end", "strand"]

gff3_sub = gff3_sub.sort_values(["contig", "strand"])
gff3_sub["next_start"] = gff3_sub.groupby(["contig","strand"]).start.shift(-1)
gff3_sub = gff3_sub.assign(igd = gff3_sub["next_start"] - gff3_sub["end"])
gff3_sub.to_csv(sys.argv[2], sep="\t", index=False)
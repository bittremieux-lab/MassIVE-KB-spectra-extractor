#!/usr/bin/env python3
import sys

import pandas as pd


def process_mzml_group(tsv_file_path):
    pd.read_csv(tsv_file_path, sep="\t")


tsv_file_path = sys.argv[1]

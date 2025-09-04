#!/usr/bin/env python3

import re
import sys
import pandas as pd
from tqdm import tqdm


def group_tsv(tsv_file_path, group_column):
    tsv_df = pd.read_csv(tsv_file_path, sep="\t", nrows=10)
    print("Opened tsv file:", tsv_file_path)
    for column_val, grp in tqdm(tsv_df.groupby(group_column)):
        grp.to_csv(
            f"group_{re.sub(r'[^A-Za-z0-9._-]', '_', str(column_val))}.tsv",
            sep="\t",
            index=False,
        )


tsv_file_path = sys.argv[1]
group_column = sys.argv[2]
group_tsv(tsv_file_path, group_column)

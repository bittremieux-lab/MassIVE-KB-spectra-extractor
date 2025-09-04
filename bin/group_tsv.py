#!/usr/bin/env python3
import os
import re
import sys

import pandas as pd
from tqdm import tqdm


def group_tsv(tsv_file_path, group_column, groups_dir):
    os.makedirs(groups_dir, exist_ok=True)
    tsv_df = pd.read_csv(tsv_file_path, sep="\t")
    print("Opened tsv file:", tsv_file_path)
    for column_val, grp in tqdm(tsv_df.groupby(group_column)):
        out_file_path = os.path.join(
            groups_dir, f"group_{re.sub(r'[^A-Za-z0-9._-]', '_', str(column_val))}.tsv"
        )
        grp.to_csv(
            out_file_path,
            sep="\t",
            index=False,
        )


tsv_file_path = sys.argv[1]
group_column = sys.argv[2]
groups_dir = sys.argv[3]
group_tsv(tsv_file_path, group_column, groups_dir)

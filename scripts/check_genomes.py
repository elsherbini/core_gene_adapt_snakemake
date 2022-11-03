import os
import re
import pandas as pd
records = []

regex = r"C\:([\d\.]+)%\[S\:([\d\.]+)%,D\:([\d\.]+)%\],F\:([\d\.]+)%,M\:-?([\d\.]+)%,n\:(\d+)"
for fpath in snakemake.input:
    isolate = os.path.splitext(os.path.basename(fpath))[0].replace("short_summary_", "")
    with open(fpath, "r") as f:
        for line in f.readlines():
            if line[0] == "#" or len(line) == 0:
                continue
            if re.search(regex, line):
                complete, single_copy, duplicated, fragmented, missing, total = re.search(regex, line).groups()
                records.append(dict(zip(["isolate", "complete", "single_copy", "duplicated", "fragmented", "missing", "total"], [isolate, complete, single_copy, duplicated, fragmented, missing, total])))
                break
pd.DataFrame(records).to_csv(snakemake.output[0], index=False, columns=["isolate", "complete", "single_copy", "duplicated", "fragmented", "missing", "total"])

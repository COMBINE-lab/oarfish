import pandas as pd
import numpy as np
import json


def main():
    # read ground truth
    t = pd.read_csv("../../data/sequin_ref_data/rnasequin_isoforms_2.4.tsv", sep="\t")
    t.columns = t.columns.str.upper()

    data_dir = "../../data/sequin_ref_data"
    d = json.loads(open(f"{data_dir}/data_map.json").read())
    for k, v in d.items():
        d[k] = v.upper()

    quant_dir = "sequin"
    # read nanocount and oarfish output
    x = {}
    y = {}
    for k, v in d.items():
        v = v.upper()
        x[v] = pd.read_csv(f"{quant_dir}/quants/oarfish/{k}_quant.tsv", sep="\t")
        y[v] = pd.read_csv(f"{quant_dir}/quants/nanocount/{k}_quant.tsv", sep="\s+")

    m = {}
    for k, v in d.items():
        v = v.upper()
        m[v] = pd.merge(x[v], t, left_on="tname", right_on="NAME", how="inner")
        m[v] = pd.merge(
            y[v], m[v], left_on="transcript_name", right_on="NAME", how="inner"
        )
        m[v] = m[v].rename(columns={"num_reads" : "oarfish", "est_count" : "NanoCount"})

    res = []
    for k, v in m.items():
        mix = "MIX_B"
        if k in ["UNDIFF1", "UNDIFF2"]:
            mix = "MIX_A"
        print(f"\n=== {k} ===\n")
        res.append(
            v.loc[:, ["oarfish", "NanoCount", mix]].corr(method="spearman")[mix]
        )
        print(res[-1])


if __name__ == "__main__":
    main()

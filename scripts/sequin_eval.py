import pandas as pd
import numpy as np
import json


def main():
    # read ground truth
    t = pd.read_excel(
        "test_data/TableS3_S4_Sequin_counts_and_diff_expression.xlsx",
        sheet_name="S3 - Sequin-transcript-count",
    )
    t.columns = t.columns.str.upper()

    data_dir = "../data"
    d = json.loads(open(f"{data_dir}/data_map.json").read())
    for k, v in d.items():
        d[k] = v.upper()

    # read nanocount and oarfish output
    x = {}
    y = {}
    for k, v in d.items():
        v = v.upper()
        x[v] = pd.read_csv(f"{data_dir}/quants/oarfish/{k}_quant.tsv", sep="\t")
        y[v] = pd.read_csv(f"{data_dir}/quants/nanocount/{k}_quant.tsv", sep="\s+")

    m = {}
    m2 = {}
    for k, v in d.items():
        v = v.upper()
        m[v] = pd.merge(x[v], t, left_on="tname", right_on="TRANSCRIPT", how="inner")
        m2[v] = pd.merge(
            y[v], t, left_on="transcript_name", right_on="TRANSCRIPT", how="inner"
        )

    res_oarfish = []
    res_nanocount = []
    for k, v in m.items():
        if k in t.columns:
            print(f"\n\n=== {k} ===\n")
            res_oarfish.append(v.corr(method="spearman").loc[k, "num_reads"])
            print(res_oarfish[-1])
            v2 = m2[k]
            res_nanocount.append(v2.corr(method="spearman").loc[k, "est_count"])
            print(res_nanocount[-1])

    print(
        f"oarfish_mean = {np.array(res_oarfish).mean()}, nanocount_mean = {np.array(res_nanocount).mean()}"
    )


if __name__ == "__main__":
    main()

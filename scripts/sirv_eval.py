import pandas as pd
import numpy as np
import json


def main():
    # read ground truth
    t = pd.read_excel("../../data/sirv_ref_data/molar_concentrations.xlsx", sheet_name="Sheet2")
    t["E0"] = t["E2"]
    t.loc[t["E0"] > 0, "E0"] = 1

    d = {"SRR6058583": "E2", "SRR6058584": "E0"}

    data_dir = "sirv"
    for a in ["c", "i", "o"]:
        merged = {}
        for k, v in d.items():
            o = pd.read_csv(
                f"{data_dir}/quants/oarfish/{k}_{a}_quant.tsv", sep="\t"
            )
            n = pd.read_csv(
                f"{data_dir}/quants/nanocount/{k}_{a}_quant.tsv", sep="\s+"
            )
            merged[v] = pd.merge(o, t, left_on="tname", right_on="Name", how="outer")
            merged[v] = pd.merge(
                n, merged[v], left_on="transcript_name", right_on="Name", how="outer"
            )
        print(f"\n============== reference {a} ==================\n")
        print(f"=== E2 ===\n")
        print(
            merged["E2"][["E2", "est_count", "num_reads"]]
            .corr(method="spearman")
            .loc["E2", :]
        )

        print(f"\n=== E0 ===\n")
        nm = merged["E0"].loc[merged["E0"].E0 == 1, "est_count"].mean()
        ns = merged["E0"].loc[merged["E0"].E0 == 1, "est_count"].std()
        om = merged["E0"].loc[merged["E0"].E0 == 1, "num_reads"].mean()
        os = merged["E0"].loc[merged["E0"].E0 == 1, "num_reads"].std()
        print(f"NanoCount CV: {ns/nm}\noarfish CV: {os/om}")

        nzs = merged["E0"].loc[merged["E0"].E0 == 0, "est_count"].sum()
        nas = merged["E0"].loc[:, "est_count"].sum()
        ozs = merged["E0"].loc[merged["E0"].E0 == 0, "num_reads"].sum()
        oas = merged["E0"].loc[:, "num_reads"].sum()
        print(f"NanoCount % reads mapped: {nzs/nas}\noarfish % reads mapped: {ozs/oas}")


if __name__ == "__main__":
    main()

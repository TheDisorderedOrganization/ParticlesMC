# /// script
# requires-python = ">=3.13"
# dependencies = [
#     "fastparquet>=2025.12.0",
#     "seaborn>=0.13.2",
# ]
# ///

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def main() -> None:
    df_here = pd.read_parquet("fskt.parquet")
    df_here["method"] = "ParticlesMC"

    temperatures = [3, 2, 1.6, 1.4, 1.25, 1.2, 1.15, 1.1, 1.05, 1]
    assert len(temperatures) == len(df_here["T"].unique())
    df_list = []
    for T in temperatures:
        df = pd.read_csv(
            f"4/fsqt_flip/temp_{T}.dat", sep=" ", skiprows=1, names=["t", "F_s(k,t)"]
        )
        df["T"] = T
        df_list.append(df)
    df_paper = pd.concat(df_list)
    df_paper["method"] = "Paper"

    df_tot = pd.concat([df_here, df_paper]).reset_index(drop=True)
    df_tot["Temperature"] = df_tot["T"].apply(lambda x: str(x))

    ax = sns.lineplot(
        data=df_tot,
        x="t",
        y="F_s(k,t)",
        hue="Temperature",
        style="method",
        markers={"ParticlesMC": "o", "Paper": ","},
    )
    ax.set_xscale("log")
    plt.savefig("fskt.png")
    plt.show()


if __name__ == "__main__":
    main()

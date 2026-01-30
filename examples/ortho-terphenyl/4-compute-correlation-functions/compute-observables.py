# /// script
# requires-python = ">=3.14"
# dependencies = [
#     "atooms-pp>=4.2.1",
#     "fastparquet>=2025.12.0",
#     "pandas>=3.0.0",
# ]
# ///

import atooms.postprocessing as pp
import numpy as np
import pandas as pd
from atooms.trajectory import Trajectory

temperatures = [3.0, 2.0, 1.6, 1.4, 1.25, 1.2, 1.15, 1.1, 1.05, 1.0]


def compute_fskt() -> pd.DataFrame:
    print("F_s(k,t)")
    df_list = []
    for T in temperatures:
        print(f"T = {T}")
        traj = Trajectory(f"../3-run-production/{T}/trajectories/1/trajectory.xyz")

        cf = pp.SelfIntermediateScatteringFast(
            traj,
            ksamples=1,
            kmin=7.4,
            kmax=7.4,
            tgrid=sorted({round(x) for x in np.logspace(0, 6, num=55)}),
        )
        cf.compute()

        print(f"k = {cf.grid[0][0]}")

        N = len(cf.grid[1])
        df_list.append(
            pd.DataFrame(
                {
                    "t": cf.grid[1],
                    cf.qualified_name: np.array(cf.values)[0],
                    "T": [T] * N,
                }
            )
        )

    return pd.concat(df_list).reset_index(drop=True)


def main() -> None:
    df = compute_fskt()
    df.to_parquet("fskt.parquet")


if __name__ == "__main__":
    main()

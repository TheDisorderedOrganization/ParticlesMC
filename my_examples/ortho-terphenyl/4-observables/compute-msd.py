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

temperatures = [3.0, 2.0, 1.5, 1.0]

def compute_msd() -> pd.DataFrame:
    print("msd(t)")
    df_list = []
    for T in temperatures:
        print(f"T = {T}")
        traj = Trajectory(f"../3-production/{T}/chains/1/trajectory.xyz")

        msd = pp.MeanSquareDisplacement(
            traj,
            rmax=-1.0
        )
        msd.compute()

        print(f"r = {msd.grid[0]}")

        N = len(msd.grid)
        df_list.append(
            pd.DataFrame(
                {
                    "t": msd.grid,
                    msd.qualified_name: np.array(msd.value),
                    "T": [T] * N,
                }
            )
        )

    return pd.concat(df_list).reset_index(drop=True)


def main() -> None:
    df = compute_msd()
    df.to_parquet("msd.parquet")


if __name__ == "__main__":
    main()
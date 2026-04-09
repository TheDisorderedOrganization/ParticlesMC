import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

df = pd.read_parquet("msd.parquet")
print(df.columns)

ax = sns.lineplot(
    data=df,
    x="t",
    y="dr^2(t)",
    hue="T",
)
# ax.set_xscale("log")
# ax.set_xlim(left=1)
# ax.set_yscale("log")
plt.xlabel("t (MC steps)")
plt.ylabel("msd(t)")
plt.savefig("msd.png")
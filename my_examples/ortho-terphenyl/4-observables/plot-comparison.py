import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

df = pd.read_parquet("fskt.parquet")
print(df.columns)

ax = sns.lineplot(
    data=df,
    x="t",
    y="F_s(k,t)",
    hue="T",
)
ax.set_xscale("log")
plt.xlabel("t (MC steps)")
plt.ylabel("F_s(k,t)")
plt.savefig("fskt.png")
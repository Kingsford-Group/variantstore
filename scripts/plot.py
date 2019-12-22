import pandas as pd
import numpy as np
import sys
import matplotlib.pyplot as plt

# plt.rcParams.update({'font.size': 28})

df = pd.read_csv(sys.argv[1], names=["count", "Time"], low_memory=True)
# print(df)

N = df["Time"].max()
print(N)
C = df["count"].max()
print(C)

#df =(df["count"].map(str) + "_" + df["time"].map(str)).value_counts()
# hist = df.plot.hist(by=df["count"])
bins = np.arange(0, 3000, 300)
print(bins)
# df = df.groupby(pd.cut(df["count"], bins, labels=labels), as_index=False)["time"].mean().plot(kind="bar", rot=0, legend=True)
df = df.groupby(pd.cut(df["count"], bins, include_lowest=True), as_index=False)["Time"].mean()
ax = df.plot.bar(rot=0)

labels=[[bins[i], bins[i+1]-1]
        for i in range(len(bins)-1)]
# list(map(list, zip(bins, bins[1:])))
print(labels)
ax.set_xticklabels(labels)
ax.tick_params(axis='both', which='major', labelsize=16)
ax.legend(loc="upper left", prop={'size': 28})

plt.ylabel("Time (secs)", fontsize=28)
plt.xlabel("Number of  variants", fontsize=28, labelpad=20)
plt.show()
# fig = hist.get_figure()
# fig.savefig(sys.argv[2])
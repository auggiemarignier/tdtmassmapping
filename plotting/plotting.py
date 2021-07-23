import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize
import sys
import os


directory = sys.argv[1]
truth = np.loadtxt(f"{directory}/truth.txt").reshape((256, 256))
mean = np.loadtxt(f"{directory}/mean.txt")
std = np.loadtxt(f"{directory}/stddev.txt")
khist = np.loadtxt(f"{directory}/khistogram.txt")
last_nonzero_k = np.argwhere(khist[:, 1]).max()
likelihoods = np.loadtxt(f"{directory}/likelihood.txt")

restarts = []
while os.path.isdir(f"{directory}/restart/"):
    restarts.append(likelihoods.shape[0])
    directory += "/restart"
    if os.path.isfile(f"{directory}/mean.txt"):
        mean += np.loadtxt(f"{directory}/mean.txt")
        std += np.loadtxt(f"{directory}/stddev.txt")  # note this assumes NO correlation between pixels
        likelihoods = np.concatenate([likelihoods, np.loadtxt(f"{directory}/likelihood.txt")])
        khist[:, 1] += np.loadtxt(f"{directory}/khistogram.txt", usecols=1)
        last_nonzero_k = np.argwhere(khist[:, 1]).max()

mean /= len(restarts) + 1
std /= len(restarts) + 1
diff = np.abs(truth - mean)

vmin = truth.min()
vmax = truth.max()
stdmin = std.min()
stdmax = std.max()
diffmin = diff.min()
diffmax = diff.max()

cmap = cm.inferno
diffcmap = cm.binary_r
stdcmap = cm.turbo

mosaic = """.ABE
            GCDF"""
fig = plt.figure(figsize=(8, 8))
axd = fig.subplot_mosaic(
    mosaic,
    gridspec_kw={"width_ratios": [1, 15, 15, 1], "wspace": 0.05},
)

axd["A"].imshow(truth, vmin=vmin, vmax=vmax, cmap=cmap)
axd["A"].set_title("Truth")
axd["A"].axis("off")

axd["B"].imshow(mean, vmin=vmin, vmax=vmax, cmap=cmap)
axd["B"].set_title("Mean TDT")
axd["B"].axis("off")

axd["C"].imshow(diff, cmap=diffcmap)
axd["C"].set_title("|Truth - Mean|")
axd["C"].axis("off")

axd["D"].imshow(std, cmap=stdcmap)
axd["D"].set_title("Standard Dev.")
axd["D"].axis("off")

fig.colorbar(
    cm.ScalarMappable(norm=Normalize(vmin=vmin, vmax=vmax), cmap=cmap),
    axd["E"],
    shrink=0.1,
)
fig.colorbar(
    cm.ScalarMappable(norm=Normalize(vmin=stdmin, vmax=stdmax), cmap=stdcmap),
    axd["F"],
    shrink=0.1,
)
fig.colorbar(
    cm.ScalarMappable(norm=Normalize(vmin=diffmin, vmax=diffmax), cmap=diffcmap),
    axd["G"],
    shrink=0.1,
)
axd["G"].yaxis.set_ticks_position('left')

plt.show()

mosaic = """A
            B"""
fig = plt.figure(figsize=(10, 5), constrained_layout=True)
axd = fig.subplot_mosaic(mosaic)

axd["A"].bar(khist[:last_nonzero_k, 0], khist[:last_nonzero_k, 1])
axd["A"].set_xlabel("Number of parameters")
axd["A"].set_ylabel("Count")

axd["B"].plot(likelihoods)
for r in restarts:
    axd["B"].axvline(r, c="k", ls="--")
axd["B"].set_yscale("log")
axd["B"].set_xlabel("Sample number")
axd["B"].set_ylabel("-log(likelihood)")

plt.show()

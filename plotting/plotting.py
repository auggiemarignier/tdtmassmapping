import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize
import sys


directory = sys.argv[1]
truth = np.loadtxt("../data/Bolshoi_7_clean_256.txt").reshape((256, 256))
mean = np.loadtxt(f"{directory}/mean.txt")
std = np.loadtxt(f"{directory}/stddev.txt")

vmin = truth.min()
vmax = truth.max()
stdmin = std.min()
stdmax = std.max()

cmap = cm.inferno
diffcmap = cm.binary
stdcmap = cm.turbo

mosaic = """ABE
            CDF"""
fig = plt.figure(figsize=(8, 8), tight_layout=True)
axd = fig.subplot_mosaic(
    mosaic,
    gridspec_kw={"width_ratios": [15, 15, 1], "wspace": 0.05},
)

axd["A"].imshow(truth, vmin=vmin, vmax=vmax, cmap=cmap)
axd["A"].set_title("Truth")
axd["A"].axis("off")

axd["B"].imshow(mean, vmin=vmin, vmax=vmax, cmap=cmap)
axd["B"].set_title("Mean TDT")
axd["B"].axis("off")

axd["C"].imshow(truth - mean, cmap=diffcmap)
axd["C"].set_title("Truth - Mean")
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

plt.show()

mosaic = """A
            B"""
fig = plt.figure(figsize=(10, 5), constrained_layout=True)
axd = fig.subplot_mosaic(mosaic)

khist = np.loadtxt(f"{directory}/khistogram.txt")
last_nonzero_k = np.argwhere(khist[:, 1]).shape[0]
axd["A"].bar(khist[:last_nonzero_k, 0], khist[:last_nonzero_k, 1])
axd["A"].set_xlabel("Number of parameters")
axd["A"].set_ylabel("Count")

likelihoods = np.loadtxt(f"{directory}/likelihood.txt")
axd["B"].plot(likelihoods)
axd["B"].set_yscale("log")
axd["B"].set_xlabel("Sample number")
axd["B"].set_ylabel("-log(likelihood)")

plt.show()

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize
import sys


directory = sys.argv[1]
truth = np.loadtxt("../data/Bolshoi_7_clean_256.txt").reshape((256, 256))
mean = np.loadtxt(f"{directory}/mean.txt")

vmin = truth.min()
vmax = truth.max()

cmap = cm.inferno

mosaic = """ABC
            D.."""
fig = plt.figure(figsize=(10, 5), constrained_layout=True)
axd = fig.subplot_mosaic(
    mosaic,
    gridspec_kw={"height_ratios": [15, 1], "wspace": 0.001},
)

axd["A"].imshow(truth, vmin=vmin, vmax=vmax, cmap=cmap)
axd["A"].set_title("Truth")
axd["A"].axis("off")

axd["B"].imshow(mean, vmin=vmin, vmax=vmax, cmap=cmap)
axd["B"].set_title("Mean TDT")
axd["B"].axis("off")

axd["C"].imshow(truth - mean, cmap="binary")
axd["C"].set_title("Truth - Mean")
axd["C"].axis("off")

fig.colorbar(
    cm.ScalarMappable(norm=Normalize(vmin=vmin, vmax=vmax), cmap=cmap),
    axd["D"],
    orientation="horizontal",
    shrink=0.1,
)
axd["D"].axis("off")
l, b, w, h = axd["D"].get_position(original=False).bounds
axd["D"].set_position([l + 0.15, 0.01, w, h])

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

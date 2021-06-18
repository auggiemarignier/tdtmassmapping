import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize

truth = np.loadtxt("../data/Bolshoi_7_clean_256.txt").reshape((256, 256))
mean = np.loadtxt("../outputs/mean.txt")

vmin = truth.min()
vmax = truth.max()

cmap = cm.inferno

mosaic = """AB
            C."""
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

fig.colorbar(
    cm.ScalarMappable(norm=Normalize(vmin=vmin, vmax=vmax), cmap=cmap),
    axd["C"],
    orientation="horizontal",
    shrink=0.1,
)
axd["C"].axis("off")
l, b, w, h = axd["C"].get_position(original=False).bounds
axd["C"].set_position([l + 0.15, 0.01, w, h])

plt.show()

mosaic = """A
            B"""
fig = plt.figure(figsize=(10, 5), constrained_layout=True)
axd = fig.subplot_mosaic(mosaic)

khist = np.loadtxt("../outputs/khistogram.txt")
axd["A"].bar(khist[:, 0], khist[:, 1])
axd["A"].set_xlabel("Number of parameters")
axd["A"].set_ylabel("Count")

likelihoods = np.loadtxt("../outputs/likelihood.txt")
axd["B"].plot(likelihoods)
axd["B"].set_xlabel("Sample number")
axd["B"].set_ylabel("-log(likelihood)")

plt.show()

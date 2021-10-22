import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
import sys
import os
from utils import meanvar_from_submeanvar


directory = sys.argv[1]
truth = np.loadtxt(f"{directory}/truth.txt").reshape((256, 256))
mean = np.loadtxt(f"{directory}/mean.txt")
std = np.loadtxt(f"{directory}/stddev.txt")
var = std ** 2
khist = np.loadtxt(f"{directory}/khistogram.txt")
best_fitting = np.loadtxt(f"{directory}/best_model.txt")
last_nonzero_k = np.argwhere(khist[:, 1]).max()
likelihoods = np.loadtxt(f"{directory}/likelihood.txt")
khistory = np.loadtxt(f"{directory}/khistory.txt")
current_n = likelihoods.shape[0]
restarts = []
while os.path.isdir(f"{directory}/restart/"):
    restarts.append(likelihoods.shape[0])
    directory += "/restart"
    if os.path.isfile(f"{directory}/mean.txt"):
        mean2 = np.loadtxt(f"{directory}/mean.txt")
        std2 = np.loadtxt(f"{directory}/stddev.txt")
        var2 = std2 ** 2
        likelihoods2 = np.loadtxt(f"{directory}/likelihood.txt")
        n_add = len(likelihoods2)

        mean, var = meanvar_from_submeanvar(
            mean,
            mean2,
            var,
            var2,
            current_n,
            n_add
        )
        current_n += n_add

        likelihoods = np.concatenate([likelihoods, likelihoods2])
        khist[:, 1] += np.loadtxt(f"{directory}/khistogram.txt", usecols=1)
        last_nonzero_k = np.argwhere(khist[:, 1]).max()
        best_fitting = np.loadtxt(f"{directory}/best_model.txt")
        khistory = np.concatenate([khistory, np.loadtxt(f"{directory}/khistory.txt")])

diff = np.abs(truth - mean)
diff_best = np.abs(truth - best_fitting)

vmin = truth.min()
vmax = truth.max()
stdmin = std.min()
stdmax = std.max()
diffmin = min([diff.min(), diff_best.min()])
diffmax = max([diff.max(), diff_best.max()])

cmap = cm.inferno
diffcmap = cm.binary_r
stdcmap = cm.turbo

mosaic = """.ABHE
            GCDIF"""
fig = plt.figure(figsize=(13, 9))
axd = fig.subplot_mosaic(
    mosaic,
    gridspec_kw={"width_ratios": [1, 15, 15, 15, 1], "wspace": 0.05},
)

axd["A"].imshow(truth, vmin=vmin, vmax=vmax, cmap=cmap)
axd["A"].set_title("Truth")
axd["A"].axis("off")

axd["B"].imshow(mean, vmin=vmin, vmax=vmax, cmap=cmap)
axd["B"].set_title("Mean TDT")
axd["B"].axis("off")

axd["D"].imshow(diff, vmin=diffmin, vmax=diffmax, cmap=diffcmap)
axd["D"].set_title("|Truth - Mean|")
axd["D"].axis("off")

axd["C"].imshow(std, cmap=stdcmap)
axd["C"].set_title("Standard Dev.")
axd["C"].axis("off")

axd["H"].imshow(best_fitting, cmap=cmap)
axd["H"].set_title("Best fitting model")
axd["H"].axis("off")

axd["I"].imshow(diff_best, vmin=diffmin, vmax=diffmax, cmap=diffcmap)
axd["I"].set_title("|Truth - Best fitting model|")
axd["I"].axis("off")

fig.colorbar(
    cm.ScalarMappable(norm=Normalize(vmin=vmin, vmax=vmax), cmap=cmap),
    axd["E"],
    shrink=0.1,
)
fig.colorbar(
    cm.ScalarMappable(norm=Normalize(vmin=stdmin, vmax=stdmax), cmap=stdcmap),
    axd["G"],
    shrink=0.1,
)
fig.colorbar(
    cm.ScalarMappable(norm=Normalize(vmin=diffmin, vmax=diffmax), cmap=diffcmap),
    axd["F"],
    shrink=0.1,
)
axd["G"].yaxis.set_ticks_position("left")

plt.show()

mosaic = """A
            B"""
fig = plt.figure(figsize=(10, 5), constrained_layout=True)
axd = fig.subplot_mosaic(mosaic)
prop_cycle = plt.rcParams["axes.prop_cycle"]
colours = prop_cycle.by_key()["color"]

axd["A"].bar(khist[:last_nonzero_k, 0], khist[:last_nonzero_k, 1])
axd["A"].set_xlabel("Number of parameters")
axd["A"].set_ylabel("Count")

axd["B"].plot(likelihoods, color=colours[0])
for r in restarts:
    axd["B"].axvline(r, c="k", ls="--")
if len(restarts) > 0:
    axd["B"].legend(
        handles=[Line2D([0], [0], linestyle="--", color="k", label="Restart")]
    )
axd["B"].set_yscale("log")
axd["B"].set_xlabel("Sample number")
axd["B"].set_ylabel("-log(likelihood)", color=colours[0])
axd["B"].spines["left"].set_color(colours[0])
axd["B"].tick_params(axis="y", color=colours[0], labelcolor=colours[0], which="both")

axd["C"] = axd["B"].twinx()
axd["C"].plot(khistory, color=colours[1])
axd["C"].set_ylabel("Number of parameters", color=colours[1])
axd["C"].spines["right"].set_color(colours[1])
axd["C"].spines["left"].set_visible(False)
axd["C"].tick_params(axis="y", color=colours[1], labelcolor=colours[1])

plt.show()

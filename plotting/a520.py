import numpy as np
import sys
import os
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1 import make_axes_locatable
from utils import meanvar_from_submeanvar


def append_colourbar(im, ax, side):
    divider = make_axes_locatable(ax)
    cax = divider.append_axes(side, size="5%", pad=0.05)
    plt.colorbar(im, cax=cax)
    cax.yaxis.set_ticks_position(side)


directory = sys.argv[1]
mean = np.loadtxt(f"{directory}/mean.txt")
std = np.loadtxt(f"{directory}/stddev.txt")
var = std ** 2
khist = np.loadtxt(f"{directory}/khistogram.txt")
best_fitting = np.loadtxt(f"{directory}/best_model.txt")
hpdrange = np.loadtxt(f"{directory}/hpdrange.txt")
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

        mean, var = meanvar_from_submeanvar(mean, mean2, var, var2, current_n, n_add)
        std = np.sqrt(var)
        current_n += n_add

        # what about the hpdrange?
        likelihoods = np.concatenate([likelihoods, likelihoods2])
        khist[:, 1] += np.loadtxt(f"{directory}/khistogram.txt", usecols=1)
        last_nonzero_k = np.argwhere(khist[:, 1]).max()
        best_fitting = np.loadtxt(f"{directory}/best_model.txt")
        khistory = np.concatenate([khistory, np.loadtxt(f"{directory}/khistory.txt")])

mosaic = """AB
            CD"""
fig = plt.figure(figsize=(10, 10))
axd = fig.subplot_mosaic(mosaic, gridspec_kw={"wspace": 0.05})

images = [mean, best_fitting, std, hpdrange - hpdrange.mean()]
sides = ["left", "right", "left", "right"]
titles = ["Mean", "Best fitting", "Standard Dev", "HPD range - <HPD range>"]
for ax, img, side, title in zip(axd, images, sides, titles):
    im = axd[ax].imshow(img, cmap="inferno")
    append_colourbar(im, axd[ax], side)
    axd[ax].axis("off")
    axd[ax].set_title(title)

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

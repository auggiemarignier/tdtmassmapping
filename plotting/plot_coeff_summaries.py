import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
import sys
import os
from utils import meanvar_from_submeanvar


def index_to_2d(index, width):
    ii = index % width
    ij = (index - ii) // width
    return ii, ij


def index_from_2d(ii, ij, width):
    return ij * width + ii


def build_mosaic_array(levels, cbar=False):
    assert levels >= 1
    mosaic = np.zeros((2 ** levels, 2 ** levels), dtype="<U2")
    mosaic[0, 0] = "0"
    for l in range(1, levels + 1):
        mosaic[: 2 ** (l - 1), 2 ** (l - 1) : 2 ** l] = f"{l}a"
        mosaic[2 ** (l - 1) : 2 ** l, : 2 ** (l - 1)] = f"{l}b"
        mosaic[2 ** (l - 1) : 2 ** l, 2 ** (l - 1) : 2 ** l] = f"{l}c"
    if cbar:
        c = np.full((2 ** levels, 1), "cbar", dtype="<U4")
        mosaic = np.hstack([mosaic, c])
    return mosaic


def img_to_mosaicaxes(img, axd, **kwargs):
    axd["0"].imshow(np.array([[img[0, 0]]]), **kwargs)
    for ax in axd:
        if ax != "cbar":
            if int(ax[0]) < 1:
                continue
            l = int(ax[:-1])
            if ax[-1] == "a":
                axd[ax].imshow(img[: 2 ** (l - 1), 2 ** (l - 1) : 2 ** l], **kwargs)
            elif ax[-1] == "b":
                axd[ax].imshow(img[2 ** (l - 1) : 2 ** l, : 2 ** (l - 1)], **kwargs)
            elif ax[-1] == "c":
                axd[ax].imshow(img[2 ** (l - 1) : 2 ** l, 2 ** (l - 1) : 2 ** l], **kwargs)
            else:
                raise KeyError(ax)
        else:
            if "vmin" in kwargs and "vmax" in kwargs:
                norm = Normalize(kwargs["vmin"], kwargs["vmax"])
            else:
                norm = Normalize()
            cmap = kwargs["cmap"] if "cmap" in kwargs else "viridis"
            plt.colorbar(ScalarMappable(norm, cmap), cax=axd[ax])
    return axd


directory = sys.argv[1]
try:
    levels = int(sys.argv[2])
except IndexError:
    levels = 8
means = np.loadtxt(f"{directory}/coeff_mean.txt")
stddevs = np.loadtxt(
    f"{directory}/coeff_std.txt"
).flatten()  # not sure why this shape is different
var = stddevs ** 2
counts = np.loadtxt(f"{directory}/coeff_n.txt")
current_n = np.loadtxt(f"{directory}/likelihood.txt").shape[0]
restarts = []
while os.path.isdir(f"{directory}/restart/"):
    restarts.append(current_n)
    directory += "/restart"
    if os.path.isfile(f"{directory}/coeff_mean.txt"):
        means2 = np.loadtxt(f"{directory}/coeff_mean.txt")
        stddevs2 = np.loadtxt(f"{directory}/coeff_std.txt").flatten()
        counts += np.loadtxt(f"{directory}/coeff_n.txt")
        var2 = stddevs2 ** 2
        n_add = np.loadtxt(f"{directory}/likelihood.txt").shape[0]

        mean, var = meanvar_from_submeanvar(means, means2, var, var2, current_n, n_add)
        current_n += n_add

for model, title in zip([means, stddevs, counts], ["Means", "Std Devs", "Counts"]):
    indexes = np.nonzero(model)[0]
    ii, ij = index_to_2d(indexes, 2 ** levels)
    wavelet_img = np.zeros((2 ** levels, 2 ** levels))
    wavelet_img[ii, ij] = model[indexes]

    mosaic = build_mosaic_array(levels, cbar=True)
    fig = plt.figure(constrained_layout=False, figsize=(10, 10))
    axd = fig.subplot_mosaic(mosaic, gridspec_kw={"wspace": 0, "hspace": 0})
    for ax in axd:
        if ax != "cbar":
            axd[ax].tick_params(
                axis="both",
                which="both",
                bottom=False,
                top=False,
                left=False,
                right=False,
                labelleft=False,
                labelbottom=False,
            )
    vmax = wavelet_img.max()
    vmin = wavelet_img.min()
    axd = img_to_mosaicaxes(
        wavelet_img,
        axd,
        cmap="cividis",
        vmin=vmin,
        vmax=vmax,
    )
    fig.suptitle(title, fontsize=24)
    plt.show()

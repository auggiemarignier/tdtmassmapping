import numpy as np
import matplotlib.pyplot as plt
import sys


def index_to_2d(index, width):
    ii = index % width
    ij = (index - ii) // width
    return ii, ij


def index_from_2d(ii, ij, width):
    return ij * width + ii


def build_mosaic_array(levels):
    assert levels >= 1
    mosaic = np.zeros((2 ** levels, 2 ** levels), dtype="<U2")
    mosaic[0, 0] = "0"
    for l in range(1, levels + 1):
        mosaic[: 2 ** (l - 1), 2 ** (l - 1) : 2 ** l] = f"{l}a"
        mosaic[2 ** (l - 1) : 2 ** l, : 2 ** (l - 1)] = f"{l}b"
        mosaic[2 ** (l - 1) : 2 ** l, 2 ** (l - 1) : 2 ** l] = f"{l}c"
    return mosaic


def img_to_mosaicaxes(img, axd, **kwargs):
    axd["0"].imshow(np.array([[img[0, 0]]]), **kwargs)
    for ax in axd:
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
    return axd


directory = sys.argv[1]
try:
    levels = int(sys.argv[2])
except IndexError:
    levels = 8
means = np.loadtxt(f"{directory}/coeff_mean.txt")
stddevs = np.loadtxt(f"{directory}/coeff_std.txt").flatten()
counts = np.loadtxt(f"{directory}/coeff_n.txt")

for model in [means, stddevs, counts]:
    indexes = np.nonzero(model)[0]
    ii, ij = index_to_2d(indexes, 256)
    wavelet_img = np.zeros((256, 256))
    wavelet_img[ii, ij] = model[indexes]

    mosaic = build_mosaic_array(levels)
    fig = plt.figure(constrained_layout=False, figsize=(10, 10))
    axd = fig.subplot_mosaic(mosaic, gridspec_kw={"wspace": 0, "hspace": 0})
    for ax in axd:
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
    cbar_end = np.abs(wavelet_img).max()
    axd = img_to_mosaicaxes(
        abs(wavelet_img), axd, cmap="cividis", vmin=0, vmax=cbar_end,
    )
    plt.show()

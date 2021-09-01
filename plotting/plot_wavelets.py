import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize


# Helper function used for visualization in the following examples
def identify_axes(ax_dict, fontsize=48):
    """
    Helper to identify the Axes in the examples below.

    Draws the label in a large font in the center of the Axes.

    Parameters
    ----------
    ax_dict : dict[str, Axes]
        Mapping between the title / label and the Axes.
    fontsize : int, optional
        How big the label should be.
    """
    kw = dict(ha="center", va="center", fontsize=fontsize, color="darkgrey")
    for k, ax in ax_dict.items():
        ax.text(0.5, 0.5, k, transform=ax.transAxes, **kw)


def index_to_2d(index, width):
    ii = index % width
    ij = (index - ii) // width
    return ii, ij


def index_from_2d(ii, ij, width):
    return ij * width + ii


def read_model(filename):
    with open(filename) as f:
        degree_width, degree_height = (int(i) for i in f.readline().split())
        width, height, size = (int(i) for i in f.readline().split())
        model = np.zeros(size)
        birth_set = np.zeros_like(model)
        next(f)
        max_depth = int(f.readline())
        Sv = True
        while Sv:
            depth, n_in_depth = (int(i) for i in f.readline().split())
            for _ in range(n_in_depth):
                st = f.readline().split()
                model[int(st[0])] = float(st[1])
            if depth == max_depth - 1:
                Sv = False
                Sb = True
                next(f)
        while Sb:
            depth, n_in_depth = (int(i) for i in f.readline().split())
            for _ in range(n_in_depth):
                st = int(f.readline())
                birth_set[st] = 1
            if depth == max_depth - 1:
                Sb = False
    return model, birth_set


def build_mosaic_array(levels):
    assert levels >= 1
    mosaic = np.zeros((2 ** levels, 2 ** levels), dtype="<U2")
    mosaic[0, 0] = "0"
    mosaic[0, 1] = "1a"
    mosaic[1, 0] = "1b"
    mosaic[1, 1] = "1c"
    for l in range(2, levels + 1):
        mosaic[: 2 ** l, 2 ** (l - 1) : 2 ** l] = f"{l}a"
        mosaic[2 ** (l - 1) : 2 ** l, : 2 ** l] = f"{l}b"
        mosaic[2 ** (l - 1) : 2 ** l, 2 ** (l - 1) : 2 ** l] = f"{l}c"
    return mosaic


model, _ = read_model("runs/checkerboard7/restart/final_model.txt")
mosaic = build_mosaic_array(8)
indexes = np.nonzero(model)[0]
ii, ij = index_to_2d(indexes, 256)
wavelet_img = np.zeros((256, 256))
wavelet_img[ii, ij] = model[indexes]
plt.imshow(wavelet_img)
plt.show()


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
identify_axes(axd, 30)
plt.show()

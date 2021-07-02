import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize

wavs_arr = np.loadtxt("../outputs/wavelets.txt")
print(wavs_arr.shape)

count_arr = (wavs_arr.astype(bool)).sum(axis=0).reshape(256, 256)
plt.imshow(count_arr, cmap="inferno")
plt.savefig("../outputs/wavelets.png")

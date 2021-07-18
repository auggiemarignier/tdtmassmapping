import networkx as nx
from networkx.drawing.nx_pydot import from_pydot
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import pydot
import numpy as np


def index_to_2d(index, width):
    ii = index % width
    ij = (index - ii) // width
    return ii, ij


def index_from_2d(ii, ij, width):
    return ij * width + ii


def parent_index(index, width):
    if index >= width ** 2 or index < 0:
        raise ValueError()
    if index == 0:
        return 0
    ii, ij = index_to_2d(index, width)
    ii = ii // 2
    ij = ij // 2
    return index_from_2d(ii, ij, width)


def depth_from_index(index, width):
    return 1 + depth_from_index(parent_index(index, width), width) if index != 0 else 0


def branch_string(index, width, depth):
    branch = np.zeros(depth + 1, dtype=int)
    branch[depth] = index
    current = index
    for dpt in range(depth - 1, -1, -1):
        parent = parent_index(current, width)
        branch[dpt] = parent
        current = parent
    return " -- ".join(branch.astype(str))


def tree_string(width, maxdepth):
    assert maxdepth <= np.log2(width)
    strings = []
    for index in range(width ** 2):
        depth = depth_from_index(index, width)
        if depth == maxdepth:
            strings.append(branch_string(index, width, maxdepth))
    return "\n".join(strings)


def generate_dot_string(width, depth):
    return f"graph my_graph {{{tree_string(width, depth)}}}"


img_width = 256
max_depth = 3
dot_string = generate_dot_string(img_width, max_depth)

usecols = []
for i in range(img_width ** 2):
    if depth_from_index(i, img_width) <= max_depth:
        usecols.append(i)
wavs_arr = np.loadtxt("outputs/wavelets.txt", usecols=usecols)

graphs = pydot.graph_from_dot_data(dot_string)
graph = graphs[0]
G = from_pydot(graph)
pos = nx.kamada_kawai_layout(G)


class UpdateGraph:
    def __init__(self, graph, pos, ax, colours_array, node_options={}, edge_options={}):
        self.graph = graph
        self.pos = pos
        self.ax = ax
        self.colours_array = colours_array  # (nframes, nnodes)
        self.node_options = node_options
        self.edge_options = edge_options
        nx.draw_networkx(
            self.graph,
            self.pos,
            ax=self.ax,
            node_color="w",
            edge_color="k",
            with_labels=False,
            **self.node_options,
            **self.edge_options,
        )

    def __call__(self, i):
        nx.draw_networkx_nodes(
            self.graph, self.pos, ax=self.ax, node_color=self.colours_array[i], **self.node_options
        )
        self.ax.set_title(f"Iteration {i}")


fig, ax = plt.subplots(figsize=(20, 8))
vmax = max([np.abs(wavs_arr.min()), wavs_arr.max()])
node_options = {"node_size": 100, "vmax": vmax, "vmin": -vmax, "cmap": "seismic"}
edge_options = {}

ud = UpdateGraph(G, pos, ax, wavs_arr, node_options, edge_options)
ani = FuncAnimation(fig, ud, frames=wavs_arr.shape[0], interval=1)
plt.show()

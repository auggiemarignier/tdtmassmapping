import networkx as nx
from networkx.drawing.nx_pydot import from_pydot
import matplotlib.pyplot as plt
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


def tree_string(width, indexes, birth_set):
    strings = []
    parents = [parent_index(child, width) for child in birth_set]
    for index in indexes:
        depth = depth_from_index(index, width)
        if index in parents:
            strings.append(branch_string(index, width, depth))
    return "\n".join(strings)


def generate_dot_string(width, indexes, birth_set):
    return f"graph my_graph {{{tree_string(width, indexes, birth_set)}}}"


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


model, birth_set = read_model("outputs/final_model.txt")

img_width = 256
max_depth = 7
dot_string = generate_dot_string(
    img_width, indexes=model.nonzero()[0], birth_set=birth_set.nonzero()[0]
)

graphs = pydot.graph_from_dot_data(dot_string)
graph = graphs[0]
G = from_pydot(graph)
pos = nx.kamada_kawai_layout(G)

fig, ax = plt.subplots(figsize=(20, 8))
node_options = {"node_size": 100, "cmap": "viridis"}
edge_options = {}
nx.draw_networkx(
    G,
    pos,
    ax=ax,
    **node_options,
    **edge_options
)
plt.show()

import networkx as nx
from networkx.drawing.nx_pydot import from_pydot
import matplotlib.pyplot as plt
import pydot


dot_string = """graph my_graph {
    a -- b1 -- c1 -- d1;
    a -- b1 -- c1 -- d2;
    a -- b1 -- c1 -- d3;
    a -- b1 -- c1 -- d4;
    a -- b1 -- c2 -- d5;
    a -- b1 -- c2 -- d6;
    a -- b1 -- c2 -- d7;
    a -- b1 -- c2 -- d8;
    a -- b1 -- c3 -- d9;
    a -- b1 -- c3 -- d10;
    a -- b1 -- c3 -- d11;
    a -- b1 -- c3 -- d12;
    a -- b1 -- c4 -- d13;
    a -- b1 -- c4 -- d14;
    a -- b1 -- c4 -- d15;
    a -- b1 -- c4 -- d16;
    a -- b2 -- c5 -- d17;
    a -- b2 -- c5 -- d18;
    a -- b2 -- c5 -- d19;
    a -- b2 -- c5 -- d20;
    a -- b2 -- c6 -- d21;
    a -- b2 -- c6 -- d22;
    a -- b2 -- c6 -- d23;
    a -- b2 -- c6 -- d24;
    a -- b2 -- c7 -- d25;
    a -- b2 -- c7 -- d26;
    a -- b2 -- c7 -- d27;
    a -- b2 -- c7 -- d28;
    a -- b2 -- c8 -- d29;
    a -- b2 -- c8 -- d30;
    a -- b2 -- c8 -- d31;
    a -- b2 -- c8 -- d32;
    a -- b3 -- c9 -- d33;
    a -- b3 -- c9 -- d34;
    a -- b3 -- c9 -- d35;
    a -- b3 -- c9 -- d36;
    a -- b3 -- c10 -- d37;
    a -- b3 -- c10 -- d38;
    a -- b3 -- c10 -- d39;
    a -- b3 -- c10 -- d40;
    a -- b3 -- c11 -- d41;
    a -- b3 -- c11 -- d42;
    a -- b3 -- c11 -- d43;
    a -- b3 -- c11 -- d44;
    a -- b3 -- c12 -- d45;
    a -- b3 -- c12 -- d46;
    a -- b3 -- c12 -- d47;
    a -- b3 -- c12 -- d48;
}
"""

graphs = pydot.graph_from_dot_data(dot_string)
graph = graphs[0]

plt.figure(figsize=(20, 8))
G = from_pydot(graph)
pos = nx.kamada_kawai_layout(G)

common_options = {"node_size": 100}
nx.draw(G, pos, **common_options)
nx.draw_networkx_nodes(G, pos, nodelist=["a"], node_color="red", **common_options)
nx.draw_networkx_nodes(
    G, pos, nodelist=["b1", "b2", "b3"], node_color="blue", **common_options
)
nx.draw_networkx_nodes(
    G, pos, nodelist=[f"c{i}" for i in range(1, 13)], node_color="green", **common_options
)
nx.draw_networkx_nodes(
    G, pos, nodelist=[f"d{i}" for i in range(1, 49)], node_color="yellow", **common_options
)
plt.show()

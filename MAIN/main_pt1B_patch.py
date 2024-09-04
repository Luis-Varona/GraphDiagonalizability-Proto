# %%
import os
import numpy as np
from sage.all import Graph, matrix

# %%
def laplacian_to_graph6(L: np.ndarray[int]) -> str:
    return Graph(matrix(np.diag(np.diag(L)) - L)).graph6_string()

# %%
def main():
    lap_int_con_graphs_order11 = np.load(
        "graph_data/LaplacianIntegralGraphs/LaplaciansIntegral_11.npy"
    ).swapaxes(0, 1).swapaxes(0, 2)

    lap_int_con_graphs_order11 = [
        laplacian_to_graph6(L) for L in lap_int_con_graphs_order11]
    np.save( # Save data to .npy file
        "graph_data/LaplacianIntegralGraphs/lap_int_con_graphs_order11.npy",
        lap_int_con_graphs_order11
    )

# %%
main()
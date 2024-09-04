# %%
import numpy as np
from itertools import chain
from sage.all import graphs
from typing import Generator

# %%
def is_laplacian_integral(L: np.ndarray[int]) -> bool:
    eigvals = np.linalg.eigvalsh(L)
    return np.allclose(eigvals, np.rint(np.real(eigvals)))

def connected_graphs(order: int) -> Generator:
    return graphs.nauty_geng(f"{order} -c -l")

def connected_regular_graphs(order: int) -> chain:
    def connected_k_regular_graphs(order: int, k: int) -> Generator:
        return graphs.nauty_geng(f"{order} -c -l -d{k} -D{k}")
    
    # An odd-order graph is k-regular only if k is even
    k_values = range(2, order, (order % 2) + 1)
    return chain.from_iterable(
        [connected_k_regular_graphs(order, k) for k in k_values])

def connected_bipartite_graphs(order: int) -> Generator:
    return graphs.nauty_geng(f"{order} -c -l -b")
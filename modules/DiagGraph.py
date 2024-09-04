# %%
import numpy as np
import networkx as nx
from sage.all import Graph as SageGraph
from typing import Union

# %%
class DiagGraph:
    def __init__(
            self, graph6_string: str,
            band_zerooneneg: int, band_oneneg: Union[float, int],
            eigvals: np.ndarray[int],
            eigvecs_zerooneneg: np.ndarray[int],
            eigvecs_oneneg: Union[None, np.ndarray[int]]) -> None:
        self.graph6_string = graph6_string
        self.band_zerooneneg = band_zerooneneg
        self.band_oneneg = band_oneneg
        self.eigvals = eigvals
        self.eigvecs_zerooneneg = eigvecs_zerooneneg
        self.eigvecs_oneneg = eigvecs_oneneg
    
    def __repr__(self) -> str:
        n = self.order()
        node_desc = "1 vertex" if n == 1 else f"{n} vertices"
        return f"{{-1,0,1}}-diagonalizable graph on {node_desc}"
    
    def order(self) -> int:
        return len(self.eigvals)
    def size(self) -> int:
        return int(np.sum(self.adjacency_matrix())/2)
    def density(self) -> float:
        n = self.order()
        return float('NaN') if n == 1 else 2*self.size()/(n*(n - 1))
    def average_degree(self) -> float:
        A = self.adjacency_matrix()
        return np.sum(A)/A.shape[0]
    
    def networkx_graph(self) -> nx.Graph:
        return nx.from_graph6_bytes(bytes(self.graph6_string, 'utf-8'))
    def sage_graph(self) -> SageGraph:
        return SageGraph(self.graph6_string)
    def adjacency_matrix(self) -> np.ndarray[int]:
        return self.sage_graph().adjacency_matrix().numpy()
    def laplacian_matrix(self) -> np.ndarray[int]:
        return self.sage_graph().kirchhoff_matrix().numpy()
    
    def is_regular(self) -> bool:
        return self.sage_graph().is_regular()
    def is_bipartite(self) -> bool:
        return self.sage_graph().is_bipartite()
    def is_cograph(self) -> bool:
        return self.sage_graph().is_cograph()
    def is_cart_prod(self) -> bool:
        return self.sage_graph().is_cartesian_product()
    def is_cart_prod_complement(self) -> bool:
        g_complement = self.sage_graph().complement()
        if not g_complement.is_connected():
            return False
        return g_complement.is_cartesian_product()
    def is_planar(self) -> bool:
        return self.sage_graph().is_planar()
    def is_outerplanar(self) -> bool:
        return self.sage_graph().is_circular_planar()
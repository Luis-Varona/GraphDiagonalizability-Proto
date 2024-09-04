# %%
import pulp as pl
import cvxpy as cp
import numpy as np
import networkx as nx

# %%
class WeightedLaplacian:
    def __init__(
            self, graph6_string: str, L: np.ndarray[float],
            weights: list[pl.LpVariable]) -> None:
        self.graph6_string = graph6_string
        self.L = L
        self.weights = weights

def adjacency_to_weighted_laplacian(
        A: np.ndarray[int], weights: list[pl.LpVariable]) -> np.ndarray[object]:
    A = A.astype(object)
    A[np.where(A)] = np.repeat(weights, 2)
    return np.diag(np.sum(A, axis = 1)) - A

def adjacency_to_weighted_laplacian_cvx(
        A: np.ndarray[int], weights: cp.Variable) -> np.ndarray[object]:
    A = A.astype(object)
    A[np.where(A)] = np.repeat([w for w in weights], 2)
    return np.diag(np.sum(A, axis = 1)) - A
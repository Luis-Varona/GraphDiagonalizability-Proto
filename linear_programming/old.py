# %%
import cvxpy as cp
import numpy as np
from itertools import chain

# %%
import sys
sys.path.append("../modules")

from weighted_laplacians import \
    WeightedLaplacian, adjacency_to_weighted_laplacian_cvx
from graph_generators import connected_graphs

# %%
weak_hadamards_orders1to4 = np.load(
    "weak_hadamards_orders1to4.npy", allow_pickle = True)

data = []

for n in range(3, 4):
    for i, g in enumerate(connected_graphs(n), 1):
        graph6 = g.graph6_string()
        A = g.adjacency_matrix().numpy()
        
        for j, W in enumerate(weak_hadamards_orders1to4[n - 1], 1):
            if j % 200 == 0:
                print(f"Graph {i}, WH {j}") # 21964
            weights = cp.Variable(g.size(), integer = True)
            L = adjacency_to_weighted_laplacian_cvx(A, weights)
            D = np.linalg.inv(W) @ L @ W
            R = D.copy()
            r = np.diagonal(R)
            r.setflags(write = True)
            r.fill(0)
            # R = np.diag(np.diag(D)) - D
            for i in range(n):
                for j in chain(range(i), range(i + 1, n)):
                    R[i, j] = cp.abs(R[i, j])
            
            off_diag = np.sum(R)
            # off_diag = cp.sum(list(chain.from_iterable([list(r) for r in R])))
            objective = cp.Minimize(off_diag)
            # objective = cp.Minimize(0)
            constraints = [weights >= 1]
            # constraints = [
            #     weights >= 1,
            #     off_diags[0] <= 0.1,
            #     off_diags[1] <= 0.1,
            #     off_diags[2] <= 0.1
            # ]
            problem = cp.Problem(objective, constraints)
            result = problem.solve()
            
            data.append((A, W, result, weights.value))
            
            if result < 0.01:
                print(A)
                print(W)
                print(result)
                print(weights.value)
                break

# Ignore everything after this point

# %%
# import pandas as pd

# p = np.array([[0, 0, 1], [0, 0, 1], [1, 1, 0]])
# k = np.array([[0, 1, 1], [1, 0, 1], [1, 1, 0]])
# one = list(map(
#     lambda x: {str(p): 'P3', str(k): 'K3'}[str(x)],
#     [x for (x, y, z) in data]))

# two = [y for (x, y, z) in data]

# three = [z for (x, y, z) in data]

# temp_data = pd.DataFrame(
#     {'Graph': one, 'Result': two, 'Weights': three})

# # %%
# P3_data = temp_data.loc[temp_data['Graph'] == 'P3']
# K3_data = temp_data.loc[temp_data['Graph'] == 'K3']

# print(list(K3_data['Result']).index(4.500000000228081))
# print(K3_data.iloc[27:28])
# print(K3_data.iloc[27]['Weights'])

# %%
# WH3 = weak_hadamards_orders1to4[2]
# unweighted_data = []
# for g in connected_graphs(3):
#     L = g.kirchhoff_matrix().numpy()
#     for W in WH3:
#         D = np.linalg.inv(W) @ L @ W
#         R = np.diag(np.diag(D)) - D
#         off_diag = np.sum([np.abs(r) for r in np.sum(R, 1)])
#         unweighted_data.append((L, off_diag))

# print(min([y for x, y in unweighted_data])) # Reversing inv yields same result

# %%
# def main():
#     weak_hadamards_orders1to5 = np.load("weak_hadamards_orders1to5.npy")
    
#     for n in range(3, 6):
#         for g in connected_graphs(n):
#             graph6 = g.graph6_string()
#             A = g.adjacency_matrix().numpy()
            
#             weights = cp.Variable(g.size())
#             L = adjacency_to_weighted_laplacian(A, weights)
#             for W in weak_hadamards_orders1to5[n]:
#                 D = np.linalg.inv(W) @ L @ W
                
#                 objective = cp.Minimize(
#                     cp.sum([cp.abs(row_sum) for row_sum in np.sum(D, 1)]))
#                 constraints = [weights > 0]
#                 prob = cp.Problem(objective, constraints)
                
#                 result = prob.solve()
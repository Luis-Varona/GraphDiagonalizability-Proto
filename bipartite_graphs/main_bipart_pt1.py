# %% Import packages
import os # For directory operations
import numpy as np # For array operations

# %% For graph generation and filtering to Laplacian integral graphs
import sys
sys.path.insert(1, "modules")
from graph_generators import is_laplacian_integral, connected_bipartite_graphs

# from modules.graph_generators \
#     import is_laplacian_integral, connected_bipartite_graphs

# %%
def main() -> None: # Define main program
    # Initialize array to store Laplacian integral bipartite graphs
    lap_int_con_bi_graphs_orders12to14 = np.empty(3, object)
    
    for order in range(12, 15): # Iterate over orders 12 to 14
        graph_list = [] # Initialize list to store graph6 strings
        
        # Iterate over all connected bipartite graphs
        for (i, g) in enumerate(connected_bipartite_graphs(order)):
        # for g in connected_bipartite_graphs(order):
            L = g.kirchhoff_matrix().numpy() # Convert to Laplacian matrix
            if is_laplacian_integral(L): # Test for integer eigenvalues
                graph_list.append(g.graph6_string()) # Store graph6 string
            
            if order == 12:
                if i % 10**4 == 0:
                    print(f"Order 12: {i/10**3} out of ~212.7 thousand done")
            elif order == 13:
                if i % 10**5 == 0:
                    print(f"Order 13: {i/10**6} out of ~2.2 million done")
            elif order == 14:
                if i % 10**5 == 0:
                    print(f"Order 14: {i/10**6} out of ~31.2 million done")
        
        # Convert native list of graph6 strings to numpy array
        lap_int_con_bi_graphs_orders12to14[order - 12] = np.array(graph_list)
    
    # Create subdirectory for saved data if it does not already exist
    if not os.path.isdir("graph_data/LaplacianIntegralGraphs/bipartite"):
        os.makedirs("graph_data/LaplacianIntegralGraphs/bipartite")
    np.save( # Save data to .npy file
        "graph_data/LaplacianIntegralGraphs/bipartite/"
        "lap_int_con_bi_graphs_orders12to14.npy",
        lap_int_con_bi_graphs_orders12to14
    )

# %%
main() # Run main program
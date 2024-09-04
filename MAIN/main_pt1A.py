"""
If you wish to error-check this program, it should take 0.5 to 2 hours to run,
depending on the speed of your computer.
"""

# %% Import packages
import os # For directory operations
import numpy as np # For array operations

# For graph generation and filtering to Laplacian integral graphs
from modules.graph_generators \
    import is_laplacian_integral, connected_graphs

# %% Save all Laplacian integral graphs of order <= 10 in graph6 format
def main() -> None: # Define main program
    # Initialize array to store Laplacian integral graphs
    lap_int_con_graphs_orders1to10 = np.empty(10, object)
    
    for order in range(1, 11): # Iterate over orders 1 to 10
        graph_list = [] # Initialize list to store graph6 strings
        
        for g in connected_graphs(order): # Iterate over all connected graphs
            L = g.kirchhoff_matrix().numpy() # Convert to Laplacian matrix
            if is_laplacian_integral(L): # Test for integral eigenvalues
                graph_list.append(g.graph6_string()) # Store graph6 string
        
        # Convert native list of graph6 strings to numpy array
        lap_int_con_graphs_orders1to10[order - 1] = np.array(graph_list)
    
    # Create subdirectory for saved data if it does not already exist
    if not os.path.isdir("graph_data/LaplacianIntegralGraphs"):
        os.makedirs("graph_data/LaplacianIntegralGraphs")
    np.save( # Save data to .npy file
        "graph_data/LaplacianIntegralGraphs/" \
            "lap_int_con_graphs_orders1to10.npy",
        lap_int_con_graphs_orders1to10
    )

main() # Run main program
"""
If you wish to error-check this program, each of the 13 'main()' calls should
take 4 to 18 hours to run, depending on the speed of your computer. (Later calls
will take longer to run due to the nature of itertools' 'islice' function.)
"""

# %% Import packages
import os # For directory operations
import numpy as np # For array operations
from itertools import islice #For taking slices of graph generators

# For graph generation and filtering to Laplacian integral graphs
from modules.graph_generators \
    import is_laplacian_integral, connected_graphs

# %% Save all Laplacian integral graphs of order <= 10 in graph6 format
# 13 slices each iterating over 77,438,505 graphs
def get_generator_slice(num_bundle: int) -> None:
    item_start = 77438505*(num_bundle - 1)
    gen = islice(connected_graphs(11), item_start, item_start + 77438505)
    return gen

def main(num_bundle: int) -> None: # Define main program
    graph_list = [] # Initialize list to store graph6 strings
    gen = get_generator_slice(num_bundle)
    
    for i, g in enumerate(gen): # Iterate over connected graphs of order 11
        L = g.kirchhoff_matrix().numpy() # Convert to Laplacian matrix
        if is_laplacian_integral(L): # Test for integral eigenvalues
            graph_list.append(g.graph6_string()) # Store graph6 string
        if i % 10**6 == 0:
            print(f"{i/10**6} out of 77.438505 million done")
    
    # Convert native list of graph6 strings to numpy array
    lap_int_con_graphs_order11_bundle = np.array(graph_list)
    
    # Create subdirectory for saved data if it does not already exist
    if not os.path.isdir("graph_data/LaplacianIntegralGraphs/order11bundles"):
        os.makedirs("graph_data/LaplacianIntegralGraphs/order11bundles")
    np.save( # Save data to .npy file
        "graph_data/LaplacianIntegralGraphs/order11bundles/"
        f"lap_int_con_graphs_order11_bundle{num_bundle:02d}.npy",
        lap_int_con_graphs_order11_bundle
    )

def main_final() -> None:
    lap_int_con_graphs_order11 = np.concatenate([
        np.load(
            f"graph_data/LaplacianIntegralGraphs/order11bundles/"
            f"lap_int_con_graphs_order11_bundle{i:02d}.npy"
        ) for i in range(1, 14)
    ])
    np.save(
        "graph_data/LaplacianIntegralGraphs/lap_int_con_graphs_order11.npy",
        lap_int_con_graphs_order11
    )

# %% Run main program (uncomment one of the following lines to run)
# main(1)
# main(2)
# main(3)
# main(4)
# main(5)
# main(6)
# main(7)
# main(8)
# main(9)
# main(10)
# main(11)
# main(12)
# main(13)
# main_final()
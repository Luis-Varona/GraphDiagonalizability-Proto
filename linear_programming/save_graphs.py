# %%
import os
import numpy as np
from scipy.io import savemat
from sage.all import graphs

def main():
    # %%
    graph6_strings_orders3to6 = np.empty(4, object)
    
    for n in range(4):
        graph6_strings_orders3to6[n] = np.array(
            [g.graph6_string() for g in graphs.nauty_geng(f"{n + 3} -c -l")])
    
    if not os.path.isdir("linear_programming/data"):
        os.makedirs("linear_programming/data")
    np.save(
        "linear_programming/data/graph6_strings_orders3to6.npy",
        graph6_strings_orders3to6)

# %%
main()
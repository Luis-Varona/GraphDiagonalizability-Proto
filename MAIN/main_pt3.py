# %% Import modules
import numpy as np
import sqlite3 as sql
from io import BytesIO

# %%
import sys
sys.path.insert(1, "modules")
from DiagGraph import DiagGraph

# %% Interoperability between Julia, Python, and SQL
# To port the 'DiagGraph' class from Julia to Python
def PyPort_to_DiagGraphs(arr: np.ndarray[object]) -> np.ndarray[DiagGraph]:
    return np.array(
        [DiagGraph(
            graph6_string = obj['graph6_string'],
            band_zerooneneg = int(obj['band_zerooneneg']),
            band_oneneg = obj['band_oneneg'],
            eigvals = obj['eigvals'],
            eigvecs_zerooneneg = obj['eigvecs_zerooneneg'],
            eigvecs_oneneg = obj['eigvecs_oneneg'])
        for obj in arr])

# To store 'DiagGraph' class properties in a tuple for SQL mapping
def DiagGraph_to_SQL(g: DiagGraph) -> tuple:
    return (
        g.order(),
        g.graph6_string,
        g.adjacency_matrix(),
        g.laplacian_matrix(),
        
        g.band_zerooneneg,
        g.band_oneneg,
        
        g.eigvals,
        g.eigvecs_zerooneneg,
        g.eigvecs_oneneg,
        
        g.is_regular(),
        g.is_bipartite(),
        g.is_cograph(),
        g.is_cart_prod(),
        g.is_cart_prod_complement(),
        g.is_planar(),
        g.is_outerplanar(),
        
        g.size(),
        g.density(),
        g.average_degree()
    )

# %% Define register methods for custom SQL types
def adapt_array(array: np.ndarray):
    out = BytesIO()
    np.save(out, array)
    out.seek(0)
    return sql.Binary(out.read())

def convert_array(obj):
    out = BytesIO(obj)
    out.seek(0)
    return np.load(out, allow_pickle = True)

# %% Declare custom SQL types and initialize a database connection
def main():
    sql.register_adapter(np.ndarray, adapt_array)
    sql.register_converter("ARRAY", convert_array)
    
    sql.register_adapter(bool, str)
    sql.register_converter("BOOLEAN", eval)
    
    conn = sql.connect(
        "graph_data/DiagonalizableGraphs/Con_ZeroOneNegDiag_Graphs.db",
        detect_types = sql.PARSE_DECLTYPES)
    cur = conn.cursor()
    
# %% Load data from Julia and convert to instances of 'DiagGraph'
    PyPort_diag_con_graphs_orders1to11 = np.load(
        "graph_data/DiagonalizableGraphs/diag_con_graphs_orders1to11.npy",
        allow_pickle = True)
    diag_con_graphs_orders1to11 = [
        PyPort_to_DiagGraphs(arr) for arr
        in PyPort_diag_con_graphs_orders1to11]
    
# %% Map all graph objects to the SQL database and commit
    for (order, graphs) in enumerate(diag_con_graphs_orders1to11, 1):
        table = f"Order{order:02d}"
        cur.execute(
            f"CREATE TABLE {table} ("
            "GraphOrder TINYINT, "
            f"graph6 VARCHAR({order}), "
            "Adjacency ARRAY, "
            "Laplacian ARRAY, "
            
            "Band01Neg TINYINT(2), "
            "Band1Neg TINYINT(2), "
            
            "Eigvals ARRAY, "
            "Eigvecs01Neg ARRAY, "
            "Eigvecs1Neg ARRAY, "
            
            "Regular BOOLEAN, "
            "Bipartite BOOLEAN, "
            "Cograph BOOLEAN, "
            "CartProd BOOLEAN, "
            "CartProdComp BOOLEAN, "
            "Planar BOOLEAN, "
            "Outerplanar BOOLEAN, "
            
            "Size TINYINT(2), "
            "Density REAL, "
            "AvgDegree REAL)"
        )
        for g in graphs:
            graph_data = DiagGraph_to_SQL(g)
            cur.execute(
                f"INSERT INTO {table} VALUES ({('?,'*19)[:-1]})", graph_data
            )
    
    conn.commit()
    conn.close()

# %% Run main program
main()
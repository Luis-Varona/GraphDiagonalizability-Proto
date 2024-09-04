# %% Import modules
import numpy as np
import pandas as pd
import sqlite3 as sql
from io import BytesIO

# %%
from sys import path
path.append("..")

from modules.DiagGraph import DiagGraph
from pusheen import filter_to

# %% Set pandas DataFrame display options
pd.set_option('display.min_rows', 20)
pd.set_option('display.max_columns', None)

# %% To convert SQL data to 'DiagGraph' instances
def SQL_to_DiagGraph(obj):
    return DiagGraph(obj[1], obj[4], obj[5], obj[6], obj[7], obj[8])

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
sql.register_adapter(np.ndarray, adapt_array)
sql.register_converter("ARRAY", convert_array)

sql.register_adapter(bool, str)
sql.register_converter("BOOLEAN", eval)

conn = sql.connect(
    "../graph_data/DiagonalizableGraphs/Con_ZeroOneNegDiag_Graphs.db",
    detect_types = sql.PARSE_DECLTYPES)
cur = conn.cursor()

# %% Load and convert graph data to numpy format from the SQL database
diag_con_graphs_orders1to10 = np.empty(10, object)

for order in range(1, 11):
    table = f'Order{order:02d}'
    cur.execute(f"SELECT * from {table}")
    data = cur.fetchall()
    diag_con_graphs_orders1to10[order - 1] = np.array(
        [SQL_to_DiagGraph(obj) for obj in data])

# %% Load and convert graph data to pandas format from the SQL database
df_list = [0]*11

for order in range(1, 12):
    table = f'Order{order:02d}'
    df_list[order - 1] = pd.read_sql_query(f"SELECT * from {table}", conn)

conn.close()

graph_data = pd.concat(
    [df.drop(df.columns[[2, 3, 6, 7, 8]], axis = 1) for df in df_list])
graph_data.reset_index(drop = True, inplace = True)

# %% Garbage collection
del(cur, data, order, table)

# %%
special_cases = filter_to(
    graph_data, (False, False, False), cols = graph_data.columns[6:9])
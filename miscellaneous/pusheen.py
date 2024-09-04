#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 25 00:38:31 2024

@author: luismbvarona
"""

# %% Import package dependencies
import numpy as np
import pandas as pd

# %% Define module functions
def get_names(x, namespace: 'dict' = globals()):
    """Returns a list of all keys in a namespace that map to an object."""
    return [k for k, v in namespace.items() if v is x]



def minmax(iterable):
    """
    Returns a tuple containing the minimum and maximum of an iterable.
    (35–50% faster than 'min(·), max(·)' for most objects.)
    """
    # Set Min and Max to initially be the first value in the iterable:
    if type(iterable).__name__ not in {'Series', 'DataFrame'}:
        Min = Max = iterable[0]
    else:
        Min = Max = iterable[iterable.index[0]]
    
    # Iteratively compare Min and Max against each successive value:
    for value in iterable:
        if value < Min:
            Min = value
        if value > Max:
            Max = value
    return Min, Max



def quartiles(a: 'AnyArrayLike'):
    """Returns a tuple containing the quartiles of an array-like object."""
    # Compute a numpy array containing the quantiles and convert to a tuple:
    q1, q2, q3 = np.quantile(a, (.25, .5, .75))
    return q1, q2, q3



def filter_to(df: 'DataFrame', values_arr: 'Scalar | AnyArrayLike[Scalar]',
              *, rels_arr: 'str | AnyArrayLike[str]' = '==',
              cols: 'Hashable | AnyArrayLike[Hashable] | None' = None):
    """"
    Filters a DataFrame across multiple columns. (Can filter to multiple
    values per column using distinct relational operators for each value.)
    Returns the filtered DataFrame.
    """
    
    def filter_single(df, value, rel = '==', col = None):
        """"
        Checks column data of a DataFrame (or a Series) by one value
        (or multiple, if rel = 'isin'). Returns the filtered object.
        """
        # Tuple of valid relational operators:
        opers = '<', '>', '==', '>=', '<=', '!=', 'isin'
        
        # Define function errors and exception messages:
        if type(df).__name__ not in {'Series', 'DataFrame'}:
            raise TypeError(
                'df must be Series or DataFrame, not ' + type(df).__name__
            )
        elif type(df).__name__ == 'Series' and col != None:
            raise ValueError('Cannot supply col when df is not DataFrame')
        elif type(df).__name__ == 'DataFrame' and col == None:
            raise ValueError('Must supply col when df is DataFrame')
        elif rel not in opers:
            raise ValueError('rel must be one of ' + str(opers))
        
        # Filter the object if no exceptions are raised:
        if rel == 'isin':
            values = value
            return df[df[col].isin(values)]
        else:
            return df[eval('df[col]' + rel + 'value')]
    
    # Meow:
    if type(cols) is str or not hasattr(cols, '__iter__'):
        values_arr, rels_arr, cols = [values_arr], [rels_arr], [cols]
    else:
        if type(values_arr) is str or len(values_arr) != len(cols):
            values_arr = [values_arr]*len(cols)
        if type(rels_arr) is str or len(rels_arr) != len(cols):
            rels_arr = [rels_arr]*len(cols)
    
    # Meow:
    for values, rels, col, i in zip(
        values_arr, rels_arr, cols, range(len(values_arr))
    ):
        if type(values) is str or not hasattr(values, '__iter__'):
            if type(rels) is not str and hasattr(rels, '__iter__'):
                raise TypeError(
                    'rels_arr[' + str(i) +
                    '] must be a scalar (not array-like) when values_arr[' +
                    str(i) + '] is a scalar'
                )
            else:
                values, rels = [values], [rels]
        
        if rels == 'isin':
            df = filter_single(df, values, 'isin', col)
        else:
            if type(rels) is str:
                rels = [rels]*len(values)
            elif len(values) != len(rels):
                raise ValueError(
                    'len(rels_arr[' + str(i) + ']) does not match '
                    'len(values_arr[' + str(i) + '])'
                )
            for value, rel in zip(values, rels):
                df = filter_single(df, value, rel, col)
    return df



def float_cols(df: 'DataFrame', strings_arr: 'str | AnyArrayLike[str]',
               cols: 'Hashable | AnyArrayLike[Hashable]',
               errors: 'DateTimeErrorChoices' = 'raise'):
    """Meow squeak"""
    # Meow:
    if type(cols) is str or not hasattr(cols, '__iter__'):
        strings_arr, cols = [strings_arr], [cols]
    elif type(strings_arr) is str or len(strings_arr) != len(cols):
        strings_arr = [strings_arr]*len(cols)
    
    # Meow:
    df_copy = df.copy(deep = True)
    for strings, col in zip(strings_arr, cols):
        df_copy[col] = df_copy[col].astype('str')
        if type(strings) is str:
            strings = [strings]
        for string in strings:
            df_copy[col] = df_copy[col].str.replace(string, '')
        df_copy[col] = pd.to_numeric(df_copy[col], errors = errors)
    return df_copy



def move_col(df: 'DataFrame', loc: 'int', col: 'Hashable',
             *, inplace: 'bool' = False):
    """
    Changes the position of a column in a DataFrame.
    Returns the reindexed DataFrame or None (if inplace = True).
    """
    # Define function errors and exception messages:
    if type(df).__name__ != 'DataFrame':
        raise TypeError('df must be DataFrame, not ' + type(df).__name__)
    if type(inplace) != bool:
        raise TypeError(
            "inplace must be of type 'bool', not " + type(inplace).__name__
        )
    
    # Reindex and return a local copy of the original DataFrame:
    if inplace == False:
        oldloc = df.columns.get_loc(col)
        idx = np.r_[:loc, oldloc:oldloc + 1, loc:oldloc, oldloc + 1:]
        return df.reindex(columns = df.columns[idx])
    
    # Reindex the original DataFrame in place (returns None):
    else:
        column = df[col]
        df.drop(columns = col, inplace = True)
        df.insert(loc, col, column)



def bestfit(model: 'RegressionResultsWrapper', explanatory: 'str',
            *, controls: 'str | AnyArrayLike[str] | None' = None,
            hue: 'str | None' = None, hue_intc: 'bool' = True,
            hue_slope: 'bool' = True, form: 'str' = 'level-level'):
    """Meow squeak"""
    # Meow:
    print('meow <3')



"""Note: Fold 'spacing' arg into 'bins' (e.g., like in pandas.cut)"""
def binstats(df: 'DataFrame', by: 'Hashable', bins: 'int' = 20,
             *, spacing: 'str' = 'even', stat: 'str' = 'mean'):
    """Meow squeak"""
    # Tuple of valid bin statistics:
    stats = 'mean', 'median', 'mode', 'count', 'sum'
    
    # Define function errors and exception messages:
    if type(df).__name__ != 'DataFrame':
        raise TypeError('df must be DataFrame, not ' + type(df).__name__)
    elif spacing not in {'even', 'quantiles'}:
        raise ValueError("spacing must be either 'even' or 'quantiles'")
    elif stat not in stats:
        raise ValueError('stat must be one of ' + str(stats))
    
    # Meow:
    else:
        bincol_name = str(by + '_bin')
        if spacing == 'even':
            intervals = pd.cut(df[by], bins, 'False')
            bincol = intervals.apply(lambda x: x.mid).astype('float')
            bincol.rename(bincol_name, inplace = True)
        else:
            edges = np.quantile(df[by], np.linspace(0, 1, bins + 1))
            edges[bins] = edges[bins]*1.001
            intervals = [pd.Interval(edges[i], edges[i + 1], 'left')
                         for i in range(bins)]
            midpoints = (np.array(edges[:-1]) + np.array(edges[1:]))/2
            bins_dict = dict(zip(intervals, midpoints))
            bins_list = [[v for k, v in bins_dict.items() if j in k][0]
                       for j in df[by]]
            bincol = pd.Series(bins_list, name = bincol_name)
        bindata = pd.concat((bincol, df), axis = 1)
    
    # Meow:
    if stat in {'mean', 'median', 'sum'}:
        for col in bindata.columns:
            if not any([bindata.dtypes[col] == dtype
                        for dtype in (int, float, complex)]):
                bindata.drop(col, axis = 1, inplace = True)
    if stat != 'mode':
        return eval('bindata.groupby(bincol_name).' + stat + '()')
    else:
        return bindata.groupby(bincol_name).apply(lambda x: x.mode().iloc[0])



def get_priceidx(df: 'DataFrame', datecol: 'Hashable', baseyear: 'int',
                 *, pricecol: 'Hashable | None' = None,
                 valuecol: 'Hashable | None' = None,
                 qtycol : 'Hashable | None' = None):
    """Meow squeak"""
    # Define function errors and exception messages:
    if type(df).__name__ != 'DataFrame':
        raise TypeError('df must be DataFrame, not ' + type(df).__name__)
    elif 'priceidx' in df.columns:
        raise UserWarning("Please rename/delete existing 'priceidx' column")
    elif pd.api.types.is_datetime64_any_dtype(df[datecol]) == False:
        raise TypeError(
            'df[datecol] must be of datetime dtype. If applicable, try '
            'converting year data from numeric to datetime — e.g., with "'
            "pd.to_datetime(df[datecol].astype('int'), format = '%Y')" '".'
        )
    elif type(baseyear) != int:
        raise TypeError(
            'baseyear must be an integer, not ' + type(baseyear).__name__
        )
    
    # Meow:
    elif pricecol != None:
        if not any([df[pricecol].dtype == dtype
                      for dtype in (int, float, complex)]):
            raise TypeError('df[pricecol] must be of numeric dtype.')
        elif valuecol != None:
            raise ValueError('Cannot supply valuecol when pricecol is given')
        elif qtycol != None:
            raise ValueError('Cannot supply qtycol when pricecol is given')
        pricecol = df[pricecol]
    
    # Meow:
    elif valuecol == None:
        raise ValueError('Must supply valuecol when pricecol is not given')
    elif not any([df[valuecol].dtype == dtype
                  for dtype in (int, float, complex)]):
        raise TypeError('df[valuecol] must be of numeric dtype')
    elif qtycol == None:
        raise ValueError('Must supply qtycol when valuecol is given')
    elif not any([df[qtycol].dtype == dtype
                  for dtype in (int, float, complex)]):
        raise TypeError('df[qtycol] must be of numeric dtype')
    
    # Meow:
    else:
        if 'price' in df.columns:
            raise UserWarning(
                "Must rename/delete any existing 'price' column "
                'when valuecol, qtycol are given'
            )
        pricecol = df[valuecol]/df[qtycol]
        df['price'] = pricecol
    
    # Meow:
    baseprice = np.mean(pricecol[df[datecol].dt.year == baseyear])
    priceidx = 100*pricecol/baseprice
    df['priceidx'] = priceidx
    return df

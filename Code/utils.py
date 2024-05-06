import requests
import openpyxl
import pandas as pd
import numpy as np
from io import BytesIO
from cnmf import cNMF
import scanpy as sc
from scipy.stats import ttest_rel
from sklearn.neighbors import NearestNeighbors


############################## Loading data from Google Drive #############################

def read_gsheet(sheet_id, tabname):
    url = "https://docs.google.com/spreadsheets/export?exportFormat=xlsx&id=" + sheet_id
    res = requests.get(url)
    data = BytesIO(res.content)
    xlsx = openpyxl.load_workbook(filename=data)
    return(pd.read_excel(data, sheet_name=tabname))


def read_dataset_log(tabname='Current Dataset Paths'):
    return(read_gsheet('1gJMdNBd7qkn_peRQEuz2ZHQaCO2EBMPCJnJZrjjzkBU', tabname))



############################## Plotting #############################

def compute_smooth_scatter_color(x, y, z, n_neighbors=5):
    '''Compute local averages of z averaging over KNN to make smoothed scatter plot'''
    nbrs = NearestNeighbors(n_neighbors=5, algorithm='auto').fit(np.column_stack([x, y]))
    distances, indices = nbrs.kneighbors(np.column_stack([x, y]))
    z_avg = np.mean(z[indices], axis=1)
    return z_avg


############################## Matching GEPs between datasets #############################
def df_col_corr(X, Y):
    '''
    Compute pairwise Pearson correlation matrix of columns of X and Y. Returns
    R which is n_X_cols X n_Y_cols
    '''
    X_norm = X.subtract(X.mean(axis=0), axis=1)
    X_norm= X_norm.divide(X_norm.std(axis=0), axis=1)
    Y_norm = Y.subtract(Y.mean(axis=0), axis=1)
    Y_norm= Y_norm.divide(Y_norm.std(axis=0), axis=1)
    
    X_mask = np.ma.array(X_norm.values, mask=np.isnan(X_norm.values))
    Y_mask = np.ma.array(Y_norm.values, mask=np.isnan(Y_norm.values))
    R = np.ma.dot(X_mask.T, Y_mask) / (X.shape[0]-1)
    
    R = pd.DataFrame(R, index=X.columns, columns=Y.columns)
    return(R)


def match_columns(X, Y):
    '''
    Assign each column in Y to a distinct column in X via greedy search of max Pearson correlation.
    The # columns in X must be >= the # columns in Y
    '''

    if X.shape[1] < Y.shape[1]:
        raise Exception('X must have more columns than Y')
    
    R = df_col_corr(X, Y)
    mapping = R.unstack().reset_index().sort_values(by=0, ascending=False)
    mapping.columns = ['Y_columns', 'X_columns', 'R']

    used_X = []
    used_Y = []
    todrop = []
    for i in mapping.index:
        if mapping.at[i, 'X_columns'] in used_X:
            todrop.append(i)
        elif mapping.at[i, 'Y_columns'] in used_Y:
            todrop.append(i)
        else:
            used_X.append(mapping.at[i, 'X_columns'])
            used_Y.append(mapping.at[i, 'Y_columns'])  
            
    mapping = mapping.drop(todrop)
    unmatched = list(set(X.columns) - set(used_X))
    return(mapping, unmatched, R)


############################## Hypothesis testing #############################

def ttest_paired_allcols(X, Y):
    '''
    Takes 2 pandas.DataFrames X and Y with the same number of rows and columns.
    Performs a paired t-test for each column of X and Y. 
    
    Returns
    ---------------------------------
    Ts - pandas.Series of T statistics
    Ps - pandas.Series of P-values
    '''
    
    Ts = []
    Ps = []
    for g in X.columns:
        T, P = ttest_rel(X[g], Y[g])
        Ts.append(T)
        Ps.append(P)

    Ts = pd.Series(Ts, index=X.columns)
    Ps = pd.Series(Ps, index=X.columns)
    return(Ts, Ps)


def permute_within_group(df, group_col, vartopermute):
    '''
    Performs permutation of the column vartopermute in the pandas.DataFrame df 
    but only permutes within the column group_col (which is assumed to contain discrete
    categories)

    Returns:
    permuted - a copy of df containing the additional column _permuted which can be used for
    downstream hypothesis testing
    '''
    permuted = df.groupby(group_col, group_keys=True).apply(lambda x: x.assign(_permuted=np.random.permutation(x[vartopermute]))).reset_index(drop=True)
    return(permuted)
    

############################## TCAT ##########################################


class starcat(cNMF):
    def __init__(self, alpha=0.0, l1_ratio=0.0,  tpm_norm=False, var_norm=True, copy=True):
        """
        Class for running TCAT on a query dataset. By default uses a standard reference
        database for human T-cells but you can optionally provide your own reference 
        object.
        """        
        self._nmf_kwargs = dict(
                        alpha_W=alpha,
                        alpha_H=0.0,
                        l1_ratio=l1_ratio,
                        beta_loss='frobenius',
                        solver='mu',
                        tol=1e-4,
                        max_iter=1000,
                        init='random',
                        update_H = False
                        )
        
        self.tpm_norm=tpm_norm
        self.var_norm=var_norm
        self.copy=copy
        
    def fit_transform(self, query, ref_spectra=None):
        """
        Takes an input data matrix and a fixed spectra and uses NNLS to find the optimal
        usage matrix. Generic kwargs for NMF are loaded from self.paths['nmf_run_parameters'].
        If input data are pandas.DataFrame, returns a DataFrame with row index matching X and
        columns index matching index of spectra

        Parameters
        ----------
        X : pandas.DataFrame or numpy.ndarray, cells X genes
            Non-negative expression data to fit spectra to

        spectra : pandas.DataFrame or numpy.ndarray, programs X genes
            Non-negative spectra of expression programs

        adata : query AnnData object
        ref_spectra : pd.DataFrame (#GEPs X #Genes) or None, optional (default=None, uses standard reference)
        """

        if ref_spectra is None:
            # Use default hard-coded reference (might live online or in the
            # TCAT code directory)
            print('Using default TCAT reference set')
            raise NotImplementedError("Not yet implemented.")
        else:
            self.ref = ref_spectra            
        self.ref = self.ref.astype(np.float32)
        
        num_nulls = self.ref.isnull().sum(axis=0)
        has_null = num_nulls>0
        nnull_genes = has_null.sum()
        if nnull_genes > 0:
            print('%d out of %d reference contain NaN in at least one GEP. Filtering them.' % (nnull_genes, self.ref.shape[1]))
            self.ref = self.ref.loc[:,~has_null]
            
        
        self.prepare_query(query)
        rf_usages = self.fit_query_usage()
        return(rf_usages)
        
    
    def prepare_query(self, adata):
        """
        Load query dataset as AnnData object and optionally perform normalization.

        adata : query AnnData object
        tp10k_norm : boolean, optional (default=True, performs TP10K normalization)
        var_norm : boolean, optional (default=True, performs gene variance normalization)
        copy : boolean, optional (default=True, stores a copy of adata rather than modifying in place)

        """

        if self.copy:
            self.query = adata.copy()
        else:
            self.query = adata
            
        overlap_genes = list(set(self.ref.columns).intersection(set(adata.var.index)))
        print('%d out of %d genes in the reference overlap with the query' % (len(overlap_genes), self.ref.shape[1]))
        self.overlap_genes = overlap_genes
        self.query = self.query[:, self.overlap_genes]
            
        if self.tpm_norm:
            sc.pp.normalize_per_cell(self.query, counts_per_cell_after=1e6)
            
        if self.var_norm:
            sc.pp.scale(self.query, zero_center=False)
         
    def fit_query_usage(self):
        rf_usages = self.refit_usage(self.query.X, self.ref[self.overlap_genes].values,
                         self._nmf_kwargs.copy())          
        rf_usages = pd.DataFrame(rf_usages, index=self.query.obs.index,
                                 columns=self.ref.index)
        return(rf_usages)
        
    def refit_usage(self, X, spectra, nmf_kwargs):
        """
        Takes an input data matrix and a fixed spectra and uses NNLS to find the optimal
        usage matrix. Generic kwargs for NMF are loaded from self.paths['nmf_run_parameters'].
        If input data are pandas.DataFrame, returns a DataFrame with row index matching X and
        columns index matching index of spectra

        Parameters
        ----------
        X : pandas.DataFrame or numpy.ndarray, cells X genes
            Non-negative expression data to fit spectra to

        spectra : pandas.DataFrame or numpy.ndarray, programs X genes
            Non-negative spectra of expression programs
        """

        nmf_kwargs['H'] = spectra
        nmf_kwargs['n_components'] = spectra.shape[0]
        _, rf_usages = self._nmf(X, nmf_kwargs=nmf_kwargs)            
        return(rf_usages)
    
import pandas as pd
import numpy as np
from scsim2 import scsim2

def isarray(x):
    return(hasattr(x, "__len__"))


def rescale(x, new_min, new_max):
    z = x - x.min()
    z = z / z.max()
    z = z*(new_max-new_min)+new_min
    return(z)

class SCsim_DiffScale(scsim2):
    
    def __init__(self, seed=757578, mean_rate=7.68, mean_shape=0.34, libloc=7.64, libscale=0.78,
                expoutprob=0.00286, expoutloc=6.15, expoutscale=0.49, bcv_dispersion=0.448, bcv_dof=22.087):
        '''
        Create an scsim2 object to simulate single-cell RNA-seq dataset.
        
        Default Single-cell parameters were estimated from an example dataset as described in the cNMF paper
        '''
        self.seed = seed
        self.mean_rate = mean_rate
        self.mean_shape = mean_shape
        self.libloc = libloc
        self.libscale = libscale
        self.expoutprob = expoutprob
        self.expoutloc = expoutloc
        self.expoutscale = expoutscale
        self.bcv_dispersion = bcv_dispersion
        self.bcv_dof = bcv_dof

    def simulate(self, ncells=10000, ngenes=10000, nidentities=8, pct_doublets=0,
                 activity_pct=[.3, .4], activity_min=[.1, .1],
                 activity_max=[.7, .7], activity_identities=[[1, 2], [2, 3]],
                 max_activity_total=.8, diffexpprob=.025, diffexpdownprob=.025,
                 diffexploc=1.0, diffexpscale=1.0, delta_range_libloc=0):
        np.random.seed(self.seed)
        
        ncells_tosim = int(ncells*(1+pct_doublets))
        ndoublets = ncells_tosim - ncells        
        nactivites = len(activity_identities)
        nprograms = nidentities + nactivites
        self.gepnames = ['Identity_%d' % i for i in range(1, nidentities+1)] + ['Activity_%d' % i for i in range(1, nactivites+1)]         
        
        self.deparams = {'diffexpprob': diffexpprob, 'diffexpdownprob':diffexpdownprob,
                        'diffexploc':diffexploc, 'diffexpscale':diffexpscale}
        for k in self.deparams.keys():
            if not isarray(self.deparams[k]):
                self.deparams[k] = [self.deparams[k]]*nprograms
        
        print('Simulating gene params')
        self.geneparams = self.get_gene_params(ngenes=ngenes)

        print('Simulating DE')
        self.spectra = self.simulate_geps(ngenes=ngenes)
        
        print('Simulating cells')
        self.cellparams, self.usage = self.get_usages(ncells=ncells_tosim, 
                                               nidentities=nidentities, activity_pct=activity_pct,
                                               activity_min=activity_min, activity_max=activity_max,
                                               activity_identities=activity_identities,
                                               max_activity_total=max_activity_total,
                                               delta_range_libloc=delta_range_libloc)

        if pct_doublets > 0:
            print('Simulating doublets')
            self.simulate_doublets(ndoublets=ndoublets)
        
        print('Simulating cell-gene means')
        self.cellgenemean = self.get_cell_gene_means()
        
        print('Adjusting means')
        self.adjust_means_bcv()
        print('Simulating counts')
        self.simulate_counts()

        
    def simulate_geps(self, ngenes=10000, activity_deprob=None, activity_deloc=None,
                     activity_descale=None, activity_dedownprob=None):
        '''Sample differentially expressed genes and the DE factor for each
        cell-type group'''
        geps = self.gepnames
        spectra = pd.DataFrame(np.nan, index=geps, columns=self.geneparams.index)
        
        print(geps)
        for i,g in enumerate(geps):
            deprob = self.deparams['diffexpprob'][i]
            deloc = self.deparams['diffexploc'][i]
            descale = self.deparams['diffexpscale'][i]
            dedownprob = self.deparams['diffexpdownprob'][i]
            
            isDE = np.random.choice([True, False], size=ngenes,
                                      p=[deprob,1-deprob])
            #isDE[proggene] = False # Program genes shouldn't be differentially expressed between groups
            DEratio = np.random.lognormal(mean=deloc, sigma=descale, size=isDE.sum())
            DEratio[DEratio<1] = 1 / DEratio[DEratio<1]
            is_downregulated = np.random.choice([True, False],
                                            size=len(DEratio),
                                            p=[dedownprob,1-dedownprob])
            DEratio[is_downregulated] = 1. / DEratio[is_downregulated]
            all_DE_ratio = np.ones(ngenes)
            all_DE_ratio[isDE] = DEratio
            gep_mean = self.geneparams['gene_mean']*all_DE_ratio

            deratiocol = '%s.DEratio' % g
            groupmeancol = '%s.genemean' % g
            self.geneparams[deratiocol] = all_DE_ratio
            self.geneparams[groupmeancol] = gep_mean
            spectra.loc[g,:] = gep_mean.values
        
        return(spectra)
    
    
    def get_usages(self, ncells=10000, nidentities=8, pct_doublets=.05,
                        activity_pct=[.3, .4], activity_min=[.1, .1],
                        activity_max=[.7, .7], activity_identities=[[1, 2], [2, 3]],
                        max_activity_total=.8, delta_range_libloc=0):
        '''Sample cell group identities and library sizes'''
        index = ['Cell_%d' % i for i in range(1, ncells+1)]
        identity = (np.random.choice(nidentities, size=ncells) + 1)
        celldata = pd.DataFrame([identity], index=['Identity'], columns=index).T
        celldata['Identity'] = celldata['Identity'].astype(int)
        identity_usage = pd.get_dummies(celldata['Identity'])
        activity_usage = pd.DataFrame(0, index=celldata.index, columns=['Activity_%d' % (i+1) for i in range(len(activity_pct))])
    
        # Calculate unnormalized activity usage
        if len(activity_pct)>0:
            for i in range(len(activity_pct)):
                activity_id = 'Activity_%d' % (i+1)
                ind = celldata['Identity'].isin(activity_identities[i])
                ind = ind.index[ind]
                num_possibly_active = len(ind)
                uniform_rand = np.random.uniform(size=num_possibly_active)
                active_cells = ind[uniform_rand <= activity_pct[i]]
                activity_pct = np.random.uniform(low=activity_min[i], high=activity_max[i], size=len(active_cells))
                activity_usage.loc[active_cells, activity_id] = activity_pct            
            
        # Renormalize activity usages when they sum to greater than a threshold
        total_activity_usage = activity_usage.sum(axis=1)
        exceeding_total = total_activity_usage>max_activity_total
        renormed = activity_usage.loc[exceeding_total, :]
        renormed = renormed.div(renormed.sum(axis=1), axis=0)*max_activity_total
        activity_usage.loc[exceeding_total, :] = renormed
        total_activity_usage = activity_usage.sum(axis=1)
    
        # Normalize identity so sum with identity activity is 1
        identity_usage = identity_usage.multiply(1-total_activity_usage, axis=0)
        celldata['Identity_Usage'] = identity_usage.max(axis=1)
        for a in activity_usage.columns:
            celldata[a] = activity_usage[a]
        
        identity_usage.columns = ['Identity_%d' % i for i in identity_usage.columns]
        usage = pd.concat([identity_usage, activity_usage], axis=1) 

        ## Update with linkage between library size and gep usage
        ecounts = usage.dot(self.spectra).sum(axis=1)        
        new_min = self.libloc*(1-delta_range_libloc)
        new_max = self.libloc*(1+delta_range_libloc)
        libloc_param = rescale(ecounts, new_min, new_max)
        libsize = np.random.lognormal(mean=libloc_param, sigma=self.libscale,
                                      size=ncells)       

        ## Original, no correlation between GEP usage and library size
        #libsize = np.random.lognormal(mean=self.libloc, sigma=self.libscale,
        #                              size=ncells)
        
        celldata['Library_Size'] = libsize   
        
        
        return(celldata, usage)
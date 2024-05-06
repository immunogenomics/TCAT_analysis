import pandas as pd
import numpy as np

class scsim2:
    def __init__(self, seed=757578, mean_rate=7.68, mean_shape=0.34, libloc=7.64, libscale=0.78,
                expoutprob=0.00286, expoutloc=6.15, expoutscale=0.49, diffexpprob=.025, diffexpdownprob=.025,
                diffexploc=1.0, diffexpscale=1.0, bcv_dispersion=0.448, bcv_dof=22.087):
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
        self.diffexpprob = diffexpprob
        self.diffexpdownprob = diffexpdownprob
        self.diffexploc = diffexploc
        self.diffexpscale = diffexpscale
        self.bcv_dispersion = bcv_dispersion
        self.bcv_dof = bcv_dof

    def simulate(self, ncells=10000, ngenes=10000, nidentities=8, pct_doublets=0,
                 activity_pct=[.3, .4], activity_min=[.1, .1],
                 activity_max=[.7, .7], activity_identities=[[1, 2], [2, 3]],
                 max_activity_total=.8):
        np.random.seed(self.seed)
        
        ncells_tosim = int(ncells*(1+pct_doublets))
        ndoublets = ncells_tosim - ncells
        
        
        print('Simulating cells')
        self.cellparams, self.usage = self.get_usages(ncells=ncells_tosim, 
                                               nidentities=nidentities, activity_pct=activity_pct,
                                               activity_min=activity_min, activity_max=activity_max,
                                               activity_identities=activity_identities,
                                               max_activity_total=max_activity_total)
        
        print('Simulating gene params')
        self.geneparams = self.get_gene_params(ngenes=ngenes)

        print('Simulating DE')
        self.spectra = self.simulate_geps(ngenes=ngenes)
        
        '''
        if (self.nproggenes is not None) and (self.nproggenes > 0):
            print('Simulating program')
            self.simulate_program()
        '''

        if pct_doublets > 0:
            print('Simulating doublets')
            self.simulate_doublets(ndoublets=ndoublets)
        
        print('Simulating cell-gene means')
        self.cellgenemean = self.get_cell_gene_means()
        
        print('Adjusting means')
        self.adjust_means_bcv()
        print('Simulating counts')
        self.simulate_counts()

    def simulate_counts(self):
        '''Sample read counts for each gene x cell from Poisson distribution
        using the variance-trend adjusted updatedmean value'''
        self.counts = pd.DataFrame(np.random.poisson(lam=self.updatedmean),
                                   index=self.updatedmean.index, columns=self.updatedmean.columns)

    def adjust_means_bcv(self):
        '''Adjust cellgenemean to follow a mean-variance trend relationship'''
        self.bcv = self.bcv_dispersion + (1 / np.sqrt(self.cellgenemean))
        chisamp = np.random.chisquare(self.bcv_dof, size=self.cellgenemean.shape[1])
        self.bcv = self.bcv*np.sqrt(self.bcv_dof / chisamp)
        self.updatedmean = np.random.gamma(shape=1/(self.bcv**2),
                                           scale=self.cellgenemean*(self.bcv**2))
        self.bcv = pd.DataFrame(self.bcv, index=self.cellgenemean.index, columns=self.cellgenemean.columns)
        self.updatedmean = pd.DataFrame(self.updatedmean, index=self.cellgenemean.index,
                                        columns=self.cellgenemean.columns)


    def simulate_doublets(self, ndoublets):
        ## Select doublet cells and determine the second cell to merge with
        ncells = self.usage.shape[0]
        nfinal = ncells - ndoublets
        d_ind = sorted(np.random.choice(self.usage.index[:nfinal], ndoublets, replace=False))
        extraind = self.usage.index[(-1*ndoublets):]        
        self.cellparams['is_doublet'] = False
        self.cellparams.loc[d_ind, 'is_doublet'] = True
        self.cellparams['Doublet_Identity2'] = -1
        self.cellparams.loc[d_ind, 'Doublet_Identity2'] = self.cellparams.loc[extraind, 'Identity'].values
        self.usage.loc[d_ind, :] = self.usage.loc[d_ind, :] + self.usage.loc[extraind, :].values
        self.usage.loc[d_ind, :] = self.usage.loc[d_ind, :].div(self.usage.loc[d_ind, :].sum(axis=1), axis=0)
        self.usage = self.usage.iloc[:nfinal, :]
        self.cellparams = self.cellparams.iloc[:nfinal, :]
        activity_geps = [x for x in self.usage.columns if 'Activity' in x]
        for g in activity_geps:
            self.cellparams[g] = self.usage[g]


    def get_cell_gene_means(self):
        '''Calculate each gene's mean expression for each cell while adjusting
        for the library size'''
        cellgenemean = self.usage.dot(self.spectra)
        cellgenemean = cellgenemean.multiply(self.cellparams['Library_Size'] / cellgenemean.sum(axis=1), axis=0)
        return(cellgenemean)


    def get_gene_params(self, ngenes=10000):
        '''Sample each genes mean expression from a gamma distribution as
        well as identifying outlier genes with expression drawn from a
        log-normal distribution'''
        basegenemean = np.random.gamma(shape=self.mean_shape,
                                       scale=1./self.mean_rate,
                                       size=ngenes)

        is_outlier = np.random.choice([True, False], size=ngenes,
                                      p=[self.expoutprob,1-self.expoutprob])
        outlier_ratio = np.ones(shape=ngenes)
        outliers = np.random.lognormal(mean=self.expoutloc,
                                       sigma=self.expoutscale,
                                       size=is_outlier.sum())
        outlier_ratio[is_outlier] = outliers
        gene_mean = basegenemean.copy()
        median = np.median(basegenemean)
        gene_mean[is_outlier] = outliers*median
        self.genenames = ['Gene%d' % i for i in range(1, ngenes+1)]
        geneparams = pd.DataFrame([basegenemean, is_outlier, outlier_ratio, gene_mean],
                                  index=['BaseGeneMean', 'is_outlier', 'outlier_ratio', 'gene_mean'],
                                 columns=self.genenames).T
        return(geneparams)


    def get_usages(self, ncells=10000, nidentities=8, pct_doublets=.05,
                        activity_pct=[.3, .4], activity_min=[.1, .1],
                        activity_max=[.7, .7], activity_identities=[[1, 2], [2, 3]],
                        max_activity_total=.8):
        '''Sample cell group identities and library sizes'''
        index = ['Cell_%d' % i for i in range(1, ncells+1)]
        libsize = np.random.lognormal(mean=self.libloc, sigma=self.libscale,
                                      size=ncells)
        identity = (np.random.choice(nidentities, size=ncells) + 1)
        celldata = pd.DataFrame([identity, libsize], index=['Identity', 'Library_Size'], columns=index).T
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
        return(celldata, usage)


    def simulate_geps(self, ngenes=10000, activity_deprob=None, activity_deloc=None,
                     activity_descale=None, activity_dedownprob=None):
        '''Sample differentially expressed genes and the DE factor for each
        cell-type group'''
        geps = list(self.usage.columns)
        spectra = pd.DataFrame(np.nan, index=geps, columns=self.geneparams.index)
        for g in geps:
            if ('Activity' in g) and (activity_deprob is not None):
                deprob = activity_deprob
                deloc = activity_deloc
                descale = activity_descale
                dedownprob = activity_dedownprob

            else:
                deprob = self.diffexpprob
                deloc = self.diffexploc
                descale = self.diffexpscale
                dedownprob = self.diffexpdownprob
            
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

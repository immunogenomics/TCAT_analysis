# scsim2
Simulate single-cell RNA-SEQ data using the [Splatter statistical framework](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1305-0) but implemented in python with additional features. scsim2 is an extension of [scsim](https://github.com/dylkot/scsim) but now allowing an arbitrary number of activity programs to be simulated. 

Added on top of this is SCsim_DiffScale within scsim2_diff_librarysize.py which allows GEPs to have different abundances of RNA content which corresponds to having cells expressing certain GEPs having higher overall library sizes. For example:

```
deprob = [.025]*7 + [.05]*2 + [.025]*2
deloc = [1.5]*7 + [2.00]*2 + [1.5]*2
simmod_catalog = SCsim_DiffScale(seed=112)
simmod_catalog.simulate(ncells=20000, ngenes=10000, nidentities=7, activity_pct=[.3]*4,
                activity_min=[.1]*4, activity_max=[.7]*4,
                activity_identities=[list(range(1,8))]*4,
                pct_doublets=.05, diffexploc=deloc, diffexpscale=1, diffexpprob=deprob,
                delta_range_libloc=0)

simmod_catalog.spectra.sum(axis=1)
```

Would cause the 8th and 9th programs to have a larger library size because of more incorporated genes and higher average differential expression for those genes.
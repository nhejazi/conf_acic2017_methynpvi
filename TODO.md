## 15 May 2017

* ...

---

## 10 May 2017

* ~make full data available for use with `r2weight` (Nima)~
* use David's `r2weight` to generate optimal IQ scores for use with the NPVI
	parameter (Ivana)
* run standard TMLE for the NPVI parameter on methylation data -- perhaps using
	just a single IQ measure as outcome -- using the 5% quantile to define the
	null range, either on a site-by-site basis or across genome (Andre)
* write script for using the data-adaptive framework to obtain the null range
	for use with the NPVI parameter -- this should be over the whole genome;
	otherwise, have to estimate cutoff for each site across each CV fold, making
	the procedure potentially computationally prohibitive (Nima)

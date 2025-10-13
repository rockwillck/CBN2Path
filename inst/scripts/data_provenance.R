## Source: 
## CT-CBN C package and the external data (in the extdata folder) are obtained from ETH Zurich (https://bsse.ethz.ch/cbg/software/ct-cbn.html)
## A collection of genetic data sets (in the RawData folder) are provided from three cancer studies conducted by the Vogelstein lab. 
## The data consists of three tables with mutations and a mapping of genes to 12 core pathways (Parsons et al, 2008). 
## The package also contains the python module cbn.py for parsing the raw data and for bootstrap and permutation analyses. The module requires the numpy python module.

##  Raw data files:
    # - RawData/JonesS2008.txt    #From Suppl. Table S3, Jones et al., Science 2008 [1]
    # - RawData/WoodS2007.txt	    #From Wood et al. Science 2007 [2]
    # - RawData/ParsonsS2008.txt  #From Suppl. Table S4, Parsons et al., Science 2008 [3]
        # Each line in each table shows a mutation.

##  RawData/CoreGroup2Gene.txt
        # contains the mapping of core pathways to genes (HUGO identifiers).

##  RawData/cbn.py 
        # is a python module for analyzing these data sets, and for manipulating the output from h-cbn. Each function contains a short help string. Requires numpy.

##  RawData/example.py
        # is a python script creating four .pat files and .poset from Jones2008.txt as input for h-cbn. To run:

## References:
   # [1] Jones S, et al. Core signaling pathways in human pancreatic cancers revealed by global genomic analyses. Science. 2008 Sep 26;321(5897):1801-6. doi: 10.1126/science.1164368. Epub 2008 Sep 4. PMID: 18772397; PMCID: PMC2848990.
   # [2] Wood LD, et al. (2007) The genomic landscapes of human breast and colorectal cancers. Science 318: 1108–1113.
   # [3] Parsons DW, et al. (2008) An integrated genomic analysis of human glioblastoma multiforme. Science 321: 1807–1812. 



## Formats
   # filestem.pat    - Mutational patterns (genotypes), unless N > 0
   # filestem.poset  - Event poset used if -e is not set; if -e is set, the file is used for determining the number of events as specified in the first row
   # filestem.lambda - Model parameters, if N > 0

## No modifications were made to these files. These files are used in examples and unit tests.

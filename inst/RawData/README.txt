Raw data and software for "The Temporal Order of Genetic and Pathway Alterations in Tumorigenesis"

Raw data files:
- JonesS2008.txt    #From Suppl. Table S3, Jones et al., Science 2008
- WoodS2007.txt	    #From Wood et al. Science 2007
- ParsonsS2008.txt  #From Suppl. Table S4, Parsons et al., Science 2008
contain tables in which each line shows a mutation.

- CoreGroup2Gene.txt
contains the mapping of core pathways to genes (HUGO identifiers).

- cbn.py 
is a python module for analyzing these data sets, and for manipulating the output from h-cbn. Each function contains a short help string. Requires numpy.

- example.py
is a python script creating four .pat files and .poset from
Jones2008.txt as input for h-cbn. To run:
$ python example.py

You can then run h-cbn on jones_genes.pat with
$ h-cbn -f jones_genes -s

See the example.py for more options.

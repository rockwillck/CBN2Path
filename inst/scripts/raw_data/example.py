#!/usr/bin/env python
import cbn

# 1. Genes
# read tumor data for Jones et al. Output: Dictionary
jones = cbn.read_tumors("JonesS2008.txt")[0]

# convert to binary pattern with gene frequencies > 5% (input to h-cbn)
jones_genes_pat = cbn.tumors2pat(jones) #numpy.array
jones_genes_pat = jones_genes_pat[:,jones_genes_pat.mean(0) > 0.05]
cbn.write_pat(jones_genes_pat, "jones_genes.pat")

# compute linear poset (use as start for simulate annealing)
poset = cbn.linear_poset(jones_genes_pat) #numpy.array
cbn.write_poset(poset, "jones_genes.poset")

# 2. Core pathways
# read core pathway mapings
C2G, G2C = cbn.read_coregroups("CoreGroup2Gene.txt") # dictionaries

# generate pattern and linear poset
jones_core = cbn.tumor2core(jones, G2C) # Tumor->core_pathway dict
jones_core_pat = cbn.tumors2pat(jones_core)
cbn.write_pat(jones_core_pat, "jones_core.pat")
cbn.write_poset(cbn.linear_poset(jones_core_pat), "jones_core.poset")

# 3. Bootstrapping
# A simple bootstrap
bones_core_pat = cbn.bootstrap(jones_core_pat)
cbn.write_pat(bones_core_pat, "jones_core_bootstrap.pat")
cbn.write_poset(cbn.linear_poset(bones_core_pat), "jones_core_bootstrap.poset")


# 4. Permutations
# Randomize tumor dictionary
rones = cbn.randomize_tumors(jones)
# Compute pattern
rones_core_pat = cbn.tumors2pat(cbn.tumor2core(jones, G2C))
cbn.write_pat(rones_core_pat, "jones_core_permuted.pat")
cbn.write_poset(cbn.linear_poset(rones_core_pat), "jones_core_permuted.poset")

note = """You can now run h-cbn on all datasets. For example
$ h-cbn -f jones_genes
Likelihood of the linear poset for the genetic data.

$ h-cbn -f jones_core -s
Simulated annealing search for the core pathway.

$ h-cbn -f jones_core_bootstrap -s -N 200
Simulated annealing search for the bootstrap with 200 steps.

$ h-cbn -f jones_core_permuted -m
Print the most likely progression in the permuted data set.

Check out h-cbn -f for more options."""

print note

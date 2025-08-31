
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CBN2Path

<!-- badges: start -->
<!-- badges: end -->

## Authors

William Choi-Kim and Sayed-Rzgar Hosseini

## Abstract

Tumorigenesis is a stepwise process that is driven by a sequence of
molecular changes forming pathways of cancer progression. Conjunctive
Bayesian Networks are probabilistic-graphical models designed for the
analysis and modeling of these pathways \[1\]. CBN models have evolved
into different varieties such as CT-CBN \[2\], H-CBN \[3\], B-CBN \[4\]
and R-CBN \[5\] each addressing different aspects of this task. However,
the software corresponding to these methods are not well-integrated as
they are implemented in different languages with heterogeneous input and
output formats. This necessitates a unifying platform that integrates
these models and enables standardization of the input and output formats
to facilitate the downstream pathway analysis and modeling. Evam-tools
\[6\] is an R package, which has taken the initial steps towards this
end. However, it partially serves this purpose, as it does not include
the B-CBN model and the recently developed R-CBN algorithm, which
focuses on the robust inference of cancer progression pathways \[5\].
Importantly, the B-CBN and R-CBN algorithms for pathway quantification
require exhaustive consideration and weighting of all the potential
dependency structures (posets) within mutational quartets. This entails
re-implementation of the CBN models and adjustment of the downstream
pathway analysis and modeling functions. Therefore, here we introduce
**CBN2Path** R package that not only includes the original
implementation of the CBN models (e.g. CT-CBN and H-CBN) in a unifying
interface, but it also accommodates the necessary modifications to
support the robust CBN algorithms (e.g. B-CBN and R-CBN). In summary,
CBN2Path is an R package that supports robust quantification, analysis
and visualization of cancer progression pathways from cross-sectional
genomic data, and so we anticipate that it will be a widely-used package
in the future.

## Installation

### GSL

To install the `CBN2Path` R package, you first need to install the
`gsl`:

**Install GSL with homebrew on Mac:**

*If you don’t have homebrew, run the following command in your
terminal/console:*

``` bash
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```

Then, also in terminal:

``` bash
brew install gsl
```

Note that if `gsl` was installed using any method other than Homebrew,
you need to uninstall `gsl`, and then reinstall it using Homebrew (see
<https://brew.sh> if you have not installed Homebrew yet).

**Install GSL on Linux:**

In your shell:

``` bash
sudo apt-get install libgsl-dev
```

**On Linux, if the `ggraph` dependency fails, run the following in your
shell:**

``` bash
sudo apt install libfontconfig1-dev
```

This appears to fix a sysfonts issue. We’re not sure why this is
necessary.

**On Windows, we suggest installing RTools (which includes a
distribution of GSL):**

Download RTools from
[here](https://cran.r-project.org/bin/windows/Rtools/) and proceed with
installation.

### Package Install

**Make sure to restart R before proceeding.**

Then, you can install the development version of `CBN2Path` by running
the following in R:

**Linux and Mac**

``` r
remotes::install_github("rockwillck/CBN2Path", build_vignettes = TRUE)
```

**Windows**

``` r
remotes::install_github("rockwillck/CBN2Path", build_vignettes = FALSE)
```

## Windows Support

Windows support for `CBN2Path` is **limited**. Functions will be missing
key functionality; the CBN models developed at ETH-Zurich that
`CBN2Path` is based on don’t support Windows inherently.

## Usage

To learn how to use different CBN models and their associated pathway
analysis and visualization functions in the `CBN2Path` R package, please
run:

``` r
vignette("CBN2Path")
```

## Cite our work

If you use the CBN2Path package, please cite the paper formally as
follows:

Choi-Kim W and Hosseini SR. CBN2Path: an R/Bioconductor package for the analysis of cancer progression pathways using Conjunctive Bayesian Networks [version 1; peer review: awaiting peer review]. F1000Research 2025, 14:834 (https://doi.org/10.12688/f1000research.168810.1)

## References

\[1\] Beerenwinkel, et al. Conjunctive Bayesian Networks. Bernoulli,
13(4):893–909, November 2007. ISSN 1350-7265. doi:
<https://doi.org/10.3150/07-BEJ6133>.

\[2\] Beerenwinkel and Sullivant. Markov models for accumulating
mutations. Biometrika, 96 (3):645–661, September 2009. ISSN 0006-3444,
1464-3510. doi: <https://doi.org/10.1093/biomet/asp023>.

\[3\] Gerstung, et al. Quantifying cancer progression with conjunctive
Bayesian networks. Bioinformatics (Oxford, England), 25(21):2809–2815,
November 2009. doi: <https://doi.org/10.1093/bioinformatics/btp505>.

\[4\] Sakoparnig and Beerenwinkel. Efficient sampling for Bayesian
inference of conjunctive Bayesian networks. Bioinformatics,
28(18):2318–2324, September 2012. ISSN 1367-4811, 1367-4803. doi:
<https://doi.org/10.1093/bioinformatics/bts433>.

\[5\] Hosseini. Robust inference of cancer progression pathways using
Conjunctive Bayesian Networks, BioRxiv. July 2025. doi:
<https://doi.org/10.1101/2025.07.15.663924>.

\[6\] Diaz-Uriarte and Herrera-Nieto. EvAM-Tools: tools for evolutionary
accumulation and cancer progression models. Bioinformatics, 38(24):
5457–5459, December 2022. ISSN 1367-4803, 1367-4811. doi:
<https://doi.org/10.1093/bioinformatics/btac710>.

\[7\] Hosseini, et al. Estimating the predictability of cancer
evolution. Bioinformatics, 35 (14):i389–i397, July 2019. ISSN 1367-4803,
1367-4811. doi: <https://doi.org/10.1093/bioinformatics/btz332>.

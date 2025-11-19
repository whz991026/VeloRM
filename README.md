
# VeloRM

<!-- badges: start -->
<!-- badges: end -->
ABSTRACT
RNA modification plays a critical role in regulating the function and fate of RNA molecules, yet its dynamics across the RNA life cycle remain under-explored. We present here VeloRM, a computational framework that disentangles the pre-splicing and post-splicing epitranscriptomes, enabling the inference of RNA modification kinetics and the prediction of future epitranscriptome states at single-nucleotide and single-cell resolution. Applying VeloRM to three single-cell epitranscriptome datasets of m6A and A-to-I, we uncover distinct modification patterns on the pre-splicing and post-splicing RNAs. Notably, N6-methyladenosine (m6A) is more heavily decorated on pre-splicing RNAs than post-splicing RNAs near the splicing junctions. Importantly, VeloRM rigorously quantifies the velocity of RNA modification, provides direct insight into the future epitranscriptome landscape of individual cells, aligning well with known data on single-cell trajectories within the cell cycle and during cell differentiation. Together, our study broadens the current understanding of the dynamic epitranscriptome by establishing a rigorous framework for tracing its kinetics across the RNA life cycle in a dynamic process spanning the pre-splicing and post-splicing RNA states.

Keywords: single-cell, epitranscriptome, RNA velocity, m6A, A-to-I, RNA modification dynamics, scDART-seq, splicing, RNA regulation


## Installation

You can install the development version of AEEIP from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("whz991026/VeloRM")
```

## Example

The example can be found in vignette/doc.

## data
(https://doi.org/10.6084/m9.figshare.28955183)


# VeloRM

<!-- badges: start -->
<!-- badges: end -->

RNA modifications play critical roles in regulating gene expression, yet their dynamic behavior at single-cell resolution remains largely unexplored. To address this, we developed VeloRM, a computational framework that models the kinetics of RNA modifications—such as m⁶A methylation and A-to-I editing—by quantifying their gain and loss over pseudotime at single-nucleotide resolution. Applying VeloRM to scDART-seq data, we uncovered distinct regulatory programs: most m⁶A modifications occur prior to splicing, though a subset arises post-splicing. In contrast, A-to-I editing cannot be clearly distinguished as occurring either prior to or after splicing. VeloRM accurately reconstructs modification-associated cellular trajectories and identifies stage-specific modification patterns that influence RNA fate. Notably, m⁶A levels correlate with RNA decay, with most regulatory effects observed in precursor RNA, although some modification sites show a positive association with mature RNA abundance. Together, these findings provide new insights into the dynamic regulation of the epitranscriptome at single-cell resolution and establish a powerful framework for studying RNA modification kinetics in development and disease.

Keywords: single-cell, epitranscriptomics, RNA velocity, m⁶A, A-to-I editing, RNA modification dynamics, scDART-seq, splicing, RNA regulation



## Installation

You can install the development version of VeloRM from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("whz991026/VeloRM")
```
## Code and data
10.6084/m9.figshare.28955183



## Example
see vignette

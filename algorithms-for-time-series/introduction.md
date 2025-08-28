# Introduction

In this section, we present the data structures and algorithms for solving the following problems on time series data:

1. computing the similarity of multiple time series,
2. finding an optimal or approximate matching of multiple time series,
3. computing the density of multiple time series, and
4. computing the diameter of multiple time series.

The last three problems are built upon the first problem, namely the similarity computation. We provide not only our own data structures and algorithms but also a review and implementation of prior results for comparison. This section, therefore, offers a comprehensive understanding of the field.



## [Data Structure](data-structures/)

### Data

* [**Polygonal curve**](data-structures/polygonal-curve.md): This representation transforms discrete time series into continuous geometric objects, enabling the use of geometric distance measures such as the Fréchet distance.

### Tools for algorithms



***

## Algorithms

### Similarity

* **Fréchet Distance**: a well-established similarity measure between two polygonal curves (time series).

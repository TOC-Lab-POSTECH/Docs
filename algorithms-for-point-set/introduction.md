# Introduction

In this section, we present the data structures and algorithms for solving the following problems on a point set:

1. computing the clustering on a point set,
2. computing the data matching on a point set

These problems arise in a wide range of applications, including trajectory analysis, pattern recognition, and large-scale data mining.

### [Data Structure](https://toc-lab-postech.gitbook.io/postech-toc-lab-sw-starlab/algorithms-for-point-set/data-structures)

#### Data

* [**Range Counting Oracle**](https://toc-lab-postech.gitbook.io/postech-toc-lab-sw-starlab/algorithms-for-point-set/data-structures/range-counting-oracle): A geometric data structure that supports efficient range counting queries on a point set. This oracle serves as a fundamental building block for sublinear-time algorithms by enabling fast access to the data set without enumerating individual points.

#### Tools for algorithms

* k-median++

### Algorithms

#### Clustering

* **k-median clustering**: A classical clustering problem that seeks to place (k) centers so as to minimize the sum of distances from each point to its nearest center. While k-median is NP-hard in general, it serves as a fundamental abstraction for facility location, summarization, and data reduction. In this work, we focus on the approximation algorithm of k-median problem, enabling efficient clustering even for very large point sets.

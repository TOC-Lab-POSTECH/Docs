# Data Structures

In this section, we present several data structures for solving some problems on time series data. We will provide the data structures for abstracting input data and the data structures that will be used as a tool or black box for the algorithms.



## Data

* [**Polygonal curve**](polygonal-curve.md): Most time series data are represented as _polygonal curves_, which are piecewise linear curves obtained by sequentially connecting data points in temporal order. Formally, given a sequence $$P = (p_1, p_2, \ldots, p_n)$$ where each $$p_i$$ is a point, the polygonal curve $$\mathcal{P}$$ is defined as the union of the line segments $$\overline{p_i p_{i+1}}$$ for $$1 \leq i < n$$. This representation transforms discrete time series into continuous geometric objects, enabling the use of geometric distance measures such as the Fréchet distance.



## Tools for algorithms

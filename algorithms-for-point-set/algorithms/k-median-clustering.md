# k-Median Clustering

This document describes the **k-median clustering pipeline** implemented in this repository. The pipeline combines a practical **k-median++ heuristic** with a **sublinear representative-set construction** based on range counting, enabling scalable clustering on large point sets.

The implementation is designed to evaluate both **baseline clustering** on the original data and a **compressed pipeline** that operates on representative sets.

***

### k-Median Clustering (Conceptual Overview)

#### Definition

Given a set of points

* $$P = \{p_1, \dots, p_n\} \subset \mathbb{R}^2$$

and an integer

* $$k,$$

the **k-median problem** asks for a set of centers

* $$C = \{c_1, \dots, c_k\}$$

minimizing the objective

* $$\sum_{p \in P} w(p) \cdot \min_{c \in C} |p - c|.$$

Here $$w(p)$$ is a non-negative weight associated with point $$p$$. The problem is NP-hard even in the plane, motivating the use of approximation algorithms and heuristics.

***

#### Sublinear Algorithm of Monemizadeh Monemizadeh

The representative-set construction implemented in this repository is **inspired by the sublinear geometric framework of Morteza Monemizadeh** for facility location.

In particular, the design closely follows the ideas in:

> _M. Monemizadeh, “Facility Location in the Sublinear Geometric Model,”_\
> &#xNAN;_&#x41;PPROX/RANDOM 2023._

At a high level, the insight is the following:

* Space is decomposed into a hierarchy of grids (quadtree levels).
* A cell is classified as **dense** or **sparse** based on a scale-dependent threshold with range counting oracle.
* Sparse cells can be safely **collapsed into single weighted representatives** without significantly affecting the objective value.
* The total number of such sparse cells is provably bounded by $$O\!\left(\frac{k \log n}{\varepsilon^3}\right).$$

In the original paper establishes that replacing each sparse cell by a single representative preserves the total cost up to a $$1 + \varepsilon$$ factor, provided the density threshold is chosen appropriately.

Our implementation adapts this idea, and replaces oracle-based sampling with explicit **range counting via a RangeTree**. Additionally, for practical efficiency, the k-median++ algorithm, which is a modified version of the k-means++ algorithm, was implemented and used.

***

### Repository

{% embed url="https://github.com/TOC-Lab-POSTECH/Approximation-k-median-clustering" %}

***



### Class Synopsis

```cpp
#include "Representative.h"

class RepresentativeSetBuilder {
 public:
  RepresentativeSetBuilder(const RangeTree& tree,
                           const RjConfig& cfg);

  // Build representative set R_j
  // Returns (success flag, representative points)
  std::pair<bool, std::vector<Point>> build();

 private:
  void processCell(const QuadCell& c);
  int countPointsInCell(const QuadCell& c) const;
};
```

The class `RepresentativeSetBuilder` implements the **representative-set construction** corresponding to a single cost guess $$r_j.$$

Given a fixed $$r_j$$, the builder performs a quadtree-based spatial decomposition and identifies a bounded collection of **sparse cells** whose contributions can be safely aggregated.

***

### Dependencies

#### STL

* `<vector>` – storage of sparse cells and representative points.
* `<utility>` – return of `(success, result)` pairs.
* `<cmath>` – logarithms and threshold computations.
* `<algorithm>` – basic numeric utilities.

#### Library

* `"Representative.h"` – representative-set construction logic.
* `"RangeTree.h"` – orthogonal range counting oracle.
* `"Point.h"` – weighted point abstraction.

No external geometry libraries are required.

***

### Design Notes (Representative-Centric)

* **Representative-based compression**
  * The algorithm constructs a weighted representative set (R\_j) for a given cost guess (r\_j), replacing many original points by a single representative per sparse region.
* **Quadtree decomposition**
  * The input domain is covered by the smallest power-of-two square enclosing the bounding box.
  * Space is recursively subdivided into quadtree cells of side length (2^i).
* **Range-counting oracle**
  * Cell densities are computed using a `RangeTree`, avoiding explicit iteration over points.
  * This enables sublinear behavior in practice.
* **Density threshold**
  * A cell at level (i) is classified as sparse if $$n_c < \frac{\delta_{\text{k-med}} \cdot r_j}{2^i},$$ where $$\delta_{\text{k-med}} = \Theta!\left(\frac{k \log n}{\varepsilon^3}\right)$$.
  * Sparse cells can be safely collapsed into weighted representatives.
* **Abort mechanism**
  * The construction aborts if the number of sparse cells exceeds $$O\!\left(\frac{k \log n}{\varepsilon^3}\right),$$ indicating that the current cost guess $$r_j$$ is too small.
* **Practical complexity**
  * Each cell is processed once and each density query costs $$O(\log n)$$.
  * In practice, runtime depends on the number of sparse cells rather than on $$n$$, yielding sublinear behavior.

***

### Overall Structure

### Example Usage

```cpp
RangeTree tree(points);
RjConfig cfg{n, k, eps, rj};

RepresentativeSetBuilder builder(tree, cfg);
auto [ok, Rj] = builder.build();

if (ok) {
  // Run k-median++ on Rj
}
```

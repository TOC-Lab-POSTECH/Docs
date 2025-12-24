# k-median++

The **k-median++** module implements a practical heuristic for the geometric **k-median clustering problem**.\
It is inspired by the k-means++ seeding strategy, but adapted to the k-median objective by sampling proportional to **linear distances** rather than squared distances.

In our framework, k-median++ is used both as:

1. a **baseline clustering algorithm** on the original point set, and
2. a **core subroutine** applied to compressed representative sets.

Conceptually, the algorithm follows a two-phase structure:

1. **k-median++ initialization**
2. **Lloyd-style local refinement using geometric medians**

***

#### Repository\*\*\*

### Problem Setting

Given a set of points $$P \subset \mathbb{R}^2$$ and an integer $$k$$, the **k-median** problem asks for a set of centers $$C$$ with $$|C| = k$$ minimizing

$$cost(P, C) = \sum_{p \in P} w(p), d(p, C), \quad d(p, C) = \min_{c \in C} |p - c|.$$

Here $$w(p)$$ is a non-negative weight associated with point $$p$$ (default $$w(p)=1$$).\
The problem is NP-hard even in the plane, motivating the use of efficient approximation algorithms and heuristics.

***

### Algorithm Overview

Our implementation consists of the following steps:

1. **k-median++ initialization**\
   Select initial centers using weighted distance-based sampling.
2. **Assignment step**\
   Assign each point to its nearest center.
3. **Update step (1-median computation)**\
   Update each center to the weighted geometric median of its assigned points.
4. **Iteration until convergence**\
   Repeat assignment and update until stabilization.

***

### k-median++ Initialization

The initialization follows the k-means++ paradigm, adapted to the k-median objective.

#### First Center

The first center is sampled randomly proportional to point weights:

$$
\Pr[p \text{ chosen}] \propto w(p).
$$

#### Remaining Centers

Each subsequent center is sampled according to

$$
\Pr[p \text{ chosen}] \propto w(p) \cdot D(p),
$$

where $$D(p)$$ is the distance from $$p$$ to its nearest already chosen center.

Unlike k-means++, **squared distances are not used**, which better aligns the sampling with the k-median objective.

If all distances are zero (degenerate case), remaining centers are chosen uniformly at random.

***

### Lloyd-Style Refinement

After initialization, the algorithm performs Lloyd-style iterations.

#### Assignment Step

Each point is assigned to the nearest center.

If assignments do not change between iterations, the algorithm terminates early.

#### Update Step (Geometric Median)

For each cluster, the center is updated to the **weighted geometric median** of its assigned points.

The geometric median is approximated using **Weiszfeld’s algorithm**, initialized at the weighted centroid of the cluster.

Empty clusters are handled explicitly by reinitializing the center to a random input point.

***

### Geometric Median via Weiszfeld’s Algorithm

Given a cluster (S), the geometric median is approximated iteratively as

$$
x^{(t+1)} = \frac{\sum_{p \in S} \frac{w(p), p}{|p - x^{(t)}|}} {\sum_{p \in S} \frac{w(p)}{|p - x^{(t)}|}}.
$$



If the current estimate becomes extremely close to a data point, the algorithm snaps to that point and terminates to ensure numerical stability.

***

### Class Synopsis

```cpp
#include "k_median++.h"

std::vector<Point> k_median_pp(
    const std::vector<Point>& pts,
    int k,
    int maxIters,
    int maxWeiszfeldIters,
    double tol,
    uint32_t seed
);

```

***

### Dependencies

#### STL

* `<vector>`
  * Used to store input points, centers, cluster assignments, and intermediate clusters.
* `<cmath>`
  * Provides `sqrt` for Euclidean distance computation and geometric median updates.
* `<random>`
  * Used for randomized sampling in k-median++ initialization (`std::mt19937`, `uniform_real_distribution`, `uniform_int_distribution`).
* `<limits>`
  * Provides `std::numeric_limits<double>::infinity()` for distance and cost initialization.
* `<algorithm>`
  * Used for `std::min`, `std::max`, and `std::lower_bound` in weighted random sampling.

***

### Design Notes

* **Weighted formulation**
  * All steps of the algorithm explicitly incorporate point weights $$w(p)$$.
  * This allows k-median++ to be applied seamlessly to compressed representative sets, where each representative aggregates multiple original points.
* **Initialization strategy**
  * The k-median++ seeding uses linear distances instead of squared distances.
  * This choice aligns the sampling distribution with the k-median objective and empirically improves stability compared to k-means++-style initialization.
* **Geometric median updates**
  * Each center update computes a weighted geometric median using Weiszfeld’s algorithm.
  * The algorithm is initialized at the weighted centroid of the cluster to accelerate convergence.
* **Robust handling of degeneracies**
  * If all points are identical during initialization, remaining centers are filled by random sampling.
  * Empty clusters arising during Lloyd iterations are explicitly handled by reinitializing the center to a random input point, preventing collapse.
* **Early termination**
  * The algorithm terminates early when assignments no longer change or when all centers move by at most a tolerance threshold.
  * This ensures good practical performance despite the worst-case theoretical bounds.

***

### Complexity

Let $$n = |P|$$ be the number of points and $$k$$ be the number of centers.

* **Initialization**
  * Each new center requires scanning all points to compute distances to the current centers.
  * Total cost: $$O(nk)$$.
* **Lloyd-style iterations**
  * Assignment step: $$O(nk)$$.
  * Update step: geometric median computation per cluster, run for a fixed number of iterations.
  * Each Lloyd iteration costs $$O(nk)$$.

Overall, the total running time is $$O(T \cdot n \cdot k)$$, where $$T$$ is the number of Lloyd iterations, which is typically small in practice.

***

{% tabs %}
{% tab title="Data Structure" %}
<table><thead><tr><th width="160">Name</th><th width="200">Type</th><th>Description</th></tr></thead><tbody><tr><td><code>Point</code></td><td><code>struct</code></td><td>Represents a weighted point in the plane, storing coordinates <code>(x, y)</code> and a non-negative weight <code>w</code>.</td></tr><tr><td><code>std::vector&#x3C;Point></code></td><td><code>container</code></td><td>Used to store input points, centers, clusters, and representative sets.</td></tr></tbody></table>
{% endtab %}

{% tab title="Public Interface" %}
<table><thead><tr><th width="430">Name</th><th width="140">Input</th><th width="140">Output</th></tr></thead><tbody><tr><td><code>std::vector&#x3C;Point> k_median_pp( const std::vector&#x3C;Point>&#x26; pts, int k, int maxIters, int maxWeiszfeldIters, double tol, uint32_t seed)</code></td><td>weighted point set, parameters</td><td>selected centers</td></tr></tbody></table>

This function executes the full k-median++ algorithm, including initialization and Lloyd-style refinement, and returns the final set of centers.
{% endtab %}

{% tab title="Internal Functions" %}
<table><thead><tr><th width="430">Name</th><th width="160">Purpose</th><th width="120">Complexity</th></tr></thead><tbody><tr><td><code>double euclidean_distance(const Point&#x26;, const Point&#x26;)</code></td><td>Compute Euclidean distance between two points</td><td><span class="math">O(1)</span></td></tr><tr><td><code>std::vector&#x3C;Point> k_median_pp_init(const std::vector&#x3C;Point>&#x26;, ...)</code></td><td>k-median++ probabilistic initialization</td><td><span class="math">O(nk)</span></td></tr><tr><td><code>Point geometric_median_weiszfeld(const std::vector&#x3C;Point>&#x26;, ...)</code></td><td>Approximate weighted geometric median of a cluster</td><td><span class="math">O(|C|)</span> per iteration</td></tr><tr><td><code>double k_median_cost(const std::vector&#x3C;Point>&#x26;, const std::vector&#x3C;Point>&#x26;)</code></td><td>Evaluate weighted k-median objective</td><td><span class="math">O(nk)</span></td></tr></tbody></table>

Here $$n = |P|$$ is the number of points and $$|C|$$ denotes the size of a cluster.
{% endtab %}
{% endtabs %}

### Usage in Our Framework

Within our framework, k-median++ is used in two distinct roles:

* **Baseline clustering**
  * Applied directly to the original input point set.
  * Serves as a reference for evaluating approximation quality and runtime improvements.
* **Clustering on representative sets**
  * Applied to compressed representative sets $$R_j$$ produced by the range-based preprocessing stage.
  * This significantly reduces the effective input size while preserving clustering cost up to a controlled approximation error.

By separating these two use cases, we can clearly distinguish the effects of data compression from the behavior of the clustering heuristic itself.

***

### Example Usage

```cpp
std::vector<Point> points = load_points("data.txt");

auto centers = k_median_pp(
    points,
    /* k = */ 50,
    /* maxIters = */ 20,
    /* maxWeiszfeldIters = */ 50,
    /* tol = */ 1e-4,
    /* seed = */ 42
);
```

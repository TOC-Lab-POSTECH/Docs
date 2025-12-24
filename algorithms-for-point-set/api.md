# API

### Point

#### Member Variable

**`int x, y`**

* Stores integer coordinates of a 2D point.

**`double w`**

* Stores the point weight.
* For original points `P`, typically `w = 1.0`.
* For representative points `R_j`, `w` is set to the number of original points in the corresponding cell (i.e., `|P ∩ cell|`).

***

#### Constructor

**`Point()`**

* Input: None
* Output: Initializes `(x, y) = (0, 0)`, `w = 1.0`
* Complexity: O(1)
* Example:

```cpp
Point p;
```

**`Point(int _x, int _y)`**

* Input: `_x`, `_y`
* Output: Initializes `(x, y) = (_x, _y)`, `w = 1.0`
* Complexity: O(1)
* Example:

```cpp
Point p(10, 20);
```

**`Point(int _x, int _y, double _w)`**

* Input: `_x`, `_y`, `_w`
* Output: Initializes `(x, y, w) = (_x, _y, _w)`
* Complexity: O(1)
* Example:

```cpp
Point rep(5, 5, 37.0);
```

***

#### End-to-End Example

```cpp
Point p0;                 // (0,0), w=1
Point p1(10, 20);         // (10,20), w=1
Point rep(5, 5, 37.0);    // representative point with weight 37
```

***

### Range Counting Oracle: `RangeTree`

#### Type Aliases

```cpp
struct RangeTreeNode {
    int x_left, x_right;
    std::vector<int> ys;
    RangeTreeNode* left  = nullptr;
    RangeTreeNode* right = nullptr;
};
```

***

#### Member Variable

**`RangeTreeNode* root`**

* Root pointer of the range tree.

**`int minX, maxX, minY, maxY`**

* Bounding box of the input point set `P`.

**`vector<Point> sortedPts`**

* Copy of input points sorted by x-coordinate.
* Used only for building the tree.

***

#### Constructor

**`explicit RangeTree(const vector<Point>& pts)`**

* Input: `pts` — the input point set `P`
* Output: Constructs the range tree oracle
* Complexity:
  * Time: O(n log n) for sorting, plus tree construction (node y-list sorting at each node)
  * Space: Stores y-lists in nodes (range-tree style overhead)
* Description:
  * Computes the bounding box of `pts`
  * Sorts points by x
  * Builds a balanced binary tree by splitting in index space
  * Each node stores sorted `ys` for efficient y-range counting
* Example:

```cpp
RangeTree oracle(P);
```

**`~RangeTree()`**

* Input: None
* Output: Frees all nodes
* Complexity: O(#nodes)
* Example:

```cpp
// destructor called automatically
```

***

#### Member Function

**`int range_count(int x1, int x2, int y1, int y2) const`**

* Input: Rectangle `[x1, x2] × [y1, y2]` (inclusive bounds in the implementation)
* Output: Number of points inside the rectangle
* Complexity (typical): O(log² n)
* Description:
  * Normalizes `(x1, x2)` and `(y1, y2)` so that `x1 ≤ x2`, `y1 ≤ y2`
  * If a node’s x-interval is fully covered, counts y via binary search on `ys`
  * Otherwise recurses into children
* Example:

```cpp
int cnt = oracle.range_count(0, 100, 0, 100);
```

**`int getMinX() const`, `int getMaxX() const`, `int getMinY() const`, `int getMaxY() const`**

* Input: None
* Output: Bounding box endpoints
* Complexity: O(1)
* Example:

```cpp
int minX = oracle.getMinX();
int maxY = oracle.getMaxY();
```

***

#### End-to-End Example

```cpp
std::vector<Point> P = /* load or generate, w=1 */;
RangeTree oracle(P);

int cnt = oracle.range_count(0, 100, 0, 100);

std::cout << "bbox x=[" << oracle.getMinX() << "," << oracle.getMaxX()
          << "], y=[" << oracle.getMinY() << "," << oracle.getMaxY() << "]\n";
```

***

### Representative Set Builder: `RepresentativeSetBuilder`

This module builds a weighted representative set (R\_j) using only range-count queries.

#### Type Aliases

```cpp
// No 'using' aliases.
// Quadtree cell:
struct QuadCell {
    int x0, y0;   // bottom-left corner (inclusive)
    int size;     // side length
    int level;    // such that size = 2^level
};

// Configuration for building R_j:
struct RjConfig {
    int n;         // |P|
    int k;         // number of centers
    double eps;    // epsilon parameter
    double rj;     // OPT guess r_j
};
```

***

#### Member Variable

**`const RangeTree& tree`**

* Range counting oracle used to compute `|P ∩ cell|`.

**`RjConfig cfg`**

* Stores `(n, k, eps, rj)` used in the build.

**`QuadCell rootCell`**

* A power-of-two square cell that covers the oracle bounding box.

**`double delta_kmed`**

* Parameter used in the dense/sparse threshold.
* Set in the constructor as:
  * `delta_kmed = 100 * (k * log2(n)) / eps^3`

**`int Kj_limit`**

* Abort threshold for the number of sparse cells `|K_j|`.
* Set in the constructor as:
  * `Kj_limit = floor(4 * k * log2(n) / eps^3) + 10`

**`bool aborted`**

* Set to `true` if `|K_j|` grows beyond `Kj_limit`.

**`vector<QuadCell> Kj`**

* The sparse-cell set `K_j`.
* Final representative set `R_j` is built from `K_j`.

***

#### Constructor

**`RepresentativeSetBuilder(const RangeTree& tree, const RjConfig& cfg)`**

* Input:
  * `tree`: range counting oracle
  * `cfg`: configuration `(n, k, eps, rj)`
* Output: Constructs the builder
* Complexity: O(1)
* Description:
  * Reads bounding box from the oracle
  * Constructs the smallest power-of-two square covering it
  * Initializes `delta_kmed` and `Kj_limit`
* Example:

```cpp
RepresentativeSetBuilder builder(oracle, cfg);
```

***

#### Member Function

**`pair<bool, vector<Point>> build()`**

* Input: None
* Output:
  * `(true, Rj)` if the run finishes without aborting
  * `(false, empty)` if aborted because `|K_j|` exceeded `Kj_limit`
* Complexity (oracle-centric):
  * One `range_count` per visited cell during recursion
  * Additional `range_count` calls when converting `K_j` into `R_j` (to assign weights)
* Description:
  1. Clears `Kj` and resets `aborted`
  2. Calls `processCell(rootCell)`
  3. If aborted, returns failure
  4. Otherwise converts each `QuadCell c ∈ Kj` into one representative point:
     * Representative point position: rounded cell center
     * Representative point weight: `w = |P ∩ c|`
* Example:

```cpp
auto [ok, Rj] = builder.build();
```

**`int countPointsInCell(const QuadCell& c) const`**

* Input: cell `c`
* Output: `|P ∩ c|`
* Complexity: one `range_count` call
* Description:
  * Queries the square region `[x0, x0+size-1] × [y0, y0+size-1]`

**`void processCell(const QuadCell& c)`**

* Input: cell `c`
* Output: None (updates `Kj` / `aborted`)
* Description:
  * Computes `n_c = |P ∩ c|`
  * If empty, returns
  * Computes threshold:
    * `threshold = delta_kmed * rj / c.size`
  * If `n_c < threshold` or `c.size == 1`, the cell is treated as **sparse** and added to `Kj`
  * Otherwise the cell is **dense** and split into 4 children (each of size `c.size/2`) and recursed

***

#### End-to-End Example

```cpp
RangeTree oracle(P);

RjConfig cfg;
cfg.n   = (int)P.size();
cfg.k   = 100;
cfg.eps = 0.2;
cfg.rj  = 1e6;  // OPT guess

RepresentativeSetBuilder builder(oracle, cfg);
auto [ok, Rj] = builder.build();

if (!ok) {
    std::cout << "R_j build aborted (|K_j| too large).\n";
} else {
    std::cout << "R_j built: |R_j|=" << Rj.size() << "\n";
    // Each point in Rj has weight w = |P ∩ cell|
}
```

***

### Weighted k-median++ (functions)

This module provides the solver run on either original `P` (weights 1) or representative `R_j` (weights = cell counts).

***

#### Member Function

**`double euclidean_distance(const Point& a, const Point& b)`**

* Input: two points
* Output: Euclidean distance
* Complexity: O(1)
* Example:

```cpp
double d = euclidean_distance(p, c);
```

**`double k_median_cost(const vector<Point>& pts, const vector<Point>& centers)`**

* Input:
  * `pts`: weighted points
  * `centers`: center points
* Output: weighted k-median cost `sum_p w(p) * min_c dist(p,c)`
* Complexity: O(|pts| · k)
* Example:

```cpp
double cost = k_median_cost(Rj, centers);
```

**`Point geometric_median_weiszfeld(const vector<Point>& clusterPoints, double x0, double y0, int maxIter=100, double tol=1e-3)`**

* Input:
  * `clusterPoints`: points assigned to a cluster (weighted)
  * `(x0, y0)`: initial position (double)
  * `maxIter`, `tol`: iteration/tolerance parameters
* Output: approximate weighted geometric median (rounded to integer coords)
* Complexity: O(maxIter · |clusterPoints|)
* Example:

```cpp
Point m = geometric_median_weiszfeld(clusterPts, initX, initY, 50, 1e-3);
```

**`vector<Point> k_median_pp_init(const vector<Point>& pts, int k, mt19937& rng)`**

* Input:
  * `pts`: weighted points
  * `k`: number of centers
  * `rng`: random generator
* Output: initial centers
* Complexity: depends on k; each new center scans points to compute nearest distance
* Description:
  * First center is sampled proportional to weight `w(p)`
  * Each additional center is sampled proportional to `w(p) * D(p)` where `D(p)` is distance to nearest chosen center
* Example:

```cpp
auto initCenters = k_median_pp_init(Rj, k, rng);
```

**`vector<Point> k_median_pp(const vector<Point>& pts, int k, int maxIters=100, int maxWeiszfeldIters=50, double tol=1e-3, uint32_t seed=12345)`**

* Input:
  * `pts`: weighted points
  * `k`: number of centers
  * `maxIters`: Lloyd iterations
  * `maxWeiszfeldIters`: iterations inside each 1-median update
  * `tol`: tolerance
  * `seed`: RNG seed
* Output: final centers
* Complexity (high-level):
  * Each Lloyd iteration: assignment O(|pts|·k), plus per-cluster Weiszfeld updates
* Example:

```cpp
auto centers = k_median_pp(Rj, k, 20, 50, 1e-3, 12345);
```

***

#### End-to-End Example

```cpp
int k = 100;

// Run solver on the representative set
auto centers = k_median_pp(Rj, k, 20, 50, 1e-3, 12345);

// Evaluate on Rj and on original P
double costRj = k_median_cost(Rj, centers);
double costP  = k_median_cost(P,  centers);

std::cout << "cost(Rj)=" << costRj << ", cost(P)=" << costP << "\n";
```

***

### End-to-End Pipeline (OPT-guess loop)

This is the full algorithm

#### End-to-End Example

```cpp
#include "point.h"
#include "RangeCountingOracle.h"
#include "Representative.h"
#include "k_median++.h"

int main() {
    // 1) Build original point set P
    std::vector<Point> P = /* generate or load, w=1 */;

    // 2) Build Range Counting Oracle
    RangeTree tree(P);

    // 3) Guess loop parameters
    const int N = (int)P.size();
    const int k = 100;
    const double eps_rep   = 0.2;
    const double eps_guess = 0.2;

    double OPT_min = 1.0;
    double OPT_max = 4.0 * N * N;
    int t_max = (int)std::ceil(std::log(OPT_max / OPT_min) / std::log(1.0 + eps_guess));

    double bestCostOnP = std::numeric_limits<double>::infinity();
    std::vector<Point> bestCenters;
    int bestJ = -1;

    // 4) For each guess r_j, build R_j and solve k-median on R_j
    for (int j = 0; j <= t_max; ++j) {
        double rj = OPT_min * std::pow(1.0 + eps_guess, j);

        RjConfig cfg;
        cfg.n   = N;
        cfg.k   = k;
        cfg.eps = eps_rep;
        cfg.rj  = rj;

        RepresentativeSetBuilder builder(tree, cfg);
        auto [ok, Rj] = builder.build();
        if (!ok) continue; // aborted (|K_j| too large)

        auto centers = k_median_pp(Rj, k, 20, 50, 1e-3, (uint32_t)(2000 + j));
        double costOnP = k_median_cost(P, centers);

        if (costOnP < bestCostOnP) {
            bestCostOnP = costOnP;
            bestCenters = centers;
            bestJ = j;
        }
    }

    // bestCenters is the selected solution
    return 0;
}
```

***

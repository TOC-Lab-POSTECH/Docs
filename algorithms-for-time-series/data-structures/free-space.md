# Free Space

Most algorithms for the continuous [_Fréchet Distance_ (FD)](../algorithms/frechet-distance.md) are built around the **free-space diagram** between two polygonal curves. &#x20;

Given two curves $$P$$ and $$Q$$ and a threshold $$\varepsilon > 0$$, the free-space diagram marks all parameter pairs $$(s, t)$$ such that the distance between $$P(s)$$ and $$Q(t)$$ is at most $$\varepsilon$$. &#x20;

The `FreeSpace` class encapsulates the construction of this diagram **on the boundary of each cell** for a fixed threshold `epsilon`. &#x20;

It takes:

* a polygonal curve `P`,
* a polygonal curve `Q`, and
* an `epsilon` value, and computes the free intervals along:
* the **vertical boundaries** (stored in `L`), associated with segments of `P` against vertices of `Q`
* the **horizontal boundaries** (stored in `B`), associated with segments of `Q` against vertices of `P`.

These boundary intervals are used as building blocks in the decision procedure for the [_Fréchet Distance_ (FD)](../algorithms/frechet-distance.md).

***

### Repository

{% embed url="https://github.com/TOC-Lab-POSTECH/prior_similarity_curves" %}

***

### Class Synopsis

```cpp
typedef std::pair<Point_2, Point_2> PointPair;
typedef std::vector<PointPair> PointPairVector;

class FreeSpace {
 public:
  // Constructor to initialize with two polygonal curves and an epsilon value
  FreeSpace(const PolygonalCurve& P, const PolygonalCurve& Q, double epsilon);

  // Destructor
  ~FreeSpace();

  // Getter
  const PolygonalCurve& getCurveP() const;
  const PolygonalCurve& getCurveQ() const;
  double getEpsilon() const;
  const PointPairVector& getL() const;
  const PointPairVector& getB() const;
  
  // Setter
  void setEpsilon(double newEpsilon);
  
  // Recomputes free-space boundary intervals for current epsilon
  void computeFreeSpace();

 private:
  PolygonalCurve P;  // Polygonal curve P
  PolygonalCurve Q;  // Polygonal curve Q
  double epsilon;    // Epsilon value
  PointPairVector L;  // Free-space intervals on vertical boundaries (segments of P vs. vertices of Q)
  PointPairVector B;  // Free-space intervals on horizontal boundaries (segments of Q vs. vertices of P)

  void processCurveForL();
  void processCurveForB();
  std::pair<int, std::vector<float>> checkPointsOnEdge(
      const Point_2& start, const Point_2& end, const Point_2& point) const;
};

```

***

### Dependencies

#### CGAL

* `Point_2` (via `polygonal_curve.h`)

#### STL

* `<utility>` for `std::pair`
* `<vector>` for `std::vector`
* `<cmath>` for `std::sqrt`, `std::fabs`
* `<stdexcept>`
* Other headers required transitively by `PolygonalCurve`

#### Library

* `"polygonal_curve.h"` – for `PolygonalCurve` and `Point_2`&#x20;

***

### Design Notes

* **Ownership of curves**  &#x20;
  * &#x20;`FreeSpace` stores **copies** of the input `PolygonalCurve` objects `P` and `Q` by value.  &#x20;
  * Modifying the original curves after construction does **not** affect the stored curves inside `FreeSpace`.
* **Free-space representation**  &#x20;
  * &#x20;`L` stores free intervals on **vertical cell boundaries**:    &#x20;
    * For each segment of `P` and each vertex of `Q`, it records the subsegment (in parameter space) of that edge that lies exactly at distance `epsilon` from the vertex. &#x20;
  * `B` stores free intervals on **horizontal cell boundaries**:    &#x20;
    * For each segment of `Q` and each vertex of `P`, it records the analogous subsegment.  &#x20;
  * Each entry in `L` or `B` is a `PointPair`:    &#x20;
    * The two `Point_2` values represent parameter coordinates in the free-space diagram, not geometric coordinates in $$\mathbb{R}^2$$. &#x20;
    * &#x20;A sentinel value `(-1, -1)` in both components indicates **no free interval** for that edge–vertex pair.
* **Epsilon updates**  &#x20;
  * Calling `setEpsilon(newEpsilon)` updates the internal `epsilon` and **recomputes** `L` and `B` immediately. &#x20;
  * &#x20;`computeFreeSpace()` can also be called directly to recompute free space after any manual modifications.
* **Degenerate segments**  &#x20;
  * If a segment has zero length (`start == end`), it is treated as **degenerate** and contributes no free-space interval on its boundary.
* **Complexity**  &#x20;
  * Let `p = P.numPoints()` and `q = Q.numPoints()`.  &#x20;
  * &#x20;`computeFreeSpace()` runs in time $$Θ(p · q)$$ and stores:    &#x20;
    * &#x20;`(p - 1) · q` entries in `L`,    &#x20;
    * `(q - 1) · p` entries in `B`.

***

### Overall Structure

{% tabs %}
{% tab title="Member Variable" %}
<table><thead><tr><th width="100">Name</th><th width="198">Type</th><th>Description</th></tr></thead><tbody><tr><td><code>P</code></td><td><code>PolygonalCurve</code></td><td>Copy of the first polygonal curve (row curve in the free-space diagram).</td></tr><tr><td><code>Q</code></td><td><code>PolygonalCurve</code></td><td>Copy of the second polygonal curve (column curve in the free-space diagram).</td></tr><tr><td><code>epsilon</code></td><td><code>double</code></td><td>Distance threshold that defines the free space (<span class="math">\|p − q\| ≤</span> <code>epsilon</code>).</td></tr><tr><td><code>L</code></td><td><code>PointPairVector</code></td><td>Free-space intervals on vertical boundaries (segments of <code>P</code> vs. vertices of <code>Q</code>).</td></tr><tr><td><code>B</code></td><td><code>PointPairVector</code></td><td>Free-space intervals on horizontal boundaries (segments of <code>Q</code> vs. vertices of <code>P</code>).</td></tr></tbody></table>
{% endtab %}

{% tab title="Constructor" %}
<table><thead><tr><th width="381">Name</th><th>Input</th><th>Complexity</th></tr></thead><tbody><tr><td><code>FreeSpace(const PolygonalCurve&#x26; P, const PolygonalCurve&#x26; Q, double epsilon)</code></td><td>two polygonal curves and a non-negative threshold</td><td><span class="math">Θ(p · q)</span></td></tr></tbody></table>

Here `p = P.numPoints()` and `q = Q.numPoints()`.\
The constructor copies the curves and immediately computes the free-space boundary intervals.
{% endtab %}

{% tab title="Member Function" %}
<table><thead><tr><th width="441">Name</th><th width="135">Input</th><th width="131">Output</th><th>Complexity</th></tr></thead><tbody><tr><td><code>const PolygonalCurve&#x26; getCurveP() const</code></td><td>—</td><td>const ref to first curve</td><td><span class="math">O(1)</span></td></tr><tr><td><code>const PolygonalCurve&#x26; getCurveQ() const</code></td><td>—</td><td>const ref to second curve</td><td><span class="math">O(1)</span></td></tr><tr><td><code>double getEpsilon() const</code></td><td>—</td><td>current threshold</td><td><span class="math">O(1)</span></td></tr><tr><td><code>const PointPairVector&#x26; getL() const</code></td><td>—</td><td>vertical boundary intervals</td><td><span class="math">O(1)</span></td></tr><tr><td><code>const PointPairVector&#x26; getB() const</code></td><td>—</td><td>horizontal boundary intervals</td><td><span class="math">O(1)</span></td></tr><tr><td><code>void setEpsilon(double newEpsilon)</code></td><td>new threshold</td><td>—</td><td><span class="math">Θ(p · q)</span></td></tr><tr><td><code>void computeFreeSpace()</code></td><td>—</td><td>—</td><td><span class="math">Θ(p · q)</span></td></tr></tbody></table>
{% endtab %}
{% endtabs %}

***

### Example Usage

```cpp
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include "polygonal_curve.h"
#include "free_space.h"

#include <iostream>
#include <vector>

using Kernel  = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_2 = Kernel::Point_2;

int main() {
  // Define two simple polygonal curves P and Q
  std::vector<Point_2> pointsP = {
      Point_2(0.0, 0.0),
      Point_2(1.0, 0.0),
      Point_2(2.0, 0.0)
  };

  std::vector<Point_2> pointsQ = {
      Point_2(0.0, 0.0),
      Point_2(0.0, 1.0),
      Point_2(0.0, 2.0)
  };


  PolygonalCurve P(pointsP);
  PolygonalCurve Q(pointsQ);

  double epsilon = 1.0;

  // Construct free-space object; free-space intervals are computed in the constructor
  FreeSpace fs(P, Q, epsilon);
  
  std::cout << "epsilon = " << fs.getEpsilon() << "\n";
  std::cout << "#L entries: " << fs.getL().size() << "\n";
  std::cout << "#B entries: " << fs.getB().size() << "\n";

  // Access some intervals in L and B
  const PointPairVector& L = fs.getL();
  const PointPairVector& B = fs.getB();
  if (!L.empty()) {
    const PointPair& interval = L[0];
    std::cout << "First L interval: ("
              << interval.first.x() << ", " << interval.first.y() << ") -> ("
              << interval.second.x() << ", " << interval.second.y() << ")\n";
  }
  // Update epsilon and recompute
  fs.setEpsilon(0.5);
  std::cout << "Updated epsilon = " << fs.getEpsilon() << "\n";

  return 0;
}

```

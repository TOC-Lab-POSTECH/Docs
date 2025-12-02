# Critical Value

In the classical Alt–Godau framework for the continuous Fréchet distance, the _decision problem_ is typically wrapped in a search over a finite set of **critical values**. &#x20;

These critical values are distances at which the combinatorial structure of the free-space diagram can change, and therefore they are the only candidates that can be the exact Fréchet distance.

For two polygonal curves $$P$$ and $$Q$$, the critical values fall into three families:

* **Type A** – distances between the endpoints of the curves $$P[0], Q[0]$$ and $$P[p-1], Q[q-1]$$,
* **Type B** – distances between a vertex of one curve and the closest point on an edge of the other curve,
* **Type C** – distances from an endpoint of a segment to the intersection of that segment with the _perpendicular bisector_ of a vertex pair on the other curve.

The `CriticalValue` class encapsulates the enumeration of these three families:

* It owns copies of the two polygonal curves `P` and `Q`.
* It computes all Type A, B, and C distances.
* It merges, sorts, and deduplicates them into a single sequence of **critical values** that can be used in a binary search for the Fréchet distance.

***

### Repository

{% embed url="https://github.com/TOC-Lab-POSTECH/prior_similarity_curves" %}

***

### Class Synopsis

```cpp
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include "polygonal_curve.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_2 Point_2;

class CriticalValue {
 public:
  // Constructor to initialize the polygonal curves P and Q
  CriticalValue(const PolygonalCurve& P, const PolygonalCurve& Q);

  // Destructor
  ~CriticalValue();

  // Functions to compute the distances for each type
  void computeTypeA();
  void computeTypeB();
  void computeTypeC();
  void computeAndSortAllTypes();

  // Getters for the computed values
  const std::vector<double>& getTypeAValues() const;
  const std::vector<double>& getTypeBValues() const;
  const std::vector<double>& getTypeCValues() const;
  const std::vector<double>& getCriticalValues() const;

 private:
  PolygonalCurve P;  // Polygonal curve P
  PolygonalCurve Q;  // Polygonal curve Q

  // Vectors to store each type of distance
  std::vector<double> typeAValues;
  std::vector<double> typeBValues;
  std::vector<double> typeCValues;
  std::vector<double> critical_values;

  // Helper functions for repeated calculations
  double distance(const Point_2& p1, const Point_2& p2) const;
  double closestPointOnEdge(const Point_2& p, const Point_2& start,
                            const Point_2& end) const;
  Point_2 findIntersectionWithPerpendicularBisector(const Point_2& p1,
                                                    const Point_2& p2,
                                                    const Point_2& start,
                                                    const Point_2& end) const;
};

```

***

### Dependencies

#### CGAL

* `Exact_predicates_inexact_constructions_kernel` – geometric kernel for 2D points.
* `Point_2` – 2D point type (through the kernel).
* `squared_distance` (from `<CGAL/squared_distance_2.h>`) – exact squared Euclidean distance.

#### STL

* `<vector>` – storage for lists of distances.
* `<algorithm>` – `std::sort` and `std::unique` for processing the value set.
* `<cmath>` – `std::sqrt`, `std::fabs` for numeric computations.
* `<iostream>` – used only in debugging / examples.

#### Project

* `"polygonal_curve.h"` – for the `PolygonalCurve` type and access to vertex coordinates.

***

### Design Notes

* **Ownership of curves**
  * `CriticalValue` stores **copies** of the polygonal curves `P` and `Q`. &#x20;
  * Modifying the original curves after constructing a `CriticalValue` instance does not affect the stored instance.
* **Lazy vs eager computation**
  * The constructor calls `computeAndSortAllTypes()`, so all value arrays are populated immediately.
  * The individual `computeTypeA/B/C()` methods are exposed mainly for clarity and potential reuse; they append to the corresponding vectors without clearing existing values.
  * `computeAndSortAllTypes()` assumes the vectors are initially empty; calling it multiple times will accumulate duplicates.
* **Distance families**
  * _Type A_ contains **exactly two distances**:
    * distance between the first vertices of `P` and `Q`,
    * distance between the last vertices of `P` and `Q`.
  * _Type B_ includes distances from:
    * each vertex of `P` to the closest point on every edge of `Q`,
    * each vertex of `Q` to the closest point on every edge of `P`.
  * &#x20; _Type C_ covers distances from a vertex to a point on an edge where that point is equidistant to two vertices of the other curve (intersection with a perpendicular bisector). Only intersections lying on the finite segment are used.
* **Sentinel values**
  * The method `findIntersectionWithPerpendicularBisector` returns the point `(-1, -1)` when no valid intersection exists on the given edge.
  * Callers check this sentinel and skip such edges.
* **Complexity**
  * Let `p = P.numPoints()` and `q = Q.numPoints()`.
  * `computeTypeA()` is $$O(1)$$.
  * `computeTypeB()` runs in $$Θ(p·(q−1) + q·(p−1)) = Θ(p·q)$$.
  * `computeTypeC()` runs in $$Θ(p²·(q−1) + q²·(p−1))$$, which is cubic in the worst case.
  * `computeAndSortAllTypes()` adds:
    * the costs above, plus
    * sorting and deduplication in $$O(N \log N)$$, where `N` is the total number of candidate values.

***

### Overall Structure

{% tabs %}
{% tab title="Member Variable" %}
<table><thead><tr><th width="100">Name</th><th width="198">Type</th><th>Description</th></tr></thead><tbody><tr><td><code>P</code></td><td><code>PolygonalCurve</code></td><td>Copy of the first polygonal curve (row curve in the Fréchet framework).  </td></tr><tr><td><code>Q</code></td><td><code>PolygonalCurve</code></td><td>Copy of the second polygonal curve (column curve in the Fréchet framework).</td></tr><tr><td><code>typeAValues</code>    </td><td><code>std::vector&#x3C;double></code></td><td>Distances of Type A (endpoints).</td></tr><tr><td><code>typeBValues</code>    </td><td><code>std::vector&#x3C;double></code></td><td>Distances of Type B (vertex-to-edge).                                      </td></tr><tr><td><code>typeCValues</code>    </td><td><code>std::vector&#x3C;double></code></td><td>Distances of Type C (bisector/edge intersections).                        </td></tr><tr><td><code>critical_values</code></td><td><code>std::vector&#x3C;double></code></td><td>All critical values combined, sorted, and deduplicated.                    </td></tr></tbody></table>
{% endtab %}

{% tab title="Constructor" %}
<table><thead><tr><th width="381">Name</th><th>Input</th><th>Complexity</th></tr></thead><tbody><tr><td><code>CriticalValue(const PolygonalCurve&#x26; P, const PolygonalCurve&#x26; Q)</code></td><td>two polygonal curves          </td><td><span class="math">Θ(p²·(q−1) + q²·(p−1)) + O(N \log N)  </span></td></tr></tbody></table>

The constructor copies the curves, computes all Type A/B/C values, merges them into `critical_values`, sorts them, and removes duplicates.
{% endtab %}

{% tab title="Member Function" %}
<table><thead><tr><th width="441">Name</th><th width="135">Input</th><th width="131">Output</th><th>Complexity</th></tr></thead><tbody><tr><td><code>void computeTypeA()</code></td><td>—</td><td>—</td><td><span class="math">O(1)</span></td></tr><tr><td><code>void computeTypeB()</code></td><td>—</td><td>—</td><td><span class="math">Θ(p·q)</span></td></tr><tr><td><code>void computeTypeC()</code></td><td>—</td><td>—</td><td><span class="math">Θ(p²·(q−1) + q²·(p−1))</span></td></tr><tr><td><code>void computeAndSortAllTypes()</code></td><td>—</td><td>—</td><td><span class="math">Θ(p²·(q−1) + q²·(p−1)) + O(N \log N)</span></td></tr><tr><td><code>const std::vector&#x3C;double>&#x26; getTypeAValues() const</code></td><td>—</td><td>list of Type A values  </td><td><span class="math">O(1)</span></td></tr><tr><td><code>const std::vector&#x3C;double>&#x26; getTypeBValues() const</code></td><td>—</td><td>list of Type B values              </td><td><span class="math">O(1)</span></td></tr><tr><td><code>const std::vector&#x3C;double>&#x26; getTypeCValues() const</code></td><td>—</td><td>list of Type C values              </td><td><span class="math">O(1)</span></td></tr><tr><td><code>const std::vector&#x3C;double>&#x26; getCriticalValues() const</code></td><td>—</td><td>all critical values (sorted unique)</td><td><span class="math">O(1)</span></td></tr></tbody></table>
{% endtab %}
{% endtabs %}

***

### Example Usage

```cpp
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include "polygonal_curve.h"
#include "critical_value.h"

#include <iostream>
#include <vector>

using Kernel  = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_2 = Kernel::Point_2;

int main() {
  // Define two simple polygonal curves P and Q
  std::vector<Point_2> pointsP = {
      Point_2(0.0, 0.0),
      Point_2(1.0, 0.0),
      Point_2(2.0, 1.0)
  };

  std::vector<Point_2> pointsQ = {
      Point_2(0.0, 0.0),
      Point_2(0.0, 1.0),
      Point_2(1.0, 2.0)
  };

  PolygonalCurve P(pointsP);
  PolygonalCurve Q(pointsQ);

  // Construct CriticalValue object (computes all families on construction)
  CriticalValue cv(P, Q);

  const auto& A = cv.getTypeAValues();
  const auto& B = cv.getTypeBValues();
  const auto& C = cv.getTypeCValues();
  const auto& all = cv.getCriticalValues();

  std::cout << "#Type A: " << A.size() << "\\n";
  std::cout << "#Type B: " << B.size() << "\\n";
  std::cout << "#Type C: " << C.size() << "\\n";
  std::cout << "Critical values (sorted, unique):\\n";
  for (double d : all) {
    std::cout << "  " << d << "\\n";
  }

  return 0;
}

```

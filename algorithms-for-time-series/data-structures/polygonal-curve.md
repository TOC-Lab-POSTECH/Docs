# Polygonal Curve

`PolygonalCurve` models a polygonal chain in the 2D Euclidean plane, represented as an ordered sequence of CGAL points (`Kernel::Point_2`). It provides basic operations such as length computation, coordinate transformations, and integer-grid projection. This class is a building block for higher-level algorithms such as Fréchet Distance (FD) and Geometric Edit Distance (GED).

## **Class Synopsis**

```cpp
class PolygonalCurve {
 public:
  PolygonalCurve(const std::vector<Point_2>& points);
  PolygonalCurve(const PolygonalCurve& P_);

  void addPoint(const Point_2& point);

  std::size_t numPoints() const;
  const Point_2& getPoint(std::size_t index) const;

  double curveLength() const;

  void printCurve() const;

  void shiftOrigin(const Point_2& newOrigin);
  void scaleGrid(double scalingFactor);
  void floorCoordinates();

 private:
  std::vector<Point_2> m_points;
};
```

***

## Dependencies

* **CGAL**
  * `Exact_predicates_inexact_constructions_kernel` (robust predicates, inexact constructions)
  * `Point_2` (2D point type)
* **STL**
  * `std::vector` for point storage
  * `std::out_of_range`, `std::invalid_argument` for error handling
  * `<cmath>` for `sqrt`, `floor`

***

## Design Notes

* **Order matters**
  * The points are stored in sequence to define the polyline.
* **Mutability**
  * Transformation functions (`shiftOrigin`, `scaleGrid`, `floorCoordinates`) mutate the curve in place.
* **Robustness**
  * Throws `std::out_of_range` on invalid indexing.
  * Throws `std::invalid_argument` if `scaleGrid` is called with factor `0`.
* **Empty curve**
  * `numPoints() == 0` is valid.
  * `curveLength()` returns `0.0` for curves with ≤1 point.

***

## Overall Structure

{% tabs %}
{% tab title="Member Variable" %}
<table><thead><tr><th width="95">Name</th><th width="198">Type</th><th>Description</th></tr></thead><tbody><tr><td><code>m_points</code></td><td><code>std::vector&#x3C;Point_2></code></td><td>Ordered vertices of the polyline (may be empty). Access via public API only.</td></tr></tbody></table>
{% endtab %}

{% tab title="Constructor" %}
<table><thead><tr><th width="381">Name</th><th>Input</th><th>Complexity</th></tr></thead><tbody><tr><td><code>PolygonalCurve(const std::vector&#x3C;Point_2>&#x26; points)</code></td><td>vector of <span class="math">2D</span> points</td><td><span class="math">O(n)</span></td></tr><tr><td><code>PolygonalCurve(const PolygonalCurve&#x26; P_)</code></td><td>existing curve</td><td><span class="math">O(n)</span></td></tr></tbody></table>
{% endtab %}

{% tab title="Member Function" %}
<table><thead><tr><th width="228">Name</th><th width="135">Input</th><th width="131">Output</th><th>Complexity</th></tr></thead><tbody><tr><td><code>void addPoint(const Point_2&#x26; point)</code></td><td>point</td><td>—</td><td>(Amortized) <span class="math">O(1)</span></td></tr><tr><td><code>std::size_t numPoints() const</code></td><td>—</td><td>count</td><td><span class="math">O(1)</span></td></tr><tr><td><code>const Point_2&#x26; getPoint(std::size_t index) const</code></td><td>index</td><td>const ref</td><td><span class="math">O(1)</span></td></tr><tr><td><code>double curveLength() const</code></td><td>—</td><td>length</td><td><span class="math">O(n)</span></td></tr><tr><td><code>void printCurve() const</code></td><td>—</td><td>—</td><td><span class="math">O(n)</span></td></tr><tr><td><code>void shiftOrigin(const Point_2&#x26; newOrigin)</code></td><td>origin point</td><td>—</td><td><span class="math">O(n)</span></td></tr><tr><td><code>void scaleGrid(double scalingFactor)</code></td><td>nonzero scalar</td><td>—</td><td><span class="math">O(n)</span></td></tr><tr><td><code>std::invalid_argument</code></td><td>—</td><td>—</td><td><span class="math">O(n)</span></td></tr></tbody></table>
{% endtab %}
{% endtabs %}

***

## Example Usage

```cpp
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include "polygonal_curve.h"
#include <vector>
#include <iostream>

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_2 = Kernel::Point_2;

int main() {
  // Create curve with 3 points
  PolygonalCurve curve({ {0.0, 0.0}, {2.0, 1.0}, {4.0, 2.0} });

  // Append new point
  curve.addPoint(Point_2(6.0, 1.0));

  std::cout << "Number of points: " << curve.numPoints() << "\n";
  std::cout << "Curve length: " << curve.curveLength() << "\n";

  // Apply transformations
  curve.shiftOrigin(Point_2(1.0, 1.0));
  curve.scaleGrid(2.0);
  curve.floorCoordinates();

  curve.printCurve(); // Debug print
}

```

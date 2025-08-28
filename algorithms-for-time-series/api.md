# API

## Polygonal Curve

### Type Aliases

```cpp
using Kernel  = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_2 = Kernel::Point_2;
```

***

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

### Member Variable

#### `std::vector<Point_2> m_points`&#x20;

* Stores the ordered sequence of vertices defining the polygonal curve.
* **Description**
  * The curve is interpreted as a polyline connecting these points in order.
  * May be empty.
  * Direct access is private; use public methods (`addPoint`, `getPoint`, etc.) to manipulate.

***

### Constructor

#### `PolygonalCurve(const std::vector<Point_2>& points)`

* **Input:** `points` — a vector of 2D points.
* **Output:** A polygonal curve initialized with these points.
* **Complexity**: $$O(n)$$
* **Description**
  * Copies the input vector into `m_points`.
* **Example**

```cpp
std::vector<Point_2> pts = { {0,0}, {1,1}, {2,1} };
PolygonalCurve P(pts);
```

&#x20;

#### `PolygonalCurve(const PolygonalCurve& P_)`

* **Input:** `P_` — an existing polygonal curve.
* **Output:** A deep copy of `P_`.
* **Complexity**: $$O(n)$$
* **Description**
  * Copies all points from `P_` into a new object.
* **Example**

```cpp
PolygonalCurve A({ {0,0}, {1,0} });
PolygonalCurve B(A); // identical copy
```

***

### Member Function

#### `void addPoint(const Point_2& point)`

* **Input:** `point` — a 2D point.
* **Output:** None.
* **Complexity**: (Amortized) $$O(1)$$
* **Description**
  * Appends the point to the end of the curve.
* **Example**

```cpp
PolygonalCurve P({ {0,0} });
P.addPoint(Point_2(1.0, 0.0));
P.addPoint(Point_2(2.0, 1.0));
// Curve now has three vertices
```



#### `std::size_t numPoints() const`

* **Input:** None.
* **Output:** Number of points in the curve.
* **Complexity**: $$O(1)$$
* **Description**
  * Returns the size of `m_points`.
* **Example**

<pre class="language-cpp"><code class="lang-cpp"><strong>PolygonalCurve P({ {0,0}, {1,1} });
</strong><strong>std::size_t n = P.numPoints(); // 2
</strong></code></pre>



#### `const Point_2& getPoint(std::size_t index) const`

* **Input:** `index` — zero-based index.
* **Output:** Immutable reference to the point at that index.
* **Complexity**: $$O(1)$$
* **Description**
  * Provides safe read-only access.
* **Example**

```cpp
PolygonalCurve P({ {0,0}, {1,1} });
const Point_2& p = P.getPoint(1); // (1,1)
```



#### `double curveLength() const`

* **Input:** None.
* **Output:** Sum of Euclidean distances between consecutive points.
* **Complexity**: $$O(n)$$
* **Description**
  * Computes polyline length. Returns 0.0 if fewer than 2 points.
* **Example**

```cpp
PolygonalCurve P({ {0,0}, {3,4}, {6,4} });
double L = P.curveLength(); // 8.0
```



#### `void printCurve() const`

* **Input:** None.
* **Output:** Prints all points to `std::cout` in `(x, y)` format.
* **Complexity**: $$O(n)$$
* **Description**
  * For debugging or inspection.
* **Example**

```cpp
PolygonalCurve P({ {0,0}, {1,2} });
P.printCurve();
// Output:
// (0, 0)
// (1, 2)
```



#### `void shiftOrigin(const Point_2& newOrigin)`

* **Input:** `newOrigin` — a 2D point.
* **Output:** None.
* **Complexity**: $$O(n)$$
* **Description**
  * Subtracts `newOrigin` coordinates from each point. Equivalent to translation by `(-newOrigin.x(), -newOrigin.y())`.
* **Example**

```cpp
PolygonalCurve P({ {2,3}, {4,6} });
P.shiftOrigin(Point_2(1,2));
// Points become (1,1) and (3,4)
```



#### `void scaleGrid(double scalingFactor)`

* **Input:** `scalingFactor` — a nonzero double.
* **Output:** None.
* **Complexity**: $$O(n)$$
* **Description**
  * Divides each coordinate by `scalingFactor`.
  * Throws `std::invalid_argument` if `scalingFactor == 0`.
* **Example**

```cpp
PolygonalCurve P({ {2,2}, {4,6} });
P.scaleGrid(2.0); // -> (1,1), (2,3)
```



#### `void floorCoordinates()`

* **Input:** None.
* **Output:** None.
* **Complexity**: $$O(n)$$
* **Description**
  * Applies `std::floor` to each coordinate. Useful for snapping to integer grid.
* **Example**

```cpp
PolygonalCurve P({ {1.9, -0.2}, {-2.1, 3.99} });
P.floorCoordinates();
// -> (1, -1), (-3, 3)
```

***

### End-to-End Example

```cpp
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include "polygonal_curve.h"
#include <iostream>

using Kernel  = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_2 = Kernel::Point_2;

int main() {
  PolygonalCurve C({ {0,0}, {3,4} });
  std::cout << "Length: " << C.curveLength() << "\n"; // 5

  C.shiftOrigin(Point_2(1,1));
  C.scaleGrid(2.0);
  C.floorCoordinates();

  try {
    const Point_2& p = C.getPoint(0);
    std::cout << "p0: (" << p.x() << "," << p.y() << ")\n";
  } catch (const std::out_of_range& e) {
    std::cerr << e.what() << "\n";
  }
}
```

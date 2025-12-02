# Free Space

### Type Aliases

```cpp
typedef std::pair<Point_2, Point_2> PointPair;
typedef std::vector<PointPair> PointPairVector;
```

```cpp

class FreeSpace {
 public:
  FreeSpace(const PolygonalCurve& P, const PolygonalCurve& Q, double epsilon);
  ~FreeSpace();

  const PolygonalCurve& getCurveP() const;
  const PolygonalCurve& getCurveQ() const;

  double getEpsilon() const;

  const PointPairVector& getL() const;
  const PointPairVector& getB() const;

  void setEpsilon(double newEpsilon);
  void computeFreeSpace();
  
 private:
  PolygonalCurve P;
  PolygonalCurve Q;
  
  double epsilon;

  PointPairVector L;
  PointPairVector B;

  void processCurveForL();
  void processCurveForB();
  
  std::pair<int, std::vector<float>> checkPointsOnEdge(
      const Point_2& start, const Point_2& end, const Point_2& point) const;
};

```

### Member Variables

#### `PolygonalCurve P`

* **Description**    &#x20;
  * Copy of the first polygonal curve. Used as the “row” curve in the free-space diagram.

#### `PolygonalCurve Q`

* **Description**    &#x20;
  * Copy of the second polygonal curve. Used as the “column” curve in the free-space diagram.

#### `double epsilon`

* **Description**    &#x20;
  * Distance threshold defining the free space: a pair of points $$(p, q)$$ is considered free if  $$|p - q| \le \epsilon$$.

#### `PointPairVector L`

* **Description**    &#x20;
  * Vector of free intervals on **vertical boundaries**. &#x20;
  * Each entry corresponds to a pair (segment of `P`, vertex of `Q`) and stores:&#x20;
    * either a valid interval as two `Point_2` parameter points in the free-space diagram,&#x20;
    * or the sentinel pair `(-1, -1) -> (-1, -1)` if there is no free interval.

#### `PointPairVector B`

* **Description**   &#x20;
  * Vector of free intervals on **horizontal boundaries**. &#x20;
  * Each entry corresponds to a pair (segment of `Q`, vertex of `P`) and uses the same convention as `L`.

***

### Constructor

#### `FreeSpace(const PolygonalCurve& P, const PolygonalCurve& Q, double epsilon)`

```cpp
FreeSpace::FreeSpace(const PolygonalCurve& P,
                     const PolygonalCurve& Q,
                     double epsilon);
```

* **Input:** &#x20;
  * &#x20;`P`: first polygonal curve.  &#x20;
  * `Q`: second polygonal curve.  &#x20;
  * &#x20;`epsilon`: non-negative distance threshold for free space.
* **Output:**    &#x20;
  * None (constructs the object).
* **Complexity:**    &#x20;
  * Time: $$Θ(p · q)$$, where `p = P.numPoints()` and `q = Q.numPoints()`, due to `computeFreeSpace()` being called in the constructor.    &#x20;
  * Space: $$O(p · q)$$ to store `L` and `B`.
* **Description:**    &#x20;
  * Initializes the `FreeSpace` object by copying the input curves and setting the distance threshold `epsilon`. &#x20;
  * Immediately computes free-space intervals on both vertical (`L`) and horizontal (`B`) boundaries.
* **Example:**  &#x20;
  * ```cpp
    PolygonalCurve P(pointsP);
    PolygonalCurve Q(pointsQ);
      
    double epsilon = 1.0;

    FreeSpace fs(P, Q, epsilon);
    ```

***

### Member Function (Public)

#### `const PolygonalCurve& getCurveP() const`

* **Input:**   &#x20;
  * None. &#x20;
* **Output:**    &#x20;
  * Const reference to the internally stored curve `P`.
* **Complexity:** $$O(1)$$
* **Description:**    &#x20;
  * Provides read-only access to the first curve used to construct the free-space diagram.
* **Example:**  &#x20;
  * ```cpp
    const PolygonalCurve& P_internal = fs.getCurveP();

    std::cout << "Curve P has " << P_internal.numPoints() << " points\n";
    ```

#### `const PolygonalCurve& getCurveQ() const`

* **Input:**   &#x20;
  * None. &#x20;
* **Output:**   &#x20;
  * Const reference to the internally stored curve `Q`.
* **Complexity:** $$O(1)$$
* **Description:**    &#x20;
  * Provides read-only access to the second curve used to construct the free-space diagram.

#### `double getEpsilon() const`

* **Input:**    &#x20;
  * None.
* **Output:**    &#x20;
  * Current epsilon value.
* **Complexity:** $$O(1)$$
* **Description:**    &#x20;
  * Returns the distance threshold currently used for computing free-space intervals.

#### `const PointPairVector& getL() const`

* **Input:**    &#x20;
  * None.
* **Output:**    &#x20;
  * Const reference to the vector of vertical boundary intervals `L`.
* **Complexity:** $$O(1)$$
* **Description:**  &#x20;
  * Gives access to the free-space intervals on vertical boundaries (segments of `P` vs. vertices of `Q`).
* **Example:**  &#x20;
  * ```cpp
    const PointPairVector& L = fs.getL();
      for (const PointPair& interval : L) {
        // skip empty intervals
        if (interval.first.x() < 0.0) continue;
        // process valid free interval
    }
    ```

#### `const PointPairVector& getB() const`

* **Input:**    &#x20;
  * None.
* **Output:**    &#x20;
  * Const reference to the vector of horizontal boundary intervals `B`.
* **Complexity:** $$O(1)$$
* **Description:**    &#x20;
  * Gives access to the free-space intervals on horizontal boundaries (segments of `Q` vs. vertices of `P`).

#### `void setEpsilon(double newEpsilon)`

* **Input:**  &#x20;
  * `newEpsilon`: new distance threshold.
* **Output:**  &#x20;
  * None. &#x20;
* **Complexity:**    &#x20;
  * Time: $$Θ(p · q)$$, as `computeFreeSpace()` is called internally.    &#x20;
  * Space: unchanged; contents of `L` and `B` are overwritten.
* **Description:**    &#x20;
  * Updates the epsilon value and **recomputes** all free-space boundary intervals. &#x20;
  * Use this when experimenting with different thresholds for the same pair of curves.
* **Example:**  &#x20;
  * ```cpp
    fs.setEpsilon(0.5);

    std::cout << "New epsilon = " << fs.getEpsilon() << "\n";
    ```

#### `void computeFreeSpace()`

* **Input:**    &#x20;
  * None.
* **Output:**    &#x20;
  * None.
* **Complexity:**    &#x20;
  * Time: $$Θ(p · q)$$  &#x20;
  * Space: $$O(p · q)$$ for storing the intervals.
* **Description:**    &#x20;
  * Clears `L` and `B` and recomputes all free-space boundary intervals using the current `epsilon`.&#x20;
  * This is automatically called by the constructor and by `setEpsilon`, but can also be invoked explicitly.

***

### Member Functions (Private)

#### `void processCurveForL()`

* **Input:**    &#x20;
  * Uses internal `P`, `Q`, and `epsilon`.
* **Output:**    &#x20;
  * Populates `L`.
* **Complexity:** $$Θ((p - 1) · q)$$
* **Description:**    &#x20;
  * Iterates over every segment of `P` and every vertex of `Q`, calling `checkPointsOnEdge` to find intersections of the edge with the circle of radius `epsilon` centered at the vertex. &#x20;
  * Converts valid parameter intervals into free-space coordinates and appends them to `L`.     If no intersection exists, appends a sentinel `(-1, -1) -> (-1, -1)` pair.

#### `void processCurveForB()`

* **Input:**    &#x20;
  * Uses internal `P`, `Q`, and `epsilon`.
* **Output:**    &#x20;
  * Populates `B`.
* **Complexity:** $$Θ((q - 1) · p)$$
* **Description:**    &#x20;
  * Symmetric to `processCurveForL`, but iterating over segments of `Q` and vertices of `P`.     The resulting intervals are pushed into `B`.

#### `std::pair<int, std::vector<float>> checkPointsOnEdge(...) const`

* **Input:**  &#x20;
  * &#x20;`start`: start point of the segment.  &#x20;
  * &#x20;`end`: end point of the segment.  &#x20;
  * &#x20;`point`: query point whose `epsilon`-circle is considered.
* **Output:**  &#x20;
  * &#x20;A pair `(count, portions)` where:    &#x20;
    * &#x20;`count` is the number of parameter values $$t$$ in $$[0, 1]$$ such that the point on the segment at parameter `t` is at distance exactly `epsilon` from `point`.      &#x20;
      * Possible values:        &#x20;
        * `0`: no intersection,       &#x20;
        * &#x20;`1`: a single intersection point on the segment,        &#x20;
        * &#x20;`2`: an interval of points on the segment at distance `≤ epsilon`, represented by two boundary parameters.&#x20;
    * &#x20;`portions` is a vector of the corresponding parameter values in $$[0, 1]$$.
* **Complexity:** $$O(1)$$
* **Description:**    &#x20;
  * Treats the segment from `start` to `end` as a straight line and considers the circle of radius `epsilon` around `point`. &#x20;
  * Computes the projection of `point` onto the line, checks whether the closest point is inside the segment, and then solves for the intersection(s) between the circle and the segment.
  * Handles several cases:&#x20;
    * No intersection (`count == 0`),&#x20;
    * Single intersection (`count == 1`),&#x20;
    * An interval of free points (`count == 2`), possibly truncated at the segment endpoints.

***

### End-to-End Example

```cpp
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include "polygonal_curve.h"
#include "free_space.h"

#include <iostream>
#include <vector>

using Kernel  = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_2 = Kernel::Point_2;

int main() {
  std::vector<Point_2> P_points = {
      Point_2(0.0, 0.0),
      Point_2(1.0, 0.0),
      Point_2(2.0, 0.0)
  };

  std::vector<Point_2> Q_points = {
      Point_2(0.0, 0.0),
      Point_2(0.0, 1.0),
      Point_2(0.0, 2.0)
  };

  PolygonalCurve P(P_points);
  PolygonalCurve Q(Q_points);

  FreeSpace fs(P, Q, 1.0);

  std::cout << "Initial epsilon: " << fs.getEpsilon() << "\n";
  std::cout << "L size: " << fs.getL().size() << "\n";
  std::cout << "B size: " << fs.getB().size() << "\n";

  // Change epsilon and recompute free-space
  fs.setEpsilon(0.5);
  std::cout << "Updated epsilon: " << fs.getEpsilon() << "\n";

  return 0;
}

```

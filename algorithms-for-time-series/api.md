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

***

## Free Space

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
```

***

## Decision Problem

### Type Aliases

`DecisionProblem` reuses the free-space data types defined in `free_space.h`:

```cpp
typedef std::pair<Point_2, Point_2> PointPair;
typedef std::vector<PointPair> PointPairVector;
```

```cpp
class DecisionProblem {
 public:
  DecisionProblem(const PolygonalCurve& P, const PolygonalCurve& Q,
                  double epsilon);

  bool doesMonotoneCurveExist() const;
  const PolygonalCurve& getCurveP() const;
  const PolygonalCurve& getCurveQ() const;
  double getEpsilon() const;
  void setEpsilon(double newEpsilon);
  void checkMonotoneCurve();
 private:
  PolygonalCurve P;
  PolygonalCurve Q;
  double epsilon;
  FreeSpace freeSpace;
  PointPairVector L_R;
  PointPairVector B_R;

  bool monotoneCurveExists;

  bool checkStartAndEndConditions();
  bool checkIfMonotoneCurveExists();
};

```

***

### Member Variables

#### `PolygonalCurve P`

* **Description:**
  * Copy of the first polygonal curve. Conceptually forms the “row” axis of the free-space diagram.

#### `PolygonalCurve Q`

* **Description:** &#x20;
  * Copy of the second polygonal curve. Conceptually forms the “column” axis of the free-space diagram.

#### `double epsilon`

* **Description:** &#x20;
  * Current distance threshold $$\varepsilon$$. The decision problem asks whether the Fréchet distance between `P` and `Q` is at most this value.

#### `FreeSpace freeSpace`

* **Description:** &#x20;
  * Stores the free-space boundary intervals for the pair $$(P, Q)$$ at distance `epsilon`. &#x20;
  * Used as the geometric input to the reachability computation.

#### `PointPairVector L_R`

* **Description:** &#x20;
  * Vector of reachable intervals on **vertical boundaries**. &#x20;
  * Each element corresponds to an entry of `freeSpace.getL()` and encodes the subset of that boundary that can be reached from $$(0, 0)$$ by a monotone path. &#x20;
  * A sentinel pair `(-1, -1) -> (-1, -1)` means the boundary is not reachable.

#### `PointPairVector B_R`

* **Description:** &#x20;
  * Vector of reachable intervals on **horizontal boundaries**. &#x20;
  * Symmetric to `L_R`, but for entries of `freeSpace.getB()`.

#### `bool monotoneCurveExists`

* **Description:** &#x20;
  * Cached Boolean result of the decision problem for the current curves `P`, `Q`, and threshold `epsilon`. &#x20;
  * It is updated by `checkMonotoneCurve()` and queried via `doesMonotoneCurveExist()`.

***

### Constructor

#### `DecisionProblem(const PolygonalCurve& P, const PolygonalCurve& Q, double epsilon)`

* **Input:**
  * &#x20;`P`: first polygonal curve.
  * `Q`: second polygonal curve.
  * `epsilon`: non-negative distance threshold for the decision.
* **Output:** &#x20;
  * None (constructs the object).
* **Complexity:** &#x20;
  * Time: $$\Theta(p \cdot q)$$, where `p = P.numPoints()` and `q = Q.numPoints()`. &#x20;
    * Constructs `freeSpace` for the given curves and `epsilon`.
    * Calls `checkMonotoneCurve()`, which performs the reachability computation.
    * Space: $$O(p \cdot q)$$ for storing free-space and reachable intervals.
  * Space: $$O(p \cdot q)$$ for storing free-space and reachable intervals.
* **Description:** &#x20;
  * Initializes the decision problem for curves `P` and `Q` at threshold `epsilon`. &#x20;
  * Immediately computes whether a monotone path exists and stores the result.
* **Example:**
  * ```cpp
    PolygonalCurve P(pointsP);
    PolygonalCurve Q(pointsQ);
    double epsilon = 1.0;

    DecisionProblem dp(P, Q, epsilon);
    if (dp.doesMonotoneCurveExist()) {
      std::cout << "Fréchet distance ≤ epsilon\n";
    }

    ```

***

### Member Functions (Public)

#### `bool doesMonotoneCurveExist() const`

* **Input:** &#x20;
  * None.
* **Output:** &#x20;
  * `true` if a monotone path exists in the free-space diagram from $$(0, 0)$$ to $$(q-1, p-1)$$ for the current `epsilon`; `false` otherwise.
* **Complexity:**  $$O(1)$$
* **Description:** &#x20;
  * Returns the cached decision result computed by the most recent call to `checkMonotoneCurve()` (which is invoked automatically by the constructor and `setEpsilon`).

#### `const PolygonalCurve& getCurveP() const`

* **Input:** &#x20;
  * None.
* **Output:** &#x20;
  * Const reference to the internally stored curve `P`.
* **Complexity:** $$O(1)$$
* **Description:** &#x20;
  * Provides read-only access to the first curve used in the decision problem.

#### `const PolygonalCurve& getCurveQ() const`

* **Input:** &#x20;
  * None.
* **Output:** &#x20;
  * Const reference to the internally stored curve `Q`.
* **Complexity:** $$O(1)$$
* **Description:**&#x20;
  * Provides read-only access to the second curve used in the decision problem.

#### `double getEpsilon() const`

* **Input:** &#x20;
  * None.
* **Output:** &#x20;
  * Current distance threshold `epsilon`.
* **Complexity:**  $$O(1)$$
* **Description:** &#x20;
  * Returns the value of $$\varepsilon$$ for which the current decision result is defined.

#### `void setEpsilon(double newEpsilon)`

* **Input:**
  * `newEpsilon`: new distance threshold.
* **Output:**&#x20;
  * None.
* **Complexity:** &#x20;
  * Time: $$\Theta(p \cdot q)$$
    * Calls `freeSpace.setEpsilon(newEpsilon)` to recompute free-space boundaries.
    * Invokes `checkMonotoneCurve()` to recompute reachability.
  * Space: unchanged; contents of `freeSpace`, `L_R`, and `B_R` are overwritten.
* **Description:** &#x20;
  * Updates the threshold `epsilon` and recomputes the decision result. &#x20;
  * Use this method to test multiple thresholds for the same pair of curves without reconstructing the object.
* **Example:**

```cpp
dp.setEpsilon(0.5);
std::cout << "Now epsilon = " << dp.getEpsilon() << "\n";
std::cout << "Monotone curve exists? "
            << (dp.doesMonotoneCurveExist() ? "yes" : "no") << "\n";
```

#### `void checkMonotoneCurve()`

* **Input:** &#x20;
  * None.
* **Output:** &#x20;
  * None (updates `monotoneCurveExists`, `L_R`, and `B_R`).
* **Complexity:**  $$\Theta(p \cdot q)$$
* **Description:**&#x20;
  * Re-evaluates the decision problem for the current curves and `epsilon`:
    * Calls `checkStartAndEndConditions()` to verify that:
      * $$(0, 0)$$ lies in the first free-space boundary cell, and
      * $$(q-1, p-1)$$ lies in the last free-space boundary cell.
      * If either fails, no monotone path can exist.
      * If either fails, no monotone path can exist.
    * If either fails, no monotone path can exist.
  * If both endpoints are feasible, calls `checkIfMonotoneCurveExists()` to propagate reachability across all boundaries and sets `monotoneCurveExists` accordingly.

***

### Member Functions (Private)

#### `bool checkStartAndEndConditions()`

* **Input:**
  * None (uses internal `P`, `Q`, and `freeSpace`).
* **Output:** &#x20;
  * `true` if the start and end parameter points are present in the free-space boundaries; `false` otherwise.
* **Complexity:**  $$O(1)$$ (only inspects the first and last entries of `freeSpace.getL()` and `freeSpace.getB()`).
* **Description:** &#x20;
  * Checks two necessary conditions for the existence of a monotone path:
    * The point `(0, 0)` must appear in the first entries of `L` or `B`.
    * The point `(q-1, p-1)` must appear in the last entries of `L` or `B`.
  * If either fails, the algorithm can immediately conclude that no monotone curve exists, avoiding the more expensive reachability computation.

#### `bool checkIfMonotoneCurveExists()`

* **Input:** &#x20;
  * None (uses internal `P`, `Q`, `freeSpace`, `L_R`, and `B_R`).
* **Output:** &#x20;
  * `true` if a monotone path exists; `false` otherwise.
* **Complexity:** $$\Theta(p \cdot q)$$
* **Description:**&#x20;
  * Implements a dynamic-programming style propagation of reachability over the free-space boundaries:&#x20;
    * Initializes `L_R` and `B_R` with sentinel intervals and copies the first column/row from `freeSpace.getL()` and `freeSpace.getB()`.
  *   Iterates over all cells $$(i, j)$$ of the free-space diagram and updates:

      * `L_R[l+1]` based on `L_R[l]`, `B_R[b]`, and `freeSpace.getL()[l+1]`.
      * `B_R[b+1]` based on `B_R[b]`, `L_R[b]`, and `freeSpace.getB()[b+1]`.

      Intersections of intervals are computed using `std::max` and `std::min` on coordinates; empty intersections produce sentinel intervals.
  * After processing all cells, checks whether the start and end parameter points are contained in the reachable boundaries; if both are present, a monotone path exists.

***

### End-to-End Example

```cpp
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include "polygonal_curve.h"
#include "free_space.h"
#include "decision_problem.h"

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

  DecisionProblem dp(P, Q, 1.0);
  if (dp.doesMonotoneCurveExist()) {
    std::cout << "There exists a monotone path for epsilon = 1.0\\n";
  } else {
    std::cout << "No monotone path for epsilon = 1.0\\n";
  }

  dp.setEpsilon(0.5);
  if (dp.doesMonotoneCurveExist()) {
    std::cout << "There exists a monotone path for epsilon = 0.5\\n";
  } else {
    std::cout << "No monotone path for epsilon = 0.5\\n";
  }

  return 0;
}

```

***

## Critical Value

### Type Aliases

```cpp
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_2 Point_2;
```

```cpp
class CriticalValue {
 public:
  CriticalValue(const PolygonalCurve& P, const PolygonalCurve& Q);
  ~CriticalValue();

  void computeTypeA();
  void computeTypeB();
  void computeTypeC();
  void computeAndSortAllTypes();

  const std::vector<double>& getTypeAValues() const;
  const std::vector<double>& getTypeBValues() const;
  const std::vector<double>& getTypeCValues() const;
  const std::vector<double>& getCriticalValues() const;

 private:
  PolygonalCurve P;
  PolygonalCurve Q;

  std::vector<double> typeAValues;
  std::vector<double> typeBValues;
  std::vector<double> typeCValues;
  std::vector<double> critical_values;

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

### Member Variables

#### `PolygonalCurve P`

* **Description** &#x20;
  * Copy of the first polygonal curve. It is treated as one of the two curves whose Fréchet distance is to be analyzed.

#### `PolygonalCurve Q`

* **Description** &#x20;
  * Copy of the second polygonal curve.

#### `std::vector<double> typeAValues`

* **Description** &#x20;
  * Stores all distances of **Type A** (endpoints of the curves). &#x20;
  * After a fresh call to `computeTypeA()` in a clean object, this vector has exactly two entries.

#### `std::vector<double> typeBValues`

* **Description** &#x20;
  * Stores all distances of **Type B** (vertex-to-edge distances between `P` and `Q` in both directions).

#### `std::vector<double> typeCValues`

* **Description** &#x20;
  * Stores all distances of **Type C** (distances from vertices to bisector–edge intersection points).

#### `std::vector<double> critical_values`

* **Description** &#x20;
  * Combined list of all critical values from types A, B, and C, after sorting in non-decreasing order and removing duplicates. &#x20;
  * This is the main output to be used by higher-level search algorithms.

***

### Constructor

#### `CriticalValue(const PolygonalCurve& P, const PolygonalCurve& Q)`

* **Input:**
  * `P`: first polygonal curve.
  * `Q`: second polygonal curve.
* **Output:** &#x20;
  * None (constructs the object).
* **Complexity:** &#x20;
  * Time: dominated by `computeAndSortAllTypes()`. &#x20;
    * For `p = P.numPoints()` and `q = Q.numPoints()` this is:
      * $$\Theta(p^2 (q-1) + q^2 (p-1)) + O(N \log N)$$,
    * where `N` is the final size of `critical_values`.
  * Space: $$O(N)$$ to store all distances.
* **Description:** &#x20;
  * Copies the input curves, then immediately computes all Type A/B/C distances and merges them into a sorted, deduplicated vector of critical values.
* **Example:**
  * ```cpp
    PolygonalCurve P(pointsP);
    PolygonalCurve Q(pointsQ);

    CriticalValue cv(P, Q);
    std::cout << "Number of critical values: "
              << cv.getCriticalValues().size() << "\\n";
    ```

***

### Member Functions (Public)

#### `void computeTypeA()`

* **Input:** &#x20;
  * None.
* **Output:** &#x20;
  * Appends Type A distances to `typeAValues`.
* **Complexity:** $$O(1)$$
* **Description:** &#x20;
  * Computes the distances:
    * between the first vertices of `P` and `Q`,
    * between the last vertices of `P` and `Q`, &#x20;

&#x20;             and appends them to `typeAValues`. The function does **not** clear the vector before appending.

#### `void computeTypeB()`

* **Input:** &#x20;
  * None.
* **Output:** &#x20;
  * Appends Type B distances to `typeBValues`.
* **Complexity:** $$Θ(p·q)$$, where `p = P.numPoints()` and `q = Q.numPoints()`.
* **Description:**&#x20;
  * For every vertex of one curve, finds the closest point on each edge of the other curve and appends that distance to `typeBValues`:&#x20;
    * For each vertex of `P`, and for each edge of `Q`,
    * For each vertex of `Q`, and for each edge of `P`.
  * The closest point on an edge is computed via projection onto the segment.

#### `void computeTypeC()`

* **Input:**&#x20;
  * None.
* **Output:** &#x20;
  * Appends Type C distances to `typeCValues`.
* **Complexity:** $$Θ(p²·(q−1) + q²·(p−1))$$
* **Description:** &#x20;
  * Enumerates all ordered pairs of distinct vertices on one curve and all edges on the other curve:
    * For pairs of vertices on `P` and edges of `Q`,
    * For pairs of vertices on `Q` and edges of `P`.
  * For each triple, it:
    * Uses `findIntersectionWithPerpendicularBisector` to find the point on the edge that lies on the perpendicular bisector of the vertex pair.
    * Skips the triple if no such point lies on the edge.
    * Otherwise, computes the distance from one of the vertices to this intersection point and appends it to `typeCValues`.

#### `void computeAndSortAllTypes()`

* **Input:** &#x20;
  * None.
* **Output:** &#x20;
  * Fills `critical_values` with all Type A/B/C distances, sorted and deduplicated.
* **Complexity:** cost of `computeTypeA/B/C()` plus sorting and deduplication:
  * $$\Theta(p^2 (q-1) + q^2 (p-1)) + O(N \log N)$$.
*   **Description:** &#x20;

    * Calls `computeTypeA()`, `computeTypeB()`, and `computeTypeC()` to populate the three type-specific vectors.
    * Concatenates all distances into `critical_values`.
    * Sorts `critical_values` in ascending order.
    * Removes duplicates using `std::unique` and `erase`.

    The method does not clear `typeAValues`, `typeBValues`, or `typeCValues`, so multiple invocations will accumulate data unless the caller clears the vectors manually.

#### `const std::vector<double>& getTypeAValues() const`

* **Input:** &#x20;
  * None.
* **Output:** &#x20;
  * Const reference to the list of Type A distances.
* **Complexity:** $$O(1)$$
* **Description:** &#x20;
  * Provides read-only access to `typeAValues`. &#x20;
  * Intended mainly for debugging and analysis; higher-level algorithms typically use `getCriticalValues()` instead.

#### `const std::vector<double>& getTypeBValues() const`

* **Input:**&#x20;
  * None.
* **Output:** &#x20;
  * Const reference to the list of Type B distances.
* **Complexity:** $$O(1)$$
* **Description:** &#x20;
  * Provides read-only access to `typeBValues`.

#### `const std::vector<double>& getTypeCValues() const`

* **Input:** &#x20;
  * None.
* **Output:** &#x20;
  * Const reference to the list of Type C distances.
* **Complexity:** $$O(1)$$
* **Description:**&#x20;
  * Provides read-only access to `typeCValues`.

#### `const std::vector<double>& getCriticalValues() const`

* **Input:** &#x20;
  * None.
* **Output:** &#x20;
  * Const reference to the sorted, deduplicated list of all critical values.
* **Complexity:** $$O(1)$$
* **Description:**&#x20;
  * Returns the primary output of the `CriticalValue` class.&#x20;
  * This sequence is suitable for binary search or parametric search to find the exact Fréchet distance.

***

### Member Functions (Private)

#### `double distance(const Point_2& p1, const Point_2& p2) const`

* **Input:**
  * `p1`, `p2`: two points in the plane.
* **Output:**&#x20;
  * Euclidean distance $$\|p1 - p2\|$$.
* **Complexity:** $$O(1)$$
* **Description:** &#x20;
  * Wraps `CGAL::squared_distance` and `std::sqrt` to compute Euclidean distance with reasonable numerical robustness.

#### `double closestPointOnEdge(const Point_2& p, const Point_2& start, const Point_2& end) const`

* **Input:**
  * `p`: query point.
  * `start`, `end`: endpoints of the segment.
* **Output:** &#x20;
  * Distance from `p` to the closest point on the segment `[start, end]`.
* **Complexity:** $$O(1)$$
* **Description:**&#x20;
  * Projects `p` onto the supporting line of the segment and clamps the parameter to $$[0, 1]$$to obtain the closest point on the finite segment. &#x20;
  * Handles degenerate segments (where `start == end`) by returning the distance from `p` to `start`.

#### `Point_2 findIntersectionWithPerpendicularBisector(const Point_2& p1, const Point_2& p2, const Point_2& start, const Point_2& end) const`

* **Input:**
  * `p1`, `p2`: two distinct vertices whose perpendicular bisector is considered.
  * `start`, `end`: endpoints of an edge on the other curve.
* **Output:** &#x20;
  * A point on the segment `[start, end]` that lies on the perpendicular bisector of `p1` and `p2`, or the sentinel point `(-1, -1)` if no such point exists on the edge.
* **Complexity:** $$O(1)$$
* **Description:** &#x20;
  * Works in the parameter space along vector $$p_1 p_2$$:
    * Projects `start` and `end` onto vector $$p_1 p_2$$ to obtain scalar parameters `a` and `b`.
    * If the interval `[a, b]` does not contain `0.5`, no intersection lies on the segment.
    * Handles special cases where one or both projections equal `0.5` exactly, returning the corresponding endpoint or the closer one to `p1`.
    * Otherwise, interpolates along the edge `[start, end]` with ratio `(0.5 - a) : (b - 0.5)` to find the intersection point.
  * This function is used in `computeTypeC()` to identify points on edges where distances to `p1` and `p2` coincide.

***

### End-to-End Example

```cpp
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include "polygonal_curve.h"
#include "critical_value.h"

#include <iostream>
#include <vector>

using Kernel  = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_2 = Kernel::Point_2;

int main() {
  std::vector<Point_2> P_points = {
      Point_2(0.0, 0.0),
      Point_2(1.0, 0.0),
      Point_2(2.0, 1.0)
  };

  std::vector<Point_2> Q_points = {
      Point_2(0.0, 0.0),
      Point_2(0.0, 1.0),
      Point_2(1.0, 2.0)
  };

  PolygonalCurve P(P_points);
  PolygonalCurve Q(Q_points);
  
  CriticalValue cv(P, Q);

  std::cout << "Total critical values: "
            << cv.getCriticalValues().size() << "\\n";

  const auto& vals = cv.getCriticalValues();
  for (double d : vals) {
    std::cout << d << "\\n";
  }

  return 0;
}

```

***

## Fréchet distance

### Type Aliases

```cpp
class FDistance {
 public:
  // Constructor to initialize with two polygonal curves
  FDistance(const PolygonalCurve& P, const PolygonalCurve& Q);

  // Getter
  double getFDistance() const;

 private:
  PolygonalCurve P;           // Polygonal curve P
  PolygonalCurve Q;           // Polygonal curve Q
  CriticalValue criticalVal;  // Critical values object
  DecisionProblem decision;   // Decision problem object

  double fDistance;           // Computed F-distance

  // Helper function to perform binary search on critical values
  void computeFDistance();
};

```

***

### Member Variables

#### `PolygonalCurve P`

* **Description** &#x20;
  * Copy of the first polygonal curve. Used as one of the two inputs for the Fréchet distance computation.

#### `PolygonalCurve Q`

* **Description** &#x20;
  * Copy of the second polygonal curve.

#### `CriticalValue criticalVal`

* **Description** &#x20;
  * Helper object that precomputes and stores all candidate critical distances between curves `P` and `Q`. &#x20;
  * Provides a sorted, unique list of candidate thresholds via `getCriticalValues()`.

#### `DecisionProblem decision`

* **Description** &#x20;
  * Decision oracle for the Fréchet distance.&#x20;
  * For a given `epsilon`, `decision.setEpsilon(epsilon)` recomputes the free-space and reachability, and `doesMonotoneCurveExist()` answers whether the Fréchet distance is at most `epsilon`.

#### `double fDistance`

* **Description** &#x20;
  * Final Fréchet distance computed by `computeFDistance()`. &#x20;
    * A value `>= 0.0` is the smallest critical value found by the binary search for which the decision is `true`. &#x20;
    * A value of `-1.0` is used as a sentinel to indicate that no valid distance was found (e.g., no critical values or inconsistent input).

***

### Constructor

#### `FDistance(const PolygonalCurve& P, const PolygonalCurve& Q)`

* **Input:**
  * `P`: first polygonal curve.
  * `Q`: second polygonal curve.
* **Output:** &#x20;
  * None (constructs the object and computes the result).
* **Complexity:** &#x20;
  * Time:&#x20;
    * Dominated by:
      * construction of `criticalVal` and its internal computations,
      * a binary search over `N` critical values, each calling the decision oracle in (\Theta(p \cdot q)). &#x20;
    * Overall: $$\Theta(p^2 (q-1) + q^2 (p-1)) + O(N \log N + \log N \cdot p \cdot q)$$.
  * Space: $$O(N + p + q)$$, where `N` is the number of critical values.
* **Description:** &#x20;
  * Initializes all internal data structures:
    * Stores copies of `P` and `Q`.
    * Constructs `criticalVal(P, Q)`, which computes all critical values.
    * Constructs `decision(P, Q, 0.0)` as a decision oracle with an initial threshold of `0.0`.
    * Calls `computeFDistance()` to compute and store the final Fréchet distance in `fDistance`.
* **Example:**
  * ```cpp
      PolygonalCurve P(pointsP);
      PolygonalCurve Q(pointsQ);

      FDistance fdist(P, Q);
      std::cout << "Fréchet distance: " << fdist.getFDistance() << "\\n";
    ```

***

### Member Functions (Public)

#### `double getFDistance() const`

* **Input:** &#x20;
  * None.
* **Output:** &#x20;
  * The computed Fréchet distance between the two curves, or `-1.0` if no valid distance was found.
* **Complexity:** $$O(1)$$
* **Description:** &#x20;
  * Returns the value stored in `fDistance`. &#x20;
  * This method does **not** perform any additional computation; all work is done in the constructor and `computeFDistance()`.

***

### Member Functions (Private)

#### `void computeFDistance()`

* **Input:** &#x20;
  * None (uses `criticalVal` and `decision`).
* **Output:** &#x20;
  * Updates `fDistance`.
* **Complexity:**  $$O(\log N \cdot p \cdot q)$$ where:
  * `N` is the number of critical values,
  * `p = P.numPoints()`, `q = Q.numPoints()`.
* **Description:** &#x20;
  * Implements a binary search over the critical values:
    * Retrieves the sorted list of critical values:
      * ```cpp
        const std::vector<double>& criticalValues = criticalVal.getCriticalValues();
        ```
      * If the list is empty, sets `fDistance = -1.0` and returns.
    * Uses indices `left = 0`, `right = criticalValues.size() - 1` and `result = -1.0`:
      * Computes `mid = left + (right - left) / 2`.
      * Sets `double currentEpsilon = criticalValues[mid];`.
      * Calls `decision.setEpsilon(currentEpsilon)` followed by `decision.doesMonotoneCurveExist()`.
    * Update rules:
      * If the decision is `true`, record `result = currentEpsilon` and continue searching in the **left half** (`right = mid - 1`) for a smaller feasible threshold.
      * If the decision is `false`, search the **right half** (`left = mid + 1`) for a larger threshold.
    * After the loop:
      * If `result != -1.0`, assign `fDistance = result`.
      * Otherwise, keep `fDistance = -1.0` as a sentinel.
  * This method encapsulates the entire numeric search phase of the Alt–Godau-style algorithm.

***

### End-to-End Example

```cpp
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include "polygonal_curve.h"
#include "fdistance.h"

#include <iostream>
#include <vector>

using Kernel  = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_2 = Kernel::Point_2;

int main() {
  std::vector<Point_2> P_points = {
      Point_2(0.0, 0.0),
      Point_2(1.0, 0.5),
      Point_2(2.0, 0.0)
  };

  std::vector<Point_2> Q_points = {
      Point_2(0.0, 0.0),
      Point_2(1.0, 1.0),
      Point_2(2.0, 0.0)
  };

  PolygonalCurve P(P_points);
  PolygonalCurve Q(Q_points);

  FDistance fdist(P, Q);

  double d = fdist.getFDistance();
  if (d < 0.0) {
    std::cout << "No valid Fréchet distance found (invalid input or empty curves)." << std::endl;
  } else {
    std::cout << "Fréchet distance between P and Q: " << d << std::endl;
  }

  return 0;
}

```

***

## Geometric Edit Distance

### Type Aliases

#### `CurveAlphabet`

* **Description** &#x20;
  * Symbol representing a grid cell after random shift and scaling. &#x20;
  * The first component is the floored x-coordinate, the second is the floored y-coordinate.

#### `CurveString`

* **Description** &#x20;
  * Sequence of `CurveAlphabet` symbols encoding a polygonal curve as a string.

#### `CurveStringPair`

* **Description** &#x20;
  * Holds the pair of encoded strings corresponding to two input curves.

#### `Matching`

* **Description** &#x20;
  * Monotone matching as a list of index pairs $$(i,j)$$, where `i` indexes into the first curve/string and `j` into the second.

***

### Function

#### `double computeSquareRootApproxGED(const PolygonalCurve& P, const PolygonalCurve& Q)`

* **Input:**
  * `P`: first polygonal curve.
  * `Q`: second polygonal curve.
* **Output:** &#x20;
  * An $$O(\sqrt{n})$$-approximation of $$\mathrm{GED}(P, Q)$$, where $$n = \min(P.\text{numPoints}(), Q.\text{numPoints}())$$.
* **Complexity:** &#x20;
  * Current implementation: $$O(n^2)$$ time, $$O(n^2)$$ space (dominated by repeated SED calls and DP tables).
* **Description:** &#x20;
  * Implements the randomized $$O(\sqrt{n})$$-approximation algorithm:
    * If the sum of Euclidean distances between corresponding points is ≤ 1, use the diagonal matching and return its cost.
    * Otherwise, for geometrically increasing grid parameters `g` and multiple random shifts:
      * Call `transformCurvesToStrings` to encode the curves as strings.
      * Run `SED` with threshold `12*sqrt(n) + 2*g`.
      * If `SED` returns a non-empty matching, compute and return its GED cost via `computeCost`.
    * If no matching is found, return the cost of the empty matching.
  * Because of randomization, calling this function multiple times may yield slightly different results, but the theoretical approximation guarantee remains.
* **Example:** &#x20;
  * See the “Example Usage” section above.

***

### Helper Functions

#### `double computeCost(const PolygonalCurve& P, const PolygonalCurve& Q, const Matching& matching)`

* **Input:**
  * `P`, `Q`: input polygonal curves.
  * `matching`: monotone matching (index pairs) between `P` and `Q`. &#x20;
    * It may be empty, in which case all points are treated as unmatched.
* **Output:** &#x20;
  * GED cost for `matching` using the library’s cost model.
* **Complexity:** &#x20;
* $$O(|M| + m + n)$$, where:
  * $$|M|$$is the size of `matching`,
  * `m = P.numPoints()`, `n = Q.numPoints()`.
* **Description:** &#x20;
  * Computes the cost $$\delta(M) = \sum_{(i,j)\in M} |p_i - q_j|_2 + \left( m - |M| \right) + \left( n - |M| \right)$$, where the last two terms count unmatched points in `P` and `Q`. &#x20;
    * This is the cost function minimized in the definition of GED.

#### `CurveStringPair transformCurvesToStrings(const PolygonalCurve& P, const PolygonalCurve& Q, int g)`

* **Input:**
  * `P`, `Q`: input curves.
  * `g`: grid parameter controlling cell size and the approximation level.
* **Output:** &#x20;
  * Pair of strings `(stringP, stringQ)` encoding `P` and `Q` as sequences of grid-cell symbols.
* **Complexity:** &#x20;
  * $$O(n)$$ in the number of points (assumes `P` and `Q` have comparable length).
*   **Description:** &#x20;

    * Set $$n = \min(P.\text{numPoints}(), Q.\text{numPoints}())$$and $$\delta = g / \sqrt{n}$$.
    * Sample random shift $$(x_0, y_0) \in [0,\delta]^2$$.
    * Copy curves `P` and `Q` into `convertedP` and `convertedQ`.
    * Shift both curves by $$(-x_0, -y_0)$$ using `PolygonalCurve::shiftOrigin`.
    * Scale coordinates by `1/δ` via `PolygonalCurve::scaleGrid(delta)`.
    * Apply `PolygonalCurve::floorCoordinates` to snap points to integer grid.
    * For each point, create a `CurveAlphabet` symbol `(floor(x), floor(y))` and append to the corresponding `CurveString`.



    * The resulting strings provide a discretized, translation-invariant representation of the curves suitable for string edit distance computation.

#### `Matching SED(const CurveString& S, const CurveString& T, double threshold)`

* **Input:**
  * `S`, `T`: encoded curve strings.
  * `threshold`: maximum allowed edit distance.
* **Output:** &#x20;
  * If the string edit distance between `S` and `T` is **≤ threshold**, returns a `Matching` recovered via `backtrace`. &#x20;
  * Otherwise, returns an **empty** matching.
* **Complexity:** &#x20;
  * Worst case $$O(n^2)$$ time and $$O(n^2)$$ space for `n = max(|S|, |T|)`.
*   **Description:** &#x20;

    * Allocates a DP table `D` of size `(n+1) × (m+1)` and initializes first row/column with insertion/deletion costs.
    * Maintains an auxiliary array `L` of length `2k+3` (where `k = floor(threshold)`) to represent the furthest-reaching positions along each diagonal within the band.
    * For each edit distance `e = 0..k` and each diagonal `h` with `|h| ≤ e` and `h ≡ e (mod 2)`:
      * Computes a candidate position `r` from neighboring `L` values.
      * Slides along the diagonal while symbols match, updating `D`.
      * Updates `L[h]` with the new furthest-reaching `r`, or sets it to “infinity” (max int) if out of range.
    * If `D[n][m]` is left at `-1`, the edit distance exceeds `threshold` → return empty matching.
    * Otherwise, call `backtrace(D)` to reconstruct the matching.



    * This function implements a banded dynamic programming SED tailored to the approximated GED setting.

#### `Matching backtrace(const std::vector<std::vector<int>>& D)`

* **Input:**
  * `D`: DP table filled by `SED`. Only the relative values matter; the function assumes that:
    * `D[i][j]` is non-negative for reachable states,
    * `D[n][m]` encodes a finite edit distance.
* **Output:** &#x20;
  * A `Matching` consisting of index pairs `(i, j)` corresponding to diagonal moves in the DP table.
* **Complexity:** &#x20;
  * $$O(m + n)$$, where `m = D.size() - 1` and `n = D[0].size() - 1`.
* **Description:** &#x20;
  * Starting from `(x, y) = (n, m)`, repeatedly:
    * If moving diagonally `(x-1, y-1)` keeps the same edit distance, interpret this as a **match** and append `(x-1, y-1)` to the matching.
    * Otherwise, move left `(x, y-1)` (skip a symbol in `T`) or up `(x-1, y)` (skip a symbol in `S`) depending on which preserves `D`’s consistency.
  * The resulting list of pairs is a **monotone matching** in reverse order; the caller may reverse the vector if needed for forward iteration.

***

### End-to-End Example

```cpp
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include "polygonal_curve.h"
#include "ged.h"

#include <iostream>
#include <vector>

using Kernel  = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_2 = Kernel::Point_2;

int main() {
  std::vector<Point_2> P_points = {
      Point_2(0.0, 0.0),
      Point_2(1.0, 0.2),
      Point_2(2.0, 0.4),
      Point_2(3.0, 0.6)
  };

  std::vector<Point_2> Q_points = {
      Point_2(0.1, 0.0),
      Point_2(0.9, 0.3),
      Point_2(2.1, 0.5),
      Point_2(3.1, 0.7)
  };
  
  PolygonalCurve P(P_points);
  PolygonalCurve Q(Q_points);

  // Main entry point: approximate geometric edit distance
  double gedValue = GED::computeSquareRootApproxGED(P, Q);
  std::cout << "Approximate GED(P, Q) = " << gedValue << std::endl;

  // (Optional) directly use helpers:
  // 1. Encode as strings with some grid parameter g
  int g = 4;
  GED::CurveStringPair strings = GED::transformCurvesToStrings(P, Q, g);

  // 2. Run SED with a threshold and inspect the matching
  double threshold = 12.0 * std::sqrt(static_cast<double>(
                        std::min(P.numPoints(), Q.numPoints())))
                     + 2.0 * g;
  GED::Matching M = GED::SED(strings.first, strings.second, threshold);
  std::cout << "Matching size from SED: " << M.size() << std::endl;

  return 0;
}

```

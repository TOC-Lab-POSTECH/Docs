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

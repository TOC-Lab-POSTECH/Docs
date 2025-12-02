# Decision Problem

In the Alt–Godau algorithm for the continuous Fréchet distance, the **decision problem** asks:

> Given two polygonal curves $$P$$ and $$Q$$ and a threshold $$\varepsilon > 0$$,&#x20;
>
> does there exist a monotone path in the free-space diagram from $$(0, 0)$$ to $$(q-1, p-1)$$?

If such a path exists, then the Fréchet distance between $$P$$ and $$Q$$ is at most $$\varepsilon$$; otherwise it is larger.

The `DecisionProblem` class encapsulates this decision step:

* It owns two polygonal curves `P` and `Q`.
* It maintains a `FreeSpace` object that stores the free-space boundary intervals for a given `epsilon`.
* It runs a dynamic-programming style propagation over these intervals to decide whether a **monotone curve** (i.e., a feasible path) exists.

This class is the core building block used by higher-level routines that perform a binary search over candidate values of $$\varepsilon$$.

***

### Repository

{% embed url="https://github.com/TOC-Lab-POSTECH/prior_similarity_curves" %}

***

### Class Synopsis

```cpp

#include "free_space.h"

class DecisionProblem {
 public:
  // Constructor to initialize with two polygonal curves and epsilon
  DecisionProblem(const PolygonalCurve& P, const PolygonalCurve& Q,
                  double epsilon);

  // Getter
  bool doesMonotoneCurveExist() const;
  const PolygonalCurve& getCurveP() const;
  const PolygonalCurve& getCurveQ() const;
  double getEpsilon() const;

  // Setter
  void setEpsilon(double newEpsilon);

  // Function to check if a monotone curve exists
  void checkMonotoneCurve();

 private:
  PolygonalCurve P;     // Polygonal curve P
  PolygonalCurve Q;     // Polygonal curve Q
  double epsilon;       // Epsilon value
  FreeSpace freeSpace;  // FreeSpace object
  PointPairVector L_R;  // Reachable L (vertical boundaries)
  PointPairVector B_R;  // Reachable B (horizontal boundaries)

  bool monotoneCurveExists;  // True if a monotone curve exists, false otherwise

  // Helper functions for checking the conditions
  bool checkStartAndEndConditions();
  bool checkIfMonotoneCurveExists();
};

```

***

### Dependencies

#### STL

* `<vector>` – for `std::vector` inside `PointPairVector`
* `<algorithm>` – for `std::max` and `std::min` (used in reachability updates)
* `<utility>` – for `std::pair` (via `PointPair`)

#### Library

* `"polygonal_curve.h"` – for `PolygonalCurve` and `Point_2`
* `"free_space.h"` – for `FreeSpace`, `PointPair`, and `PointPairVector`

***

### Design Notes

* **Ownership of curves**
  * `DecisionProblem` stores **copies** of the input curves `P` and `Q` by value.
  * Changing the original curves after constructing a `DecisionProblem` has no effect on the stored curves.
* **Coupling with FreeSpace**
  * The class maintains a `FreeSpace` object that is always kept consistent with the current `epsilon`.
  * The constructor and `setEpsilon` both recompute free-space boundary intervals and re-evaluate the decision.
* **Reachability representation**
  * `L_R` and `B_R` represent **reachable** portions of the free-space diagram:
    * `L_R` mirrors the structure of `freeSpace.getL()` but keeps only the portions that are reachable by a monotone path from $$(0,0)$$.
    * `B_R` mirrors `freeSpace.getB()` for horizontal boundaries.
    * Each entry is a `PointPair`:
  * Each entry is a `PointPair`:
    * The two `Point_2` values encode parameter coordinates in the free-space diagram.
    * The sentinel pair `(-1, -1) -> (-1, -1)` indicates that the corresponding boundary is not reachable.
* **Decision semantics**&#x20;
  * `doesMonotoneCurveExist()` returns the result of the **last call** to `checkMonotoneCurve()`.
  * The constructor and `setEpsilon` call `checkMonotoneCurve()` automatically, so the result is always up-to-date with the current `epsilon`.
* **Complexity**
  * Let `p = P.numPoints()` and `q = Q.numPoints()`.
  * The decision routine runs in time $$\Theta(p \cdot q)$$:
    * Free-space construction via `FreeSpace` is $$\Theta(p \cdot q)$$.
    * The reachability propagation in `checkIfMonotoneCurveExists()` is also $$\Theta(p \cdot q)$$.

***

### Overall Structure

{% tabs %}
{% tab title="Member Variable" %}
<table><thead><tr><th width="100">Name</th><th width="198">Type</th><th>Description</th></tr></thead><tbody><tr><td><code>P</code></td><td><code>PolygonalCurve</code></td><td>Copy of the first polygonal curve (row curve in the free-space diagram).</td></tr><tr><td><code>Q</code></td><td><code>PolygonalCurve</code></td><td>Copy of the second polygonal curve (column curve in the free-space diagram).</td></tr><tr><td><code>epsilon</code></td><td><code>double</code></td><td>Current distance threshold for the decision problem. </td></tr><tr><td><code>freeSpace</code>          </td><td><code>FreeSpace</code>    </td><td>Free-space boundary representation for curves <code>P</code> and <code>Q</code> at distance <code>epsilon</code>.                    </td></tr><tr><td><code>L_R</code>                </td><td><code>PointPairVector</code></td><td>Reachable intervals on vertical boundaries (subset of <code>freeSpace.getL()</code>).                          </td></tr><tr><td><code>B_R</code>                </td><td><code>PointPairVector</code></td><td>Reachable intervals on horizontal boundaries (subset of <code>freeSpace.getB()</code>).                        </td></tr><tr><td><code>monotoneCurveExists</code></td><td><code>bool</code>          </td><td>Cached result of the decision problem: <code>true</code> if a monotone curve exists, <code>false</code> otherwise.        </td></tr></tbody></table>
{% endtab %}

{% tab title="Constructor" %}
<table><thead><tr><th width="381">Name</th><th>Input</th><th>Complexity</th></tr></thead><tbody><tr><td><code>DecisionProblem(const PolygonalCurve&#x26; P, const PolygonalCurve&#x26; Q, double epsilon)</code></td><td>two polygonal curves and a non-negative threshold</td><td><span class="math">Θ(p · q)</span></td></tr></tbody></table>

Here `p = P.numPoints()` and `q = Q.numPoints()`. &#x20;

The constructor copies the curves, constructs the `FreeSpace` object, and immediately evaluates whether a monotone curve exists.
{% endtab %}

{% tab title="Member Function" %}
<table><thead><tr><th width="441">Name</th><th width="135">Input</th><th width="131">Output</th><th>Complexity</th></tr></thead><tbody><tr><td><code>bool doesMonotoneCurveExist() const</code>    </td><td>—</td><td>decision result    </td><td><span class="math">O(1)</span></td></tr><tr><td><code>const PolygonalCurve&#x26; getCurveP() const</code></td><td>—</td><td>const ref to <code>P</code>  </td><td><span class="math">O(1)</span></td></tr><tr><td><code>const PolygonalCurve&#x26; getCurveQ() const</code></td><td>—</td><td>const ref to <code>Q</code>  </td><td><span class="math">O(1)</span></td></tr><tr><td><code>double getEpsilon() const</code>              </td><td>—</td><td>current threshold  </td><td><span class="math">O(1)</span></td></tr><tr><td><code>void setEpsilon(double newEpsilon)</code>      </td><td>new threshold</td><td>—</td><td><span class="math">Θ(p · q)</span></td></tr><tr><td><code>void checkMonotoneCurve()</code></td><td>—</td><td>—</td><td><span class="math">Θ(p · q)</span></td></tr></tbody></table>
{% endtab %}
{% endtabs %}

***

### Example Usage

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

  // Construct the decision problem for (P, Q, epsilon)
  DecisionProblem dp(P, Q, epsilon);
  std::cout << "epsilon = " << dp.getEpsilon() << "\n";
  std::cout << "Monotone curve exists? "
            << (dp.doesMonotoneCurveExist() ? "yes" : "no") << "\n";

  // Update epsilon and recompute
  dp.setEpsilon(0.5);
  std::cout << "Updated epsilon = " << dp.getEpsilon() << "\n";
  std::cout << "Monotone curve exists now?"
            << (dp.doesMonotoneCurveExist() ? "yes" : "no") << "\n";

  return 0;
}

```

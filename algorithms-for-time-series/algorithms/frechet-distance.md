# Fréchet Distance

The **FDistance** class implements the top-level computation of the continuous Fréchet distance between two polygonal curves $$P$$ and $$Q$$.

Conceptually, the algorithm follows the classical **Alt–Godau** framework:

1. **Enumerate critical values** &#x20;

&#x20;  Use `CriticalValue` to generate all candidate distances at which the structure of the free-space diagram can change (Types A, B, and C).

2. **Decision oracle** &#x20;

&#x20;  For a given threshold $$\varepsilon$$, use `DecisionProblem` to decide whether there exists a **monotone path** in the free-space diagram from $$(0, 0)$$ to $$(q-1, p-1)$$. &#x20;

&#x20;  If such a path exists, then $$\text{Fréchet}(P, Q) \le \varepsilon$$.

3. **Binary search on critical values** &#x20;

&#x20;  Sort the critical values once, and then perform a binary search over this discrete set using the decision oracle:

* If the decision is `true` at index `mid`, search the **left half** (smaller $$\varepsilon$$).
* Otherwise, search the **right half** (larger $$\varepsilon$$).

The smallest critical value for which the decision problem is `true` is returned as the **Fréchet distance** of $$P$$ and $$Q$$. &#x20;

If no critical value satisfies the decision oracle, the implementation stores `-1.0` as a sentinel value.

***

### Repository

{% embed url="https://github.com/TOC-Lab-POSTECH/prior_similarity_curves" %}

***

### Class Synopsis

```cpp
#include "critical_value.h"
#include "decision_problem.h"

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

### Dependencies

#### STL

* `<vector>` – implicit through `CriticalValue` and `DecisionProblem`.
* `<algorithm>` – implicit through `CriticalValue` (sorting critical values).

No additional CGAL headers are included directly in `fdistance.h`; geometric operations are handled in dependent classes.

#### Library

* `"polygonal_curve.h"` – definition of `PolygonalCurve`, used for curves `P` and `Q`.
* `"critical_value.h"` – definition of `CriticalValue`, used to generate critical distances.
* `"decision_problem.h"` – definition of `DecisionProblem`, used as the decision oracle.

***

### Design Notes

* **Ownership and copies**
  * `FDistance` stores copies of the input polygonal curves `P` and `Q`.
  * It owns a `CriticalValue` instance and a `DecisionProblem` instance:
    * `criticalVal` is constructed once for `(P, Q)`.
    * `decision` is instantiated with `(P, Q, 0.0)` and then re-used for different values of `epsilon` via `setEpsilon`.
* **Single-pass computation**
  * The constructor calls `computeFDistance()` immediately, so the Fréchet distance is fully computed at construction time.
  * The result is stored in `fDistance` and can be queried multiple times at $$O(1)$$ cost.
* **Binary search over critical values**
  * The vector of critical values is obtained from `CriticalValue::getCriticalValues()` (already sorted, unique).
  * A standard binary search is performed over indices `[0, N - 1]`, where each probe:
    * Sets a candidate threshold `epsilon = criticalValues[mid]`.
    * Calls `decision.setEpsilon(epsilon)` and queries `doesMonotoneCurveExist()`.
* **Sentinel value**
  * When `criticalValues` is empty, or when no threshold passes the decision test (which mathematically should not happen for finite curves), the class sets:
  * `fDistance = -1.0`.
* **Complexity**
  * Let `p = P.numPoints()`, `q = Q.numPoints()`, and `N` be the number of critical values.
  * Constructing `criticalVal` costs:&#x20;
    * $$\Theta(p^2 (q-1) + q^2 (p-1)) + O(N \log N)$$ (dominated by Type C enumeration and sorting).
  * Each decision call (for a fixed $$\varepsilon$$) runs in $$\Theta(p \cdot q)$$.
  * The binary search performs $$O(\log N)$$ decision calls, for a total of:
    * $$O(\log N \cdot p \cdot q)$$.
  * Overall, `FDistance` construction is dominated by:
    * $$\Theta(p^2 (q-1) + q^2 (p-1)) + O(N \log N + \log N \cdot p \cdot q)$$.

***

### Overall Structure

{% tabs %}
{% tab title="Member Variable" %}
<table><thead><tr><th width="100">Name</th><th width="198">Type</th><th>Description</th></tr></thead><tbody><tr><td><code>P</code></td><td><code>PolygonalCurve</code></td><td>Copy of the first polygonal curve.                                </td></tr><tr><td><code>Q</code></td><td><code>PolygonalCurve</code></td><td>Copy of the second polygonal curve.                                </td></tr><tr><td><code>criticalVal</code></td><td><code>CriticalValue</code></td><td>Helper object that enumerates and stores all critical values.      </td></tr><tr><td><code>decision</code>  </td><td><code>DecisionProblem</code></td><td>Decision oracle to test if Fréchet distance is ≤ given threshold.</td></tr><tr><td><code>fDistance</code></td><td><code>double</code>        </td><td>Final computed Fréchet distance; <code>-1.0</code> is used as a sentinel.    </td></tr></tbody></table>
{% endtab %}

{% tab title="Constructor" %}


<table><thead><tr><th width="381">Name</th><th>Input</th><th>Complexity</th></tr></thead><tbody><tr><td><code>FDistance(const PolygonalCurve&#x26; P, const PolygonalCurve&#x26; Q)</code></td><td>two polygonal curves</td><td><span class="math">\Theta(p^2 (q-1) + q^2 (p-1)) + O(N \log N + \log N \cdot p \cdot q)</span></td></tr></tbody></table>

The constructor copies the curves, constructs `criticalVal` and `decision`, and then computes the Fréchet distance with a binary search over the critical values.
{% endtab %}

{% tab title="Member Function" %}
<table><thead><tr><th width="441">Name</th><th width="135">Input</th><th width="131">Output</th><th>Complexity</th></tr></thead><tbody><tr><td><code>double getFDistance() const</code>        </td><td>—</td><td>Fréchet distance      </td><td><span class="math">O(1)</span></td></tr><tr><td><code>void computeFDistance()</code> (private)</td><td>—</td><td>— (updates <code>fDistance</code>)</td><td><span class="math">O(\log N \cdot p \cdot q)</span>  </td></tr></tbody></table>

Here `N` is the number of critical values, `p = P.numPoints()`, and `q = Q.numPoints()`.
{% endtab %}
{% endtabs %}

***

### Example Usage

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

  FDistance fdist(P, Q);
  std::cout << "Fréchet distance: " << fdist.getFDistance() << "\\n";

  return 0;
}

```

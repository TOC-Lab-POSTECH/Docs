# Geometric Edit Distance

The **GED** namespace provides an implementation of an $$O(\sqrt{n})$$-approximation algorithm for the **Geometric Edit Distance (GED)** between two polygonal curves.

***

### Geometric Edit Distance (Conceptual Overview)

Given two point sequences

* P = $$\langle p_1,\dots,p_m\rangle,\qquad$$
* Q = $$\langle q_1,\dots,q_n\rangle,$$

a **monotone matching** is a set of index pairs

* $$M = {(i_1, j_1), \dots, (i_k, j_k)}$$ such that:
  * indices are strictly increasing in both sequences, and &#x20;
  * order is preserved
    * if $$(i,j)$$ and $$(i',j')$$ are in $$M$$ and $$i < i'$$, then $$j < j'$$.

Any point that does not appear in a pair of $$M$$ is **unmatched**. Let $$\Gamma(M)$$ denote this unmatched set.

For a fixed matching $$M$$, the GED **cost function** implemented in this library is:

* $$\delta(M) = \sum_{(i,j)\in M} \operatorname{dist}(p_i, q_j) + |\Gamma(M)|$$
  * The first term sums Euclidean distances between matched points.
  * The second term is a **linear gap penalty**: each unmatched point contributes a cost of 1.

The **Geometric Edit Distance** between (P) and (Q) is then

* $$\mathrm{GED}(P, Q) = \min_M \delta(M)$$,

where the minimum ranges over all monotone matchings $$M$$.

***

### Approximate GED via Random Grid Shift and SED

Computing GED exactly is expensive, so this implementation follows the $$O(\sqrt{n})$$-approximation framework from the literature:

1. **Diagonal sanity check** &#x20;

* Let $$n = \min(\lvert P\rvert, \lvert Q\rvert)$$. &#x20;
* Compute the sum of Euclidean distances between corresponding points $$(p_i, q_i)$$for $$i = 0,\dots,n-1$$. &#x20;
* If this sum is $$\le 1$$, the algorithm uses the **diagonal matching** $$(0,0), (1,1), \dots, (n-1, n-1)$$and returns its cost. &#x20;
  * This covers the case where the two sequences are already very similar.

2. **Random grid shift and string encoding** &#x20;

* For increasing grid parameters $$g = 2^i$$(for $$i = 0,\dots,\lceil\log_2 n\rceil$$) and several random repetitions:
  * Compute $$\delta = g / \sqrt{n}$$.
  * Choose a random shift $$(x_0, y_0) \in [0,\delta]^2$$.
  * Translate both curves by $$(-x_0, -y_0)$$, scale by $$1/\delta$$, and **floor** all coordinates. &#x20;
    * Each point now lies on an integer grid cell.
  * Interpret each grid cell $$(\lfloor x\rfloor, \lfloor y\rfloor)$$ as a **symbol** in an alphabet, and encode each curve as a **CurveString**.

3. **String Edit Distance (SED) with threshold** &#x20;

* For each random shift:
  * Run a banded **string edit distance** algorithm `SED` on the two CurveStrings,
  * With an edit-distance threshold $$k = \lfloor 12\sqrt{n} + 2g \rfloor$$.
  * If the edit distance exceeds the threshold, `SED` returns an empty matching.
  * Otherwise, it backtracks the DP table to recover a monotone matching in the original point sequences.

4. **Cost evaluation and approximation guarantee** &#x20;

* As soon as a non-empty matching is found, its cost is computed by `computeCost` and returned. &#x20;
  * The resulting value is an $$O(\sqrt{n})$$-approximation of the true GED (under the theoretical assumptions of the algorithm).
* If **no** matching is found across all grid scales and random shifts, the algorithm falls back to the cost of the **empty matching**, i.e., only gap penalties.

The current implementation uses a straightforward dynamic programming SED (without suffix-tree acceleration), so the practical time complexity is $$O(n^2)$$for sequences of length $$n$$.

***

### Repository

{% embed url="https://github.com/TOC-Lab-POSTECH/prior_similarity_curves" %}

***

### Namespace Synopsis

```cpp
#include <utility>
#include <vector>

#include "polygonal_curve.h"

typedef std::pair<int, int> CurveAlphabet;       // Alphabet symbol (grid cell)
typedef std::vector<CurveAlphabet> CurveString;  // Encoded curve as a string
typedef std::pair<CurveString, CurveString> CurveStringPair;
typedef std::vector<std::pair<int, int>> Matching;  // Monotone matching

namespace GED {
// Computes O(n^(1/2))-approximation of GED
double computeSquareRootApproxGED(const PolygonalCurve& P,
                                  const PolygonalCurve& Q);

// Computes the GED cost for given matching for two polygonal curves P, Q
double computeCost(const PolygonalCurve& P, const PolygonalCurve& Q,
                   const Matching& matching);

// Transforms the curves into strings by randomly shifted grid
CurveStringPair transformCurvesToStrings(const PolygonalCurve& P,
                                         const PolygonalCurve& Q, int g);

// Computes the String Edit Distance (SED) for GED
Matching SED(const CurveString& S, const CurveString& T, double threshold);

// Backtrace the DP table and return matching
Matching backtrace(const std::vector<std::vector<int>>& D);
}  // namespace GED

```

***

### Dependencies

#### STL&#x20;

* `<vector>` – storage for strings, matchings, and DP tables.
* `<utility>` – `std::pair` for symbols and index pairs.
* `<cmath>` – `sqrt`, `pow`, `log`, `log2`, `ceil`, `floor`.
* `<random>` – `std::random_device`, `std::mt19937`, `std::uniform_real_distribution` for random grid shifts.
* `<limits>` – `std::numeric_limits<int>::max()` used as an “infinity” sentinel.
* `<cstdlib>`, `<ctime>`, `<iostream>` – used internally for debugging/IO or compatibility.

#### Library

* `"polygonal_curve.h"` – representation of input curves.

***

### Design Notes

* **Namespace-based design** &#x20;
  * All GED-related functionality lives in the `GED` namespace rather than a class:
    * `computeSquareRootApproxGED` is the main entry point.
    * `computeCost`, `transformCurvesToStrings`, `SED`, and `backtrace` are helper routines that can be reused if necessary.
* **Type aliases**
  * `CurveAlphabet` encodes a grid cell as an integer pair $$(x,y)$$.
  * `CurveString` is a sequence of `CurveAlphabet` values – the string representation of a curve.
  * `Matching` is a list of index pairs $$(i,j)$$describing a monotone matching.
* **Randomization**
  * `transformCurvesToStrings` uses a fresh random shift $$(x_0,y_0)$$ in every call.
  * Consequently, `computeSquareRootApproxGED` is **randomized**: running it multiple times on the same curves can yield different matchings and slightly different approximate GED values.
  * The theoretical approximation factor holds in expectation (and with high probability under sufficient repetitions).
* **Cost model**
  * The cost function implemented in `computeCost` matches the theoretical definition:
    * Sum of Euclidean distances for matched pairs.
    * Plus a unit penalty for every unmatched point in either sequence.
* **Thresholded SED**
  * `SED` maintains a DP table `D` of size `(n+1) × (m+1)` but only explores a **band** of width `2k+1` around the main diagonal (where `k = floor(threshold)`).
  * If the string edit distance exceeds `threshold`, the bottom-right cell `D[n][m]` remains `-1` and an empty matching is returned.
  * Otherwise, `backtrace` reconstructs a monotone matching from `D`.
* **Complexity (current implementation)**
  * Let $$n = |S|$$, $$m = |T|$$, and assume both are $$O(n)$$.
  * `transformCurvesToStrings`: $$O(n)$$.
  * `computeCost`: $$O(|M| + m + n)$$; in the worst case, $$|M| = O(n)$$.
  * `SED`: worst-case $$O(n^2)$$time and $$O(n^2)$$space for the DP table.
  * `computeSquareRootApproxGED`:
    * Initial diagonal check: $$O(n)$$.
    * Main loops: $$O(\log n)$$ grid levels, $$O(\log n)$$ repetitions per level, each calling `transformCurvesToStrings` and `SED`.
    * Practically dominated by the SED calls → **overall** $$O(n^2)$$ time for the current implementation.
  * A more advanced implementation with preprocessed suffix trees could reduce the SED part to $$O(n(\log n)^2)$$, but this is **not** implemented here.

***

### Overall Structure

{% tabs %}
{% tab title="Type Alias" %}
<table><thead><tr><th width="100">Name</th><th width="198">Type</th><th>Description</th></tr></thead><tbody><tr><td><code>CurveAlphabet</code>  </td><td><code>std::pair&#x3C;int, int></code>              </td><td>Symbol representing a grid cell <span class="math">(x, y)</span> after discretization.</td></tr><tr><td><code>CurveString</code>    </td><td><code>std::vector&#x3C;CurveAlphabet></code>      </td><td>Encoded curve as a string of grid-cell symbols.              </td></tr><tr><td><code>CurveStringPair</code></td><td><code>std::pair&#x3C;CurveString, CurveString></code></td><td>Pair of encoded strings derived from curves <span class="math">P</span> and <span class="math">Q</span>.</td></tr><tr><td><code>Matching</code></td><td><code>std::vector&#x3C;std::pair&#x3C;int, int>></code></td><td>Monotone matching: list of index pairs <span class="math">(i, j)</span>.          </td></tr></tbody></table>
{% endtab %}

{% tab title="Function" %}
<table><thead><tr><th width="381">Name</th><th>Input</th><th>Output</th><th>Time-Complexity</th></tr></thead><tbody><tr><td><code>computeSquareRootApproxGED(const PolygonalCurve&#x26; P, const PolygonalCurve&#x26; Q)</code></td><td>two polygonal curves</td><td>approximate GED value</td><td>(O(n^2))    </td></tr></tbody></table>
{% endtab %}

{% tab title="Helper Function" %}
<table><thead><tr><th width="441">Name</th><th width="135">Input</th><th width="131">Output</th><th>Complexity</th></tr></thead><tbody><tr><td><code>computeCost(const PolygonalCurve&#x26; P, const PolygonalCurve&#x26; Q, const Matching&#x26; matching)</code></td><td>two curves, monotone matching</td><td>GED cost for this matching</td><td><span class="math">O(|M| + m + n)</span></td></tr><tr><td><code>transformCurvesToStrings(const PolygonalCurve&#x26; P, const PolygonalCurve&#x26; Q, int g)</code></td><td>two curves, grid parameter <code>g</code>  </td><td><code>CurveStringPair</code>          </td><td><span class="math">O(n)</span>  </td></tr><tr><td><code>SED(const CurveString&#x26; S, const CurveString&#x26; T, double threshold)</code></td><td>two encoded strings, threshold  </td><td><code>Matching</code> (possibly empty)</td><td><span class="math">O(n^2)</span></td></tr><tr><td><code>backtrace(const std::vector&#x3C;std::vector&#x3C;int>>&#x26; D)</code></td><td>DP table <code>D</code> from <code>SED</code>        </td><td><code>Matching</code>                </td><td><span class="math">O(m+n)\</span></td></tr></tbody></table>

Here $$n = \min(P.\text{numPoints}(), Q.\text{numPoints}())$$, $$m$$ and $$n$$ denote string lengths, and $$|M|$$ is the size of the matching.
{% endtab %}
{% endtabs %}

***

### Example Usage

```cpp
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include "polygonal_curve.h"
#include "ged.h"

#include <iostream>
#include <vector>

using Kernel  = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_2 = Kernel::Point_2;

int main() {
  // Define two polygonal curves P and Q
  std::vector<Point_2> P_points = {
      Point_2(0.0, 0.0),
      Point_2(1.0, 0.5),
      Point_2(2.0, 0.0),
      Point_2(3.0, 0.5)
  };

  std::vector<Point_2> Q_points = {
      Point_2(0.0, 0.0),
      Point_2(0.8, 0.4),
      Point_2(1.9, 0.1),
      Point_2(3.1, 0.6)
  };

  PolygonalCurve P(P_points);
  PolygonalCurve Q(Q_points);

  // Compute O(sqrt(n))-approximate GED
  double approxGed = GED::computeSquareRootApproxGED(P, Q);
  std::cout << "Approximate Geometric Edit Distance: " << approxGed << std::endl;

  return 0;
}

```

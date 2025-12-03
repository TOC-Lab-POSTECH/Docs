# Range Counting Oracle

The **Range Counting Oracle** provides an efficient way to access a large 2D point set without examining every point directly. Instead of accessing the entire dataset, the oracle answers simple queries that count how many points lie within an axis-aligned rectangular range.

Formally, given a point set  $$P={p_1, p_2, \dots, p_n}$$  and a query rectangle $$Q = [x_{min}, x_{max}] \times [y_{min}, y_{max}]$$, the range counting oracle returns  $$|P \cap Q|$$, the number of points contained in the given range. This interface is widely supported in modern databases, enabling algorithms to capture essential geometric information while operating in sub-linear time.

The `RangeTree` class implements this oracle by pre-processing input points into a static **2D Range Tree**. After construction, the structure can quickly return the number of points inside any query rectangle, allowing higher-level algorithms—such as clustering or geometric optimization—to work efficiently using only range-counting queries.

### Repository

{% embed url="https://github.com/TOC-Lab-POSTECH/Approximation-k-median-clustering" %}

***

### **Class Synopsis**

```cpp
class RangeTree {
 public:
  // Constructor: Builds the 2D Range Tree from a set of points
  RangeTree(const vector<Point>& points);

  // Query: Returns the number of points within [x1, x2] x [y1, y2]
  int range_count(int x1, int x2, int y1, int y2) const;

  // Getters for the bounding box of the dataset
  int getMinX() const;
  int getMaxX() const;
  int getMinY() const;
  int getMaxY() const;

 private:
  struct RangeTreeNode {
      int x_left, x_right;
      vector<int> ys;
      RangeTreeNode *left = nullptr, *right = nullptr;
  };
  
  RangeTreeNode* root;
  int minX, maxX, minY, maxY;
  
  // Internal helper functions
  RangeTreeNode* build(const vector<Point>& pts, int l, int r);
  int query(RangeTreeNode* node, int x1, int x2, int y1, int y2) const;
};

```

***

### Dependencies

* **STL**
  * `std::vector` for point and coordinate storage.
  * `<algorithm>` for `std::sort`, `std::lower_bound`, and `std::upper_bound`.
  * `<iostream>` for input/output operations.
* **Custom Headers**
  * `point.h` defining the `Point` structure with integer coordinates and weights.

***

### Design Notes

* **Static Structure**
  * The data structure is immutable after initialization. The range tree is fully constructed within the constructor, and dynamic insertions or deletions of points are not supported.
* **Integer Coordinates**
  * The system is designed for integer coordinates (`int`), aligning with the discrete space $$[\Delta]^d$$ assumption often used in sub-linear geometric algorithms. Floating-point data requires coordinate compression before usage.
* **Space-Time Trade-off**
  * Achieves $$O(\log^2 N)$$ query time by utilizing $$O(N \log N)$$ space. Each node in the primary x-tree stores a sorted vector of y-coordinates for its sub-tree to enable binary search.
* **Robustness**
  * Handles out-of-bound queries gracefully by returning 0 for non-overlapping regions.
  * Automatically calculates and stores the global bounding box (`minX`, `maxX`, `minY`, `maxY`) upon construction for quick reference.

***

### Overall Structure

{% tabs %}
{% tab title="Member Variable" %}
| Name          | Type             | Description                                        |
| ------------- | ---------------- | -------------------------------------------------- |
| `root`        | `RangeTreeNode*` | Pointer to the root node of the primary range tree |
| `minX`,`maxX` | `int`            | Minimum and maximum x-coordinate in the dataset    |
| `minY`,`maxY` | `int`            | Minimum and maximum y-coordinate in the dataset    |
{% endtab %}

{% tab title="Constructor" %}
| Name                                     | Input                     | Complexity     | Description                                                   |
| ---------------------------------------- | ------------------------- | -------------- | ------------------------------------------------------------- |
| `RangeTree(const vector<Point>& points)` | list of 2D integer points | $$O(n log n)$$ | Builds the 2D range tree and computes the global bounding box |


{% endtab %}

{% tab title="Member Function" %}


| Name                                          | Input            | Output                 | Time Complexity | Description                                             |
| --------------------------------------------- | ---------------- | ---------------------- | --------------- | ------------------------------------------------------- |
| `range_count(int x1, int x2, int y1, int y2)` | rectangle bounds | count of points(`int`) | $$O(log^2 n)$$  | Returns the number of points inside the query rectangle |
| `getMinX()`/`getMaxX()`                       | —                | `minX`/`maxX`          | $$O(1)$$        | Returns the minimum x-value/Returns the maximum x-value |
| `getMinY()`/`getMaxY()`                       | —                | coordinate(`int`)      | $$O(1)$$        | Returns the minimum y-value/Returns the maximum y-value |
{% endtab %}
{% endtabs %}

***

### Example Usage

```cpp
#include <iostream>
#include <vector>
#include "point.h"
#include "RangeCountingOracle.h" 

int main() {
    // 1. Prepare Data Points
    std::vector<Point> data;
    data.push_back(Point(1, 3));
    data.push_back(Point(2, 5));
    data.push_back(Point(4, 2));
    data.push_back(Point(5, 7));
    data.push_back(Point(7, 1));
    
    // 2. Build Range Tree (Preprocessing)
    // Time Complexity: O(N log N)
    RangeTree oracle(data);

    std::cout << "Data range X: [" << oracle.getMinX() << ", " << oracle.getMaxX() << "]\n";
    std::cout << "Data range Y: [" << oracle.getMinY() << ", " << oracle.getMaxY() << "]\n";

    // 3. Query Range Counting
    // Query: Count points in rectangle x:[1, 5], y:[1, 4]
    // Points within range: (1,3), (4,2)
    int x1 = 1, x2 = 5;
    int y1 = 1, y2 = 4;
    
    int count = oracle.range_count(x1, x2, y1, y2);
    
    std::cout << "Points in [" << x1 << "," << x2 << "]x[" 
              << y1 << "," << y2 << "]: " << count << std::endl;
    // Expected Output: 2

    return 0;
}

```

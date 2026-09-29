# Create Knots from pairs of binary trees

In 2014, Vaughan Jones developed a method to produce knots/links
from elements of Thompson's group *F*. Such elements can be represented
as pairs of binary trees, and showed that all links arise in this way.

Jones's proof that all links come from trees suggests
investigating which tree pairs give the same link.
It would be desirable to find a set of "moves" on tree-pairs
that that generate all other tree-pairs giving the same link,
analogous to Markov's theorem on braids.

This repository aims to computationally explore this space.
In particular, we can compute knot invariants from a tree pair:

```
~/treeknots$ conda activate sage
(sage) ~/treeknots$ sage
sage: from treeknots import TreePair
sage: tp = TreePair("(o((oo)(oo)))", "(((oo)(oo))o)")
sage: L = tp.get_link()
sage: L.homfly_polynomial()
L^-2*M^2 - 2*L^-2 - L^-4
sage: L.jones_polynomial()
-t^4 + t^3 + t
sage: L.get_knotinfo()
KnotInfo['K3_1']
```

The `data` folder in this repository has homfly polynomials computed
for tree-pairs with small numbers of leaves. These were generated with
commands like following:

```
(sage) ~/treeknots$ python -m treeknots.generate_data 7 -p
```
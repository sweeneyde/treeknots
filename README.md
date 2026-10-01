# Create Knots from pairs of binary trees

Much of this repository is based on
[this repository from 2020-2021](https://github.com/sweeneyde/thompson_knots/),
but rewritten and reorganized in September 2026.

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
python -m treeknots.generate_data 7 -p
```

The files `prime_knots_summary.md` and `prime_links_summary.md` list
some minimal tree-pairs representing each knot/link found.
These were generated with commands similar to the following:
```
python -m treeknots.get_knotinfo ./data/treepairs11op_sample.txt.gz
```

The `images` folder contains images of link diagrams generated like so:
```
python -m treeknots.drawing
```

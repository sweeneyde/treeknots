from collections import deque
from .tree_pair import TreePair, Vertex
from pathlib import Path

def tree_node_to_position(tree: Vertex):
    """
    >>> tree = Vertex.parse_from_word("(o(oo))")
    >>> n2p = tree_node_to_position(tree)
    >>> list(map(n2p.get, tree.in_order_vertices()))
    [(0, 0), (2, 2), (2, 0), (3, 1), (4, 0)]
    >>> tree = Vertex.parse_from_word("((oo)(oo))")
    >>> n2p = tree_node_to_position(tree)
    >>> list(map(n2p.get, tree.in_order_vertices()))
    [(0, 0), (1, 1), (2, 0), (3, 3), (4, 0), (5, 1), (6, 0)]
    """
    n2p = {}
    q = deque()
    leaves = []
    for v in tree.in_order_vertices():
        if v.is_leaf():
            leaves.append(v)
        else:
            q.append(v)
    for i, leaf in enumerate(leaves):
        n2p[leaf] = (2*i, 0)
    while q:
        v = q.popleft()
        if v.left in n2p and v.right in n2p:
            (x1, y1) = n2p[v.left]
            (x2, y2) = n2p[v.right]
            x = (x2 + x1 + y2 - y1) // 2
            y = y1 + (x - x1)
            assert y == y2 - (x - x2)
            n2p[v] = (x, y)
        else:
            q.append(v)
    return n2p

def assign_components(tp: TreePair):
    # Map (vertex, +1 or -1 for over.under) --> color index
    assignment = {}
    ([gc, signs], label_to_vertex) = tp.oriented_gauss_code_and_label_to_vertex()
    for color, comp in enumerate(gc):
        for signed_index in comp:
            sign = +1 if signed_index > 0 else -1
            v = label_to_vertex[abs(signed_index)]
            assignment[v, sign] = color
    return assignment

DIST = 100
SLOPE = 1.0
IMAGE_MARGIN = 0.75 * DIST
LEAF_SEP = 0.5 * DIST
STRAIGHT_SEGMENT = 0.6 * DIST
GAP = 0.3 * DIST
MIDLEAF_SEGMENT = LEAF_SEP + STRAIGHT_SEGMENT
STROKE_WIDTH = 0.25*DIST
OUTER_SEGMENT = 1.5 * DIST

LEAF_DOT_RADIUS = STROKE_WIDTH

def draw_treepair(
    tp: TreePair,
    *,
    filename=None,
    colors=("blue", "red", "green", "purple"),
):
    n = tp.tree1.to_word().count("o")
    if filename is None:
        filename = str(tp).replace(" / ", "-")
    top_verices_list = list(tp.tree1.in_order_vertices())
    top_vertices = set(top_verices_list)
    width = IMAGE_MARGIN+DIST*2*(n-1)+IMAGE_MARGIN
    height = OUTER_SEGMENT+DIST*2*(n-1)+OUTER_SEGMENT
    x0 = IMAGE_MARGIN
    y0 = height//2
    n2p1 = tree_node_to_position(tp.tree1)
    n2p2 = tree_node_to_position(tp.tree2)
    n2p = {}
    components = assign_components(tp)
    for v, (x, y) in n2p1.items():
        n2p[v] = (x0+DIST*x, y0 - LEAF_SEP - SLOPE*DIST*y)
    for v, (x, y) in n2p2.items():
        n2p[v] = (x0+DIST*x, y0 + LEAF_SEP + SLOPE*DIST*y)
    midleaf_position = {}
    for i, v in enumerate(top_verices_list):
        if not v.is_leaf():
            prev = top_verices_list[i - 1]
            assert prev.is_leaf()
            x, y = n2p1[top_verices_list[i - 1]]
            midleaf_position[v] = (x0+DIST*x + DIST, y0)

    import drawsvg as draw
    d = draw.Drawing(width=width, height=height)
    for v, (x, y) in n2p.items():
        if v.is_leaf():
            if v in top_vertices:
                if v.parent is None:
                    color = colors[0]
                else:
                    color = colors[components[v.parent, +1]]
                cx, cy = n2p[v.corresponding]
                d.append(draw.Line(
                    x, y, cx, cy,
                    stroke=color,
                    stroke_width=STROKE_WIDTH,
                ))
        else:
            over_color = colors[components[v, +1]]
            under_color = colors[components[v, -1]]
            lx, ly = n2p[v.left]
            rx, ry = n2p[v.right]
            plane = 1 if v in top_vertices else -1
            d.append(draw.Lines(
                lx, ly - (not v.left.is_leaf()) * GAP * plane,
                lx, ly - STRAIGHT_SEGMENT * plane,
                x - STRAIGHT_SEGMENT, y,
                x + STRAIGHT_SEGMENT, y,
                rx, ry - STRAIGHT_SEGMENT * plane,
                rx, ry - (not v.right.is_leaf()) * GAP * plane,
                stroke=over_color,
                stroke_width=STROKE_WIDTH,
                close=False,
                fill='none',
            ))
        if not v.is_leaf() and v in top_vertices:
            mx, my = midleaf_position[v]
            x1, y1 = n2p[v.corresponding]
            d.append(draw.Lines(
                x, y + GAP,
                x, y + STRAIGHT_SEGMENT,
                mx, my - MIDLEAF_SEGMENT,
                mx, my + MIDLEAF_SEGMENT,
                x1, y1 - STRAIGHT_SEGMENT,
                x1, y1 - GAP,
                stroke=under_color,
                stroke_width=STROKE_WIDTH,
                close=False,
                fill='none',
                stroke_dasharray="30,10",
            ))
    x, y = n2p[tp.tree1]
    d.append(draw.Line(
        x, y - GAP,
        x, y - OUTER_SEGMENT,
        stroke=colors[0] if tp.tree1.is_leaf() else colors[components[tp.tree1, -1]],
        stroke_width=STROKE_WIDTH,
        fill='none',
    ))
    x, y = n2p[tp.tree2]
    d.append(draw.Line(
        x, y + GAP,
        x, y + OUTER_SEGMENT,
        stroke=colors[0] if tp.tree2.is_leaf() else colors[components[tp.tree2, -1]],
        stroke_width=STROKE_WIDTH,
        fill='none',
    ))

    d.save_svg(
        Path(__file__).parent.parent
        / "images"
        / f"{filename}.svg"
    )
    print(f"{filename} ok")


def draw_treepair_trees_only(
    tp: TreePair,
    *,
    filename=None,
):
    n = tp.tree1.to_word().count("o")
    if filename is None:
        filename = str(tp).replace(" / ", "-")
    width = IMAGE_MARGIN+DIST*2*(n-1)+IMAGE_MARGIN
    height = OUTER_SEGMENT+DIST*2*(n-1)+OUTER_SEGMENT
    x0 = IMAGE_MARGIN
    y0 = height//2
    n2p1 = tree_node_to_position(tp.tree1)
    n2p2 = tree_node_to_position(tp.tree2)
    n2p = {}
    components = assign_components(tp)
    for v, (x, y) in n2p1.items():
        n2p[v] = (x0+DIST*x, y0 - LEAF_SEP - SLOPE*DIST*y)
    for v, (x, y) in n2p2.items():
        n2p[v] = (x0+DIST*x, y0 + LEAF_SEP + SLOPE*DIST*y)

    import drawsvg as draw
    d = draw.Drawing(width=width, height=height)
    for v, (x, y) in n2p.items():
        if v.is_leaf():
            d.append(draw.Circle(x, y, LEAF_DOT_RADIUS, fill="black"))
        else:
            lx, ly = n2p[v.left]
            rx, ry = n2p[v.right]
            d.append(draw.Lines(
                lx, ly, x, y, rx, ry,
                stroke="black",
                stroke_width=STROKE_WIDTH,
                close=False,
                fill='none',
            ))

    d.save_svg(
        Path(__file__).parent.parent
        / "images"
        / f"{filename}_trees.svg"
    )
    print(f"{filename} ok")


def main():
    def draw(filename, word1, word2, **kwargs):
        tp = TreePair(word1, word2)
        draw_treepair_trees_only(tp, filename=filename)
        draw_treepair(tp, filename=filename, **kwargs)
    draw("1leaf",
         "o",
         "o")
    draw("2leaf",
         "(oo)",
         "(oo)")
    draw("3leaf_1comp",
         "(o(oo))",
         "((oo)o)")
    draw("3leaf_3comp",
         "(o(oo))",
         "(o(oo))")
    draw("hopf",
         "((oo)((oo)o))",
         "(o((oo)(oo)))")
    draw("trefoil",
         "(((oo)(oo))o)",
         "(o((oo)(oo)))")
    draw("borromean",
         "(o(o(o(((o(oo))(oo))o))))",
         "((((oo)((oo)o))(oo))(oo))")

if __name__ == "__main__":
    import sys
    if len(sys.argv) == 2:
        wordpair = sys.argv[1]
        tp = TreePair(*wordpair.split(" / "))
        draw_treepair(tp)
    elif len(sys.argv) == 1:
        main()
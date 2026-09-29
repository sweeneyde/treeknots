from pathlib import Path
import gzip
from .tree_pair import Vertex, TreePair

KNOWN_TRIVIAL_WORD_PAIR_SET = None

def is_known_trivial_word_pair(word1, word2):
    global KNOWN_TRIVIAL_WORD_PAIR_SET
    if KNOWN_TRIVIAL_WORD_PAIR_SET is None:
        path = Path(__file__).parent.parent / "data" / "known_unknots.txt.gz"
        with gzip.open(path, "rt", encoding="ascii") as f:
            KNOWN_TRIVIAL_WORD_PAIR_SET = set(f.read().splitlines())
    return f"{word1} / {word2}" in KNOWN_TRIVIAL_WORD_PAIR_SET

def mirror_across_trees(tp):
    """
    >>> tp = TreePair("((o(oo))(oo))", "(((oo)o)(oo))")
    >>> mirror_across_trees(tp).to_words()
    ('(((oo)o)(oo))', '((o(oo))(oo))')
    """
    tp = tp.copy()
    tp.tree1, tp.tree2 = tp.tree2, tp.tree1
    # "corresponding" links stay the same.
    return tp

def mirror_left_right(tp):
    """
    >>> tp = TreePair("((o(oo))(oo))", "(((oo)o)(oo))")
    >>> mirror_left_right(tp).to_words()
    ('((oo)((oo)o))', '((oo)(o(oo)))')
    """
    tp = tp.copy()
    for t in tp.tree1, tp.tree2:
        for v in list(t.in_order_vertices()):
            v.left, v.right = v.right, v.left
    # "corresponding" links stay the same.
    return tp

def subsquare_decompositions(tp):
    """
    >>> tp = TreePair("(o((oo)o))", "(o(o(oo)))")
    >>> list(subsquare_decompositions(tp))
    [(('(oo)', '(oo)'), ('((oo)o)', '(o(oo))'))]
    >>> tp = TreePair("((o(oo))o)", "((o(oo))o)")
    >>> list(subsquare_decompositions(tp))
    [(('(oo)', '(oo)'), ('(o(oo))', '(o(oo))')), (('((oo)o)', '((oo)o)'), ('(oo)', '(oo)'))]
    """
    tree1, tree2 = tp.tree1, tp.tree2
    left_right_pairs_1 = {
        (v.leftmost_descendent(), v.rightmost_descendent()) : v
        for v in tree1.in_order_vertices()
        if v is not tree1
        if not v.is_leaf()
    }
    left_right_pairs_2 = {
        (v.leftmost_descendent().corresponding,
            v.rightmost_descendent().corresponding) : v
        for v in tree2.in_order_vertices()
        if v is not tree2
        if not v.is_leaf()
    }
    for pair, v1 in left_right_pairs_1.items():
        if (v2 := left_right_pairs_2.get(pair)) is None:
            continue
        subsquare_word_pair = (v1.to_word(), v2.to_word())
        sans_subsquare_word_pair = (tree1.to_word(omitted=v1),
                                    tree2.to_word(omitted=v2))
        yield sans_subsquare_word_pair, subsquare_word_pair

def subsquare_moves(tp):
    # labelled #1 in my google doc
    for sans_subsquare_word_pair, subsquare_word_pair in subsquare_decompositions(tp):
        if is_known_trivial_word_pair(*sans_subsquare_word_pair):
            yield TreePair(*subsquare_word_pair)
        if is_known_trivial_word_pair(*subsquare_word_pair):
            yield TreePair(*sans_subsquare_word_pair)

def local_move_1(tp):
    """
    >>> tp = TreePair("(((oo)o)((o(oo))o))", "(o((o(o(oo)))(oo)))")
    >>> [tp1.to_words() for tp1 in local_move_1(tp)]
    [('(((oo)o)((oo)o))', '(o((o(oo))(oo)))')]
    >>> tp = TreePair("(o((o(oo))o))", "((o((oo)o))o)")
    >>> [tp1.to_words() for tp1 in local_move_1(tp)]
    []
    """
    # labelled #2 in my google doc
    for v in tp.tree1.in_order_vertices():
        match v:
            case Vertex(
                left=Vertex(
                    left=Vertex(left=None, right=None) as leaf1,
                    right=Vertex(
                        left=Vertex(left=None, right=None) as leaf2,
                        right=Vertex(left=None, right=None) as leaf3,
                    )
                ),
                right=Vertex(left=None, right=None) as leaf4,
            ) if (
                leaf1.corresponding.parent is leaf2.corresponding.parent
                and leaf3.corresponding.parent is leaf4.corresponding.parent
            ):
                copy, mapping = tp.copy_with_mapping()
                v_copy = mapping[v]
                LR = v_copy.left.right
                Lc = v_copy.left.corresponding
                LR.left = LR.right = None
                Lc.left = Lc.right = None
                copy.relink_corresponding_and_parents()
                yield copy

def local_move_2(tp):
    """
    >>> tp = TreePair("((o((oo)(oo)))o)", "((((oo)o)o)(oo))")
    >>> [tp1.to_words() for tp1 in local_move_2(tp)]
    [('((o(oo))o)', '((oo)(oo))')]
    >>> tp = TreePair("(o((((oo)o)(oo))o))", "(((oo)o)((oo)(oo)))")
    >>> [tp1.to_words() for tp1 in local_move_2(tp)]
    [('(((oo)(oo))o)', '(o((oo)(oo)))')]
    >>> tp = TreePair("(o(oo))", "((oo)o)")
    >>> [tp1.to_words() for tp1 in local_move_2(tp)]
    [('o', 'o')]
    >>> tp = TreePair("((oo)(oo))", "(o((oo)o))")
    >>> [tp1.to_words() for tp1 in local_move_2(tp)]
    []
    """
    # labelled #5 in my google doc
    for v2 in tp.tree2.in_order_vertices():
        match v2:
            case Vertex(
                left=Vertex(
                    left=Vertex(left=None, right=None) as leaf1,
                    right=Vertex(left=None, right=None) as leaf2,
                ),
                right=Vertex(left=None, right=None) as leaf3,
            ) if (
                leaf1.corresponding is leaf1.corresponding.parent.left
                and leaf2.corresponding.parent is leaf3.corresponding.parent
            ):
                copy, mapping = tp.copy_with_mapping()
                v2_copy = mapping[v2]
                v1_copy = mapping[leaf1.corresponding.parent]
                v2c = v2_copy.corresponding
                v1R = v1_copy.right
                v2_copy.left = v2_copy.right = None
                v2c.left = v2c.right = None
                v1_copy.left, v1_copy.right = v1R.left, v1R.right
                copy.relink_corresponding_and_parents()
                yield copy

def local_move_3(tp):
    """
    >>> tp = TreePair("(o((o(o(oo)))o))", "((o((oo)o))(oo))")
    >>> [tp1.to_words() for tp1 in local_move_3(tp)]
    [('(o(((oo)o)o))', '((oo)((oo)o))')]
    """
    # labelled #7 in my google doc
    for v1 in tp.tree1.in_order_vertices():
        match v1:
            case Vertex(
                left=Vertex(left=None, right=None) as leaf1,
                right=Vertex(
                    left=Vertex(left=None, right=None) as leaf2,
                    right=Vertex(
                        left=Vertex(left=None, right=None) as leaf3,
                        right=Vertex(left=None, right=None) as leaf4,
                    ),
                )
            ) if (
                leaf1.corresponding.parent is leaf2.corresponding.parent
                and leaf1.corresponding.parent.parent is leaf3.corresponding.parent
            ):
                copy, mapping = tp.copy_with_mapping()
                v1_copy = mapping[v1]
                v2a = mapping[leaf3].corresponding.parent
                v2b = mapping[leaf4].corresponding
                v1_copy.left = Vertex.parent_of(Vertex.leaf(), Vertex.leaf())
                v1_copy.right = Vertex.leaf()
                v2a.left = v2a.right = None
                v2b.left, v2b.right = Vertex.leaf(), Vertex.leaf()
                copy.relink_corresponding_and_parents()
                yield copy

def local_move_4(tp):
    """
    >>> tp = TreePair('(o(oo))', '(o(oo))')
    >>> [tp1.to_words() for tp1 in local_move_4(tp)]
    [('((oo)o)', '((oo)o)')]
    >>> tp = TreePair("((o(o(oo)))o)", "((oo)(o(oo)))")
    >>> [tp1.to_words() for tp1 in local_move_4(tp)]
    [('((o((oo)o))o)', '(((oo)o)(oo))')]
    """
    # labelled #10 in my google doc
    for leaf1 in tp.tree1.in_order_vertices():
        leaf2 = leaf1.corresponding
        if (
            leaf1.is_leaf()
            and (p1 := leaf1.parent) is not None
            and leaf1 is p1.left
            and (pp1 := p1.parent) is not None
            and p1 is pp1.right
            and (p2 := leaf2.parent) is not None
            and leaf2 is p2.left
            and (pp2 := p2.parent) is not None
            and p2 is pp2.right
        ):
            copy, mapping = tp.copy_with_mapping()
            p3 = mapping[p1]
            p4 = mapping[p2]
            pp3 = mapping[pp1]
            pp4 = mapping[pp2]
            pp3.left, pp3.right, p3.left, p3.right = p3, p3.right, pp3.left, p3.left
            pp4.left, pp4.right, p4.left, p4.right = p4, p4.right, pp4.left, p4.left
            copy.relink_corresponding_and_parents()
            yield copy

def local_move_5(tp):
    """
    >>> tp = TreePair('(o((o(oo))(oo)))', '(o((oo)((oo)o)))')
    >>> [tp1.to_words() for tp1 in local_move_5(tp)]
    [('(o(oo))', '(o(oo))')]
    """
    # labelled #9 in my google doc
    for leaf_a1 in tp.tree1.in_order_vertices():
        if not leaf_a1.is_leaf(): continue
        leaf_a2 = leaf_a1.corresponding
        if (q2 := leaf_a2.parent) is None: continue
        if leaf_a2 is not q2.left: continue
        if not (leaf_b2 := q2.right).is_leaf(): continue
        leaf_b1 = leaf_b2.corresponding
        if (p1 := leaf_b1.parent) is None: continue
        if leaf_b1 is not p1.left: continue
        if not (leaf_c1 := p1.right).is_leaf(): continue
        leaf_c2 = leaf_c1.corresponding
        if (p2 := leaf_c2.parent) is None: continue
        if leaf_c2 is not p2.left: continue
        if not (leaf_d2 := p2.right).is_leaf(): continue
        leaf_d1 = leaf_d2.corresponding
        if (q1 := leaf_d1.parent) is None: continue
        if leaf_d1 is not q1.left: continue
        if not (leaf_e1 := q1.right).is_leaf(): continue
        leaf_e2 = leaf_e1.corresponding
        if (pp1 := leaf_a1.parent) is None: continue
        if pp1.right is not p1: continue
        if (pp2 := leaf_e2.parent) is None: continue
        if pp2.left is not p2: continue
        copy, mapping = tp.copy_with_mapping()
        pp1 = mapping[pp1]
        q1 = mapping[q1]
        pp2 = mapping[pp2]
        q2 = mapping[q2]
        pp1.left = pp1.right = None
        q1.left = q1.right = None
        pp2.left = pp2.right = None
        q2.left = q2.right = None
        copy.relink_corresponding_and_parents()
        yield copy

def global_move_1(tp):
    """
    >>> tp = TreePair('(o(((oo)o)((oo)o)))', '((oo)((o(oo))(oo)))')
    >>> [tp1.to_words() for tp1 in global_move_1(tp)]
    [('((((oo)o)o)((oo)o))', '(o(o((o(oo))(oo))))')]
    """
    # labelled #3 in my google doc
    tree1 = tp.tree1
    if tree1.is_leaf() or not tree1.left.is_leaf():
        return
    tree2 = tp.tree2
    if tree2.left.is_leaf():
        return
    if not tree2.left.left.is_leaf():
        return
    if not tree2.left.right.is_leaf():
        return
    copy = tp.copy()
    v1 = copy.tree1
    v2 = copy.tree2
    L = v2.left
    LL = L.left
    LR = L.right
    LRc = LR.corresponding
    R = v2.right
    copy.tree1 = v1.right
    leaf1 = Vertex.leaf()
    leaf2 = Vertex.leaf()
    LRc.left, LRc.right = leaf1, leaf2
    v2.left, v2.right = LL, L
    L.left, L.right = LR, R
    copy.relink_corresponding_and_parents()
    yield copy

def global_move_2(tp):
    """
    >>> tp = TreePair("(o((o((oo)o))(oo)))", "(((o(oo))(o(oo)))o)")
    >>> [tp1.to_words() for tp1 in global_move_2(tp)]
    [('((o((oo)o))(oo))', '(((oo)(o(oo)))o)')]
    """
    # labelled #6 in my google doc
    match [tp.tree1, tp.tree2]:
        case [
            Vertex(
                left=Vertex(left=None, right=None) as leaf_a1,
                right=Vertex(
                    left=Vertex(
                        left=Vertex(left=None, right=None) as leaf_b1
                    ),
                    right=Vertex(
                        right=Vertex(left=None, right=None) as leaf_c1
                    )
                )
            ),
            Vertex(
                left=Vertex(
                    left=Vertex(
                        left=Vertex(left=None, right=None) as leaf_a2,
                        right=Vertex(
                            left=Vertex(left=None, right=None) as leaf_b2
                        )
                    ),
                ),
                right=Vertex(left=None, right=None) as leaf_c2,
            )
        ]:
            assert leaf_a2 is leaf_a1.corresponding
            assert leaf_b2 is leaf_b1.corresponding
            assert leaf_c2 is leaf_c1.corresponding
            copy = tp.copy()
            copy.tree1 = copy.tree1.right
            left = copy.tree2.left
            left.left = left.left.right
            copy.relink_corresponding_and_parents()
            yield copy

MOVES = [
    global_move_1,
    global_move_2,
    local_move_1,
    local_move_2,
    local_move_3,
    local_move_4,
    local_move_5,
    local_move_5,
]

def all_moves_with_reflections(tp):
    yield mirror_across_trees(tp)
    yield mirror_left_right(tp)
    yield from subsquare_moves(tp)
    for moves in MOVES:
        yield from moves(tp)

def all_moves_no_reflections(tp):
    yield mirror_left_right(mirror_across_trees(tp))
    yield from subsquare_moves(tp)
    # So we only have to write each move once,
    # we will apply reflections as long as we promise
    # to reflect back afterwards.
    for reflection in [lambda x: x,
                        mirror_across_trees,
                        mirror_left_right]:
        tp0 = reflection(tp)
        for moves in MOVES:
            for tp1 in moves(tp0):
                yield reflection(tp1)

if __name__ == "__main__":
    import doctest
    doctest.testmod()
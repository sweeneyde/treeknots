from collections import defaultdict

class Vertex:
    __slots__ = "parent", "left", "right", "corresponding"

    def is_leaf(self):
        return self.left is None

    @classmethod
    def parent_of(cls, left : "Vertex", right: "Vertex") -> "Vertex":
        parent = cls()
        left.parent = parent
        right.parent = parent
        parent.parent = None
        parent.left = left
        parent.right = right
        return parent

    @classmethod
    def leaf(cls) -> "Vertex":
        leaf = cls()
        leaf.left = leaf.right = leaf.parent = None
        return leaf

    @classmethod
    def parse_from_word(cls, word : str) -> "Vertex":
        """
        >>> Vertex.parse_from_word('(o(oo))').to_word()
        '(o(oo))'
        """
        OPEN = object()
        stack = []
        for c in word:
            if c == '(':
                stack.append(OPEN)
            elif c == 'o':
                stack.append(cls.leaf())
            elif c == ')':
                if len(stack) < 3:
                    raise ValueError("unmatched ')'")
                matching, left, right = stack[-3:]
                if left is OPEN or right is OPEN:
                    raise ValueError("too few children")
                if matching is not OPEN:
                    raise ValueError("too many children")
                stack[-3:] = [cls.parent_of(left, right)]
        if len(stack) != 1:
            raise ValueError("malformed word")
        [result] = stack
        assert result.to_word() == word
        return result

    def to_word(self, *, omitted=object()):
        NEEDS_CLOSED = object()
        result = []
        stack = [self]
        while stack:
            v = stack.pop()
            if v is NEEDS_CLOSED:
                result.append(")")
            elif v is omitted or v.is_leaf():
                result.append("o")
            else:
                result.append("(")
                stack.append(NEEDS_CLOSED)
                stack.append(v.right)
                stack.append(v.left)
        return "".join(result)

    __str__ = to_word

    def in_order_vertices(self):
        if self.is_leaf():
            yield self
        else:
            yield from self.left.in_order_vertices()
            yield self
            yield from self.right.in_order_vertices()

    def leftmost_descendent(self):
        v = self
        while (left := v.left) is not None:
            v = left
        return v

    def rightmost_descendent(self):
        v = self
        while (right := v.right) is not None:
            v = right
        return v

    def copy(self):
        if self.is_leaf():
            return Vertex.leaf()
        else:
            return Vertex.parent_of(self.left.copy(), self.right.copy())


class TreePair:
    __slots__ = "tree1", "tree2"

    def __init__(self, word1, word2):
        c1 = word1.count("o")
        c2 = word2.count("o")
        if c1 != c2:
            raise ValueError("Different sized trees.")
        tree1 = Vertex.parse_from_word(word1)
        tree2 = Vertex.parse_from_word(word2)
        self.tree1 = tree1
        self.tree2 = tree2
        self.link_corresponding()

    def link_corresponding(self):
        vertices1 = list(self.tree1.in_order_vertices())
        vertices2 = list(self.tree2.in_order_vertices())
        for v1, v2 in zip(vertices1, vertices2, strict=True):
            v1.corresponding = v2
            v2.corresponding = v1

    def relink_corresponding_and_parents(self):
        self.tree1.parent = None
        self.tree2.parent = None
        vertices1 = list(self.tree1.in_order_vertices())
        vertices2 = list(self.tree2.in_order_vertices())
        for v1, v2 in zip(vertices1, vertices2, strict=True):
            v1.corresponding = v2
            v2.corresponding = v1
        for v in vertices1 + vertices2:
            if not v.is_leaf():
                v.left.parent = v
                v.right.parent = v

    def copy(self):
        return __class__(self.tree1.to_word(), self.tree2.to_word())

    def copy_with_mapping(self):
        other = self.copy()
        self_vertices = list(self.tree1.in_order_vertices()) + list(self.tree2.in_order_vertices())
        other_vertices = list(other.tree1.in_order_vertices()) + list(other.tree2.in_order_vertices())
        d = dict(zip(self_vertices, other_vertices))
        return (other, d)

    def __str__(self):
        return f"{self.tree1} / {self.tree2}"

    def to_words(self):
        return (self.tree1.to_word(), self.tree2.to_word())

    def oriented_gauss_code(self):
        """
        >>> tp = TreePair("(o(o((oo)o)))", "(((oo)o)(oo))")
        >>> ((comp1, comp2), signs) = tp.oriented_gauss_code()
        >>> comp1
        [1, 5, 2, -4, -8, 7, -6, -2]
        >>> comp2
        [-1, -7, -3, 4, 8, 3, 6, -5]
        >>> signs
        [-1, 1, -1, -1, 1, -1, -1, 1]
        """
        tree1, tree2 = self.tree1, self.tree2
        LEFT, RIGHT, CORRESPONDING, PARENT = object(), object(), object(), object()
        vertices1 = list(tree1.in_order_vertices())
        vertices2 = list(tree2.in_order_vertices())

        def get_next_node_and_direction(v, d):
            # Given a non-leaf node v and a direction d,
            # follow the string along to the next vertex and direciton.
            if d is PARENT:
                parent = v.parent
                if parent is None:
                    return (tree2 if v is tree1 else tree1, CORRESPONDING)
                else:
                    return (parent, LEFT if v is parent.right else RIGHT)
            elif d is CORRESPONDING:
                return (v.corresponding, PARENT)
            elif d is LEFT:
                left = v.left
                if left.is_leaf():
                    other_leaf = left.corresponding
                    other_parent = other_leaf.parent
                    return (other_parent, LEFT if other_leaf is other_parent.right else RIGHT)
                else:
                    return (left, CORRESPONDING)
            elif d is RIGHT:
                right = v.right
                if right.is_leaf():
                    other_leaf = right.corresponding
                    other_parent = other_leaf.parent
                    return (other_parent, LEFT if other_leaf is other_parent.right else RIGHT)
                else:
                    return (right, CORRESPONDING)
            else:
                raise AssertionError()

        vertex_to_label = {}
        for v in vertices1 + vertices2:
            if not v.is_leaf():
                label = len(vertex_to_label) + 1
                vertex_to_label[v] = label

        OVER, UNDER = +1, -1
        DIRECTON_TO_OVER_UNDER = {
            LEFT: OVER, RIGHT: OVER,
            PARENT: UNDER, CORRESPONDING: UNDER
        }

        # Labels Over/undercrossings visited, negative if undercrossings
        seen = set()
        directions_taken = defaultdict(set)
        gauss_components = []

        def traverse(v0, d0):
            v, d = v0, d0
            while True:
                signed_label = vertex_to_label[v] * DIRECTON_TO_OVER_UNDER[d]
                if signed_label in seen:
                    return
                yield signed_label
                seen.add(signed_label)
                directions_taken[v].add(d)
                v, d = get_next_node_and_direction(v, d)

        for v0 in vertex_to_label:
            for d0 in (LEFT, PARENT):
                component = list(traverse(v0, d0))
                if component:
                    gauss_components.append(component)
        assert len(seen) == 2*len(vertex_to_label)

        bottom_halfplane = set(vertices2)
        crossing_signs = []
        for v in vertex_to_label:
            taken = directions_taken[v]
            sign = +1
            if CORRESPONDING in taken:
                sign = -sign
            if LEFT in taken:
                sign = -sign
            if v in bottom_halfplane:
                sign = -sign
            crossing_signs.append(sign)
        return [gauss_components, crossing_signs]

    def get_link(self):
        from sage.all import Link
        return Link(self.oriented_gauss_code())


if __name__ == "__main__":
    import doctest
    doctest.testmod()

from functools import cache
from itertools import product, chain

@cache
def tree_words(num_leaves):
    if num_leaves == 1:
        return ('o',)
    all_words = []
    for k in range(1, num_leaves):
        lefts = tree_words(k)
        rights = tree_words(num_leaves - k)
        for left, right in product(lefts, rights):
            all_words.append(f"({left}{right})")
    return tuple(all_words)

def all_word_pairs(num_leaves):
    return map(" / ".join, product(tree_words(num_leaves), repeat=2))

def all_word_pairs_up_to(num_leaves):
    return chain.from_iterable(map(all_word_pairs, range(1, num_leaves + 1)))

def prime_word_pairs(num_leaves):
    """
    >>> prime_word_pairs(num_leaves)
    """
    tws = tree_words(num_leaves)
    subtrees = list(map(subtree_ranges, tws))
    for word1, subtrees1 in zip(tws, subtrees):
        for word2, subtrees2 in zip(tws, subtrees):
            if subtrees1.isdisjoint(subtrees2):
                yield f"{word1} / {word2}"

def prime_word_pairs_up_to(num_leaves):
    for n in range(3, num_leaves + 1):
        yield from prime_word_pairs(n)

def subtree_ranges(word):
    result = []
    stack = []
    count = 0
    for c in word:
        if c == "(":
            stack.append(count)
        elif c == ")":
            result.append((stack.pop(), count))
        elif c == "o":
            count += 1
    assert not stack
    assert result[-1] == (0, count)
    result.pop()
    return frozenset(result)


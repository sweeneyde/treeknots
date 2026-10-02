"""
Try to evaluate how good our set of moves is.
"""

from collections import deque, defaultdict, Counter
import functools
import gzip
from math import comb
from pathlib import Path

from tqdm import tqdm

from .tree_pair import TreePair

def group_by_reachability(word_pairs, moves):
    # union/find data structure
    parent = {x: x for x in word_pairs}
    size = {x: 1 for x in word_pairs}
    def find(x):
        while (ppx := parent[parent[x]]) != x:
            parent[x] = ppx
            x = ppx
        return x
    def union(x, y):
        x = find(x)
        y = find(y)
        if x == y:
            return
        if size[x] < size[y]:
            x, y = y, x
        parent[y] = x
        size[x] += size[y]
    num_added_to_q = 0
    q = deque()
    done = set()
    def link_neighbors(x : str):
        if x in done:
            return
        done.add(x)
        tp = TreePair(*x.split(" / "))
        for move in moves:
            for tp1 in move(tp):
                y = str(tp1)
                if y not in parent:
                    # print(f"{move.__name__}: {x} --> {y} not in parent")
                    size[y] = 1
                    parent[y] = y
                    q.append(y)
                    nonlocal num_added_to_q
                    num_added_to_q += 1
                union(x, y)
    for x in tqdm(word_pairs, "uniting moves"):
        link_neighbors(x)
    print(f"Queue has {len(q)} more to check")
    while q:
        link_neighbors(q.popleft())
    bins = defaultdict(list)
    for x in word_pairs:
        bins[find(x)].append(x)
    result = list(bins.values())
    for bin in result:
        bin.sort(key=lambda s: (len(s), s))
    result.sort(key=lambda arr:(len(arr), arr))
    return result, num_added_to_q

def make_flipped_move(flip, move):
    @functools.wraps(move)
    def flipped_move(tp):
        return map(flip, move(flip(tp)))
    flipped_move.__name__ = f"{flip.__name__}__{move.__name__}"
    return flipped_move

def do_evaluation(
    infile_path,
    knot_type="0_1",
    restrict_to_max_leaves=None,
    allow_reflections=False,
    allow_pastings=False,
):
    from .moves import (
        MOVES,
        mirror_left_right,
        mirror_across_trees,
        subsquare_moves,
    )
    def move_mirror_left_right(tp):
        yield mirror_left_right(tp)
    def move_mirror_across_trees(tp):
        yield mirror_across_trees(tp)
    def move_rotate(tp):
        yield mirror_left_right(mirror_across_trees(tp))
    moves = list(MOVES)
    if allow_reflections:
        moves.append(move_mirror_left_right)
        moves.append(move_mirror_across_trees)
    else:
        moves.append(move_rotate)
        for move in MOVES:
            moves.append(make_flipped_move(mirror_left_right, move))
            moves.append(make_flipped_move(mirror_across_trees, move))
    if allow_pastings:
        moves.append(subsquare_moves)
    homfly = {
        "0_1": "1",
        "4_1": "-L^2 + M^2 - 1 - L^-2",
        "3_1": "-L^4 + L^2*M^2 - 2*L^2",
        "L2a1": "L^3*M^-1 - L*M + L*M^-1",
    }[knot_type]
    homfly = f"homfly: {homfly}\n"
    word_pairs = []
    with gzip.open(infile_path, "rt") as f:
        lines = iter(f)
        for line in lines:
            if line == homfly:
                break
        for line in lines:
            if line.startswith("homfly"):
                break
            word_pairs.append(line.strip())
    bins, num_extra = group_by_reachability(word_pairs, moves)
    if restrict_to_max_leaves is not None:
        for bin in bins:
            bin[:] = [pair for pair in bin
                      if pair.count("o")//2 <= restrict_to_max_leaves]
        bins = list(filter(None, bins))
        bins.sort(key=lambda x: (len(x), x))

    namecode = Path(infile_path).name.removeprefix("treepairs").removesuffix(".txt.gz")
    result_name = f"({knot_type}){'' if restrict_to_max_leaves is None else f'_max{restrict_to_max_leaves}'}_{namecode}"

    with open(Path(__file__).parent.parent / "evaluations" / f"comp{result_name}.txt", "w") as f:
        for bin in bins:
            for word_pair in bin:
                print(word_pair, file=f)
            print("="*80, file=f)

    lengths = list(map(len, bins))
    total = sum(lengths)
    connections = sum(comb(L, 2) for L in lengths)
    with open(Path(__file__).parent.parent / "evaluations" / f"eval{result_name}.md", "w") as f:
        print("# TreePair Move Evaluations", file=f)
        print("", file=f)
        print(f"* knot type: {knot_type} (`{homfly}`)", file=f)
        print("", file=f)
        print(f"* using data from file: {Path(infile_path).name}", file=f)
        print("", file=f)
        print(f"* restrict to max leaves: {restrict_to_max_leaves}", file=f)
        print("", file=f)
        print(f"* Also reached {num_extra} unlisted treepairs.", file=f)
        print("", file=f)
        print("* Moves used:", file=f)
        for move in moves:
            print(f"    * `{move.__name__}`", file=f)
        print(f"* Largest component has size {max(lengths):,}/{total:,}={max(lengths)/total:.5%}", file=f)
        print(f"* Connection probability: {connections:,}/{comb(total, 2):,}={connections/comb(total,2):.5%}", file=f)
        print(f"* Component sizes:", file=f)
        print(f"```", file=f)
        for length, count in Counter(lengths).items():
            print(f"{length}," * count, file=f)
        print(f"```", file=f)
        print(file=f)
    return (max(lengths)/total, connections / comb(total, 2))

def multiple_evaluations(
        start,
        max_intermediate_leaves,
        knot_type,
        flags,
        allow_reflections,
        allow_pastings,

):
    results = []
    for n in range(start, max_intermediate_leaves+1):
        max_score, prob_score = do_evaluation(
            infile_path=Path(__file__).parent.parent / "data" / f"treepairs{n}{flags}.txt.gz",
            knot_type=knot_type,
            restrict_to_max_leaves=start,
            allow_reflections=allow_reflections,
            allow_pastings=allow_pastings,
        )
        results.append((n, max_score, prob_score))
        if max_score == 1.0:
            break
    print(f"{knot_type=}, {flags=}, {allow_pastings=}, {allow_reflections=}")
    for n, max_score, prob_score in results:
        print(f"{start} leaf treepairs using up to {n} leaves: {prob_score:.5%} reachability (largest={max_score:.5%})")

if __name__ == "__main__":
    multiple_evaluations(
        start=5,
        max_intermediate_leaves=8,
        knot_type="0_1",
        flags="",
        allow_reflections=False,
        allow_pastings=False,
    )
    multiple_evaluations(
        start=6,
        max_intermediate_leaves=8,
        knot_type="0_1",
        flags="",
        allow_reflections=False,
        allow_pastings=False,
    )
    multiple_evaluations(
        start=5,
        max_intermediate_leaves=8,
        knot_type="3_1",
        flags="",
        allow_reflections=False,
        allow_pastings=False,
    )
    multiple_evaluations(
        start=5,
        max_intermediate_leaves=8,
        knot_type="L2a1",
        flags="",
        allow_reflections=False,
        allow_pastings=False,
    )

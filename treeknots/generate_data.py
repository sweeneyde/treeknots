from collections import defaultdict, deque
from pathlib import Path
import gzip

from tqdm import tqdm

from .tree_pair import TreePair
from .generate_word_pairs import (
    all_word_pairs_up_to,
    prime_word_pairs_up_to,
)
from .moves import (
    all_moves_no_reflections,
    all_moves_with_reflections,
    quick_test_maximality,
)

def partition_basic(word_pairs, **_):
    result = defaultdict(list)
    for x in tqdm(list(word_pairs)):
        L = TreePair(*x.split(" / ")).get_link()
        result[str(L.homfly_polynomial())].append(x)
    return result

def partition_using_moves(word_pairs, moves):
    # Union-find data structure
    word_pairs = list(word_pairs)
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

    q = deque()
    done = set()

    def link_neighbors(x : str):
        if x in done:
            return
        done.add(x)
        tp = TreePair(*x.split(" / "))
        for tp1 in moves(tp):
            y = str(tp1)
            if y not in parent:
                size[y] = 1
                parent[y] = y
                q.append(y)
            union(x, y)

    for x in tqdm(word_pairs, "uniting moves"):
        link_neighbors(x)
    print(f"Queue has {len(q)} more to check")
    while q:
        link_neighbors(q.popleft())
    bins = defaultdict(list)
    for x in word_pairs:
        bins[find(x)].append(x)
    bins = list(bins.values())
    del parent, size
    print(f"reduced from {len(word_pairs)=} to {len(bins)=}")

    result = defaultdict(list)
    for bin in tqdm(bins, "calc on each bin"):
        rep = min(bin, key=lambda x: (len(x), x))
        L = TreePair(*rep.split(" / ")).get_link()
        p = str(L.homfly_polynomial())
        result[p].extend(bin)

    return result

def homfly_worker(x : str):
    tp = TreePair(*x.split(" / "))
    L = tp.get_link()
    h = str(L.homfly_polynomial())
    return x, h

def construct_one_example_each(word_pairs, **_):
    result = {}
    def iterator():
        yield "o / o"
        for x in tqdm(word_pairs):
            tp = TreePair(*x.split(" / "))
            if quick_test_maximality(tp):
                yield x
    import multiprocessing as mp
    mp.set_start_method("spawn")
    with mp.Pool(8) as pool:
        for x, h in pool.imap_unordered(homfly_worker, iterator()):
            if h not in result or len(x) < len(result[h]) or (len(x) == len(h) and x > result[h]):
                result[h] = x
    return {h: [x] for h, x in result.items()}

def main():
    import argparse
    parser = argparse.ArgumentParser(prog="generate data",
                                     description="generate trees pairs and find their")
    parser.add_argument("n", help="the number of leaves", type=int)
    parser.add_argument("-p",
                        help="set this flag to only consider tree pairs with no inserted sub pairs",
                        action="store_true")
    parser.add_argument("-m",
                        help="set this flag to use moves before computing link invariants",
                        action="store_true")
    parser.add_argument("-o",
                        help="set this flag to only consider 1-component links (knots)",
                        action="store_true")
    parser.add_argument("-r",
                        help="set this flag to allow moves that reflect",
                        action="store_true")
    parser.add_argument("--sample",
                        help="set this flag to produce only one per",
                        action="store_true")
    args = parser.parse_args()
    if args.r:
        assert args.m
    if args.sample and args.m:
        raise ValueError()
    if args.sample:
        args.p = True

    word_pairs_func = prime_word_pairs_up_to if args.p else all_word_pairs_up_to
    word_pairs = word_pairs_func(args.n)
    if args.o:
        word_pairs = (wordpair
                      for wordpair in word_pairs
                      if len(TreePair(*wordpair.split(" / ")).oriented_gauss_code()[0]) == 1)
    partition_func = (partition_using_moves if args.m
                      else construct_one_example_each if args.sample
                      else partition_basic)
    moves = all_moves_with_reflections if args.r else all_moves_no_reflections

    result_dd = partition_func(word_pairs, moves=moves)

    result = list(result_dd.items())
    for k, arr in result:
        arr.sort(key=lambda s: (len(s), s))
    result.sort()
    result.sort(key=lambda kv: (len(kv[1]), len(kv[0]), kv[0], kv[1]))

    tags = ""
    for tag in "ompr":
        if getattr(args, tag):
            tags += tag
    if args.sample:
        tags += "_sample"
    filename = f"treepairs{args.n}{tags}.txt.gz"
    path = Path(__file__).parent.parent / "data" / filename
    with gzip.open(path, 'wt', encoding='ascii') as f:
        for k, arr in result:
            print(f"homfly: {k}", file=f)
            for wordpair in arr:
                print(wordpair, file=f)

    print(f"wrote to {path}")

    if args.n == 7 and args.p and not args.m and not args.sample:
        new_path = path.parent / "known_unknots.txt.gz"
        [arr] = [arr for (k, arr) in result if k == "1"]
        with gzip.open(new_path, "wt", encoding="ascii") as f:
            for wordpair in arr:
                print(wordpair, file=f)
        print(f"wrote to {new_path}")

if __name__ == "__main__":
    main()
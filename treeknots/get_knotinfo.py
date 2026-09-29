from collections import Counter
import gzip
from itertools import batched
from pathlib import Path
import sys

from tqdm import tqdm

from .tree_pair import TreePair


def simplify_one_component(L):
    import snappy
    sL = snappy.Link(L.pd_code())
    sL.simplify("basic")
    L1 = type(L)([[x+1 for x in tup] for tup in sL.PD_code()])
    return L1

def prime_knots_summary(data_filepath):
    with gzip.open(data_filepath, "rt", encoding="ascii") as f:
        lines = f.read().splitlines()
    data = {k.removeprefix("homfly: "): x for k, x in batched(lines, 2)}
    data["1"] = "o / o"
    max_leaves = max(rep.count("o") // 2 for rep in data.values())
    results_by_num_leaves = [[] for _ in range(1 + max_leaves)]
    pbar = tqdm(total=len(data), dynamic_ncols=True)
    for rep in data.values():
        pbar.update(1)
        tp = TreePair(*rep.split(" / "))
        L = tp.get_link()
        L = simplify_one_component(L)
        try:
            possible = L.get_knotinfo(unique=False)
        except NotImplementedError:
            knotinfo = "???"
        else:
            if all(sum(p.dict().values()) >= 2 for p in possible):
                # not prime
                continue
            knotinfo = "|".join(
                "*".join(sorted(Counter(p.dict()).elements()))
                for p in possible)
        leaves = rep.count("o") // 2
        results_by_num_leaves[leaves].append((knotinfo, rep))
        pbar.write(f"{knotinfo}: {rep}")
    pbar.close()
    with open(data_filepath.parent / "prime_knots_summary.md", "w") as f:
        print("# Prime Knots from Small Tree-Pairs", file=f)
        print(file=f)
        print(f"Generated using treeknots.get_knotinfo from {data_filepath.name}", file=f)
        print(file=f)
        print("Warning: in cases where homfly polynomials coincide, the results here may be incorrect.", file=f)
        print("In particular if `Link.get_knotinfo()` is more refined than the homfly polynomial,", file=f)
        print("then only a subset of the knots with that homfly polynomial will be listed here.", file=f)
        for num_leaves, arr in enumerate(results_by_num_leaves):
            if not arr:
                continue
            if num_leaves == 1:
                print(file=f)
                print("The unknot takes 1 leaf: `K0_1: o / o`", file=f)
                continue
            def sortkey(knotinfo):
                ki = knotinfo.split("|")[0].split("*")[0]
                a, b = ki.strip("Kmc").split("_")
                return (int(a.strip("n")), a, int(b), knotinfo)
            arr.sort(key=lambda kv: (float("inf"),) if kv[0] == "???" else sortkey(kv[0]))
            print(file=f)
            print("<details>", file=f)
            print(f"<summary>Prime Knots Requiring {num_leaves} {'leaf' if num_leaves == 1 else 'leaves'}</summary>", file=f)
            print(file=f)
            print("```", file=f)
            for knotinfo, rep in arr:
                print(f"{knotinfo}:\t{rep}", file=f)
            print("```", file=f)
            print("</details>", file=f)

if __name__ == "__main__":
    [path] = sys.argv[1:]
    prime_knots_summary(Path(path))
from collections import Counter, defaultdict
import gzip
from itertools import batched
from pathlib import Path
import sys

from tqdm import tqdm

from .tree_pair import TreePair

def get_name_from_wordpair(x: str):
    from sage.all import Link as sage_Link
    tp = TreePair(*x.split(" / "))
    [gc, _] = ogc = tp.oriented_gauss_code()
    L = sage_Link(ogc)

    import snappy
    sL = snappy.Link(L.pd_code())
    sL.simplify("global")
    L = sage_Link([[x+1 for x in tup] for tup in sL.PD_code()])
    if L.number_of_components() < len(gc):
        # Had unlinked unknots that got removed
        return None

    try:
        possible = L.get_knotinfo(unique=False)
    except NotImplementedError:
        return "???"
    if all(sum(p.dict().values()) >= 2 for p in possible):
        # not prime
        return None
    return "|".join(
        sorted(
            ["*".join(
                sorted(Counter(p.dict()).elements(),
                       key=lambda s: (int(s.split("_")[0].strip("Kn")),
                                      int(s.split("_")[1].strip("mrc")),
                                      s))
            ) for p in possible],
            key=lambda s: (int(s.split("*")[0].split("_")[0].strip("Kn")),
                            int(s.split("*")[0].split("_")[1].strip("mrc")),
                            s)
        )
    )

def get_all_names_from_worpairs(wordpairs):
    pbar = tqdm(total=len(wordpairs), dynamic_ncols=True)
    results = []
    for rep in wordpairs:
        knotinfo = get_name_from_wordpair(rep)
        pbar.update(1)
        if knotinfo is None:
            continue
        pbar.write(f"{knotinfo}: {rep}")
        results.append((knotinfo, rep))
    pbar.close()
    return results

def prime_knots_summary(data_filepath):
    with gzip.open(data_filepath, "rt", encoding="ascii") as f:
        lines = f.read().splitlines()
    data = ["o / o"] + [x for k, x in batched(lines, 2)]
    reps_by_knotinfo = defaultdict(list)
    for knotinfo, rep in get_all_names_from_worpairs(data):
        if knotinfo == "???":
            reps_by_knotinfo[f"???{rep}"].append(rep)
        else:
            knotinfo = "|".join(dict.fromkeys(
                p if "*" in p
                else p.replace("m", "").replace("r", "").replace("c","")
                for p in knotinfo.split("|")
            ))
            reps_by_knotinfo[knotinfo].append(rep)

    results_by_num_leaves = [[] for _ in range(30)]
    for knotinfo, reps in reps_by_knotinfo.items():
        rep = min(reps, key=len)
        leaves = rep.count("o") // 2
        results_by_num_leaves[leaves].append((knotinfo, rep))

    with open(data_filepath.parent / "prime_knots_summary.md", "w") as f:
        print("# Prime Knots from Small Tree-Pairs", file=f)
        print(file=f)
        print(f"Generated using treeknots.get_knotinfo from {data_filepath.name}", file=f)
        print(file=f)
        print("Warning: in cases where homfly polynomials coincide, the results here may be incorrect.", file=f)
        print("In particular if `Link.get_knotinfo()` is more refined than the homfly polynomial,", file=f)
        print("then only a subset of the knots with that homfly polynomial will be listed here.", file=f)
        print(file=f)

        for num_leaves, arr in enumerate(results_by_num_leaves):
            if not arr:
                continue
            if num_leaves == 1:
                print(file=f)
                print("The unknot takes 1 leaf: `K0_1: o / o`", file=f)
                continue
            def sortkey(knotinfo):
                if knotinfo.startswith("???"):
                    return (float("inf"),)
                ki = knotinfo.split("|")[0].split("*")[0]
                a, b = ki.strip("Kmcr").split("_")
                return (int(a.strip("n")), a, int(b), knotinfo)
            arr.sort(key=lambda kv: sortkey(kv[0]))
            print(file=f)
            drop_down = len(arr) > 10
            if drop_down:
                print("<details>", file=f)
                print(f"<summary>Prime Knots Requiring {num_leaves} leaves</summary>", file=f)
            else:
                print(f"Prime Knots Requiring {num_leaves} leaves:", file=f)
            print(file=f)
            print("```", file=f)
            for knotinfo, rep in arr:
                if knotinfo.startswith("???"):
                    print(f"???:\t{rep}", file=f)
                else:
                    print(f"{knotinfo}:\t{rep}", file=f)
            print("```", file=f)
            if drop_down:
                print("</details>", file=f)

if __name__ == "__main__":
    [path] = sys.argv[1:]
    path = Path(path)
    assert "sample" in path.name
    if "o" in path.name:
        prime_knots_summary(path)
    else:
        prime_links_summary(path)
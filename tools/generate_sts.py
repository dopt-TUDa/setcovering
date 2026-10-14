# This code produces the Steiner Triple System (STS) Set Covering instances.
# It is based on code written by Jasper van Doornmaalen and Christopher Hojny.
# https://github.com/JasperNL/scip-unified/blob/unified/problem_instances/covering_designs/generate.py
#
# Checking function produced by Claude.

import pyscipopt as ps
import pathlib
import numpy as np
import itertools
from collections import Counter
from itertools import combinations
from typing import List, Tuple


# Adpated from Sage code: https://github.com/sagemath/sage/blob/develop/src/sage/combinat/designs/bibd.py
#
#    if n % 6 == 3:
#        t = (n-3) // 6
#        Z = list(range(2 * t + 1))
#
#        T = lambda x_y : x_y[0] + (2*t+1)*x_y[1]
#
#        sts = [[(i,0),(i,1),(i,2)] for i in Z] + \
#            [[(i,k),(j,k),(((t+1)*(i+j)) % (2*t+1),(k+1) % 3)] for k in range(3) for i in Z for j in Z if i != j]
#    elif n % 6 == 1:
#
#        t = (n-1) // 6
#        N = list(range(2 * t))
#        T = lambda x_y : x_y[0]+x_y[1]*t*2 if x_y != (-1,-1) else n-1
#
#        L1 = lambda i,j : (i+j) % ((n-1)//3)
#        L = lambda i,j : L1(i,j)//2 if L1(i,j) % 2 == 0 else t+(L1(i,j)-1)//2
#
#        sts = [[(i,0),(i,1),(i,2)] for i in range(t)] + \
#            [[(-1,-1),(i,k),(i-t,(k+1) % 3)] for i in range(t,2*t) for k in [0,1,2]] + \
#            [[(i,k),(j,k),(L(i,j),(k+1) % 3)] for k in [0,1,2] for i in N for j in N if i < j]

def generate_triples(n: int) -> Tuple[List[tuple], List[tuple], list]:
    """
    Generating Steiner Triple Systems of size n
    """

    # construction for v = 6l + 3
    if ( n % 6 == 3 ):
        t = (n - 3) // 6
        Z = list(range(2 * t + 1))
        V = Z

        triples = []
        for i in Z:
            triples.append(((i,0),(i,1),(i,2)))

        for i in Z:
            for j in range(i):
                for k in range(3):
                     triples.append(((i,k), (j,k), ((t + 1) * (i + j) % (2 * t + 1), (k + 1) % 3)))

        points = [(v, k) for v in V for k in range(3)]
                     
    elif ( n % 6 == 1 ):
        t = (n - 1) // 6
        N = list(range(2 * t))
        V = N

        triples = []
        for i in range(t):
            triples.append(((i,0),(i,1),(i,2)))

        for i in range(t, 2 * t):
            for k in [0,1,2]:
                triples.append(((-1,-1), (i,k), (i - t, (k + 1) % 3)))

        L1 = lambda i,j: (i + j) % ((n-1) // 3)
        L = lambda i,j: L1(i,j) // 2 if L1(i,j) % 2 == 0 else t + (L1(i,j) - 1) // 2

        for i in N:
            for j in range(i):
                for k in [0,1,2]:
                    triples.append(((i,k), (j,k), (L(i,j),(k+1) % 3)))

        points = [(v, k) for v in V for k in range(3)] + [(-1, -1)]

    return triples, points, V


def verify_sts(triples: List[tuple], points: List[tuple]) -> bool:
    """
    Check that (points, triples) is a valid Steiner Triple System:
      1. every triple consists of 3 distinct points, all of which lie in `points`
      2. every pair of distinct points occurs in exactly one triple
      3. the number of triples matches n(n-1)/6, as required for an STS(n)
 
    Raises a ValueError describing the first violation found; returns True
    if the system passes the tests.
    """

    point_set = set(points)
    if len(point_set) != len(points):
        raise ValueError("Duplicate points found in the point set")
 
    pair_counts = Counter()
 
    for t_idx, triple in enumerate(triples):
        if len(set(triple)) != 3:
            raise ValueError(f"Triple {t_idx} does not consist of 3 distinct points: {triple}")
        for p in triple:
            if p not in point_set:
                raise ValueError(f"Triple {t_idx} contains point {p}, which is not in the point set")
        for pair in combinations(sorted(triple), 2):
            pair_counts[pair] += 1
 
    expected_pairs = set(combinations(sorted(point_set), 2))
    found_pairs = set(pair_counts.keys())
 
    missing = expected_pairs - found_pairs
    if missing:
        example = sorted(missing)[:5]
        raise ValueError(f"{len(missing)} pair(s) not covered by any triple, e.g. {example}")
 
    extra = found_pairs - expected_pairs
    if extra:
        example = sorted(extra)[:5]
        raise ValueError(f"{len(extra)} pair(s) found involving points outside the point set, e.g. {example}")
 
    bad = {pair: c for pair, c in pair_counts.items() if c != 1}
    if bad:
        example = list(bad.items())[:5]
        raise ValueError(f"{len(bad)} pair(s) covered more than once, e.g. {example}")
 
    n_points = len(point_set)
    expected_num_triples = n_points * (n_points - 1) // 6
    if len(triples) != expected_num_triples:
        raise ValueError(
            f"Number of triples ({len(triples)}) does not match the expected "
            f"n(n-1)/6 = {expected_num_triples} for n={n_points}"
        )
 
    return True


def generate_STS_model(n: int, verify: bool = True):
    """
    Generating the set covering problem for Steiner Triple Systems.
 
    If `verify` is True (default), the underlying triples are checked to
    actually form a Steiner Triple System before the ILP model is built.
    """
 
    model = ps.Model(problemName=f"STS{n}")
 
    triples, points, V = generate_triples(n)
 
    if verify:
        if verify_sts(triples, points):
            print('Test passed.')
        else:
             print('Test failed.')

    vars = []
    for v in V:
        row = []
        for i in range(3):
            row.append(model.addVar(name=f"x#{v}#{i}", vtype="B", obj=1.0))
        vars.append(row)
 
    if ( n % 6 == 1 ):
        specialvar = model.addVar(name=f"xx", vtype="B", obj=1.0)
 
    for i, S in enumerate(triples):
        model.addCons(name=f"set#{i}", cons=ps.quicksum(vars[x][i] if i != -1 else specialvar for j, (x,i) in enumerate(S)) >= 1)
 
    return model


if __name__ == "__main__":
    output_path = pathlib.Path(__file__).parent / pathlib.Path("instances/")

    output_path.mkdir(parents=True, exist_ok=True)

    for n in range(27,136):
        if ( n % 6 == 1 or n % 6 == 3 ):
            writepath = output_path / pathlib.Path(f"sts{n}.lp")

            model : ps.Model = generate_STS_model(n, True)
            model.writeProblem(writepath)

import numpy as np
import itertools



def id_info(id):
    sl = id.split('.')
    subset = sl[0]
    iscrm = int(sl[2])
    seqnum = int(sl[1])
    return subset, iscrm, seqnum


def validate_dpair(dpair, D):
    l = []
    for t in D:
        for s in D[t]:
            for pair in s:
                if pair not in dpair:
                    e = "Input file lacks some pairwise comparisons."
                    return e, False
                l.append('{}\t{}\t{}'.format(pair[0], pair[1], dpair[pair]))
    return False, "\n".join(l)

def benchmark(d, D):
    # d - user dict
    # D - ref CRM dict
    for t in D:
        l = []
        for pair in D[t][0]:
            l.append((d[pair], 1))
        for pair in D[t][1]:
            l.append((d[pair], 0))
        l.sort()
        # Establish n
        n = len(D[t][0])
        if n > 300: n = 300

        k = 0 
        s = []
        for i in range(n):
            k += l[i][1]
            s.append(str(k))
        yield t, k, n, k/n*100, ",".join(s)


def stats(perc_list, n_list):
    ps = np.array(perc_list, dtype='float')
    ns = np.array(n_list)
    average = np.mean(ps)
    if np.any(ns):
        weighted_average = np.average(ps, weights=ns)
    else:
        weighted_average = 0.0
    std = np.std(ps)
    return average, weighted_average, std

def ref_list(id_list):
    ids = id_list
    seq_num = len(ids)
    h = int(len(ids)/2)
    pos = ids[:h]
    neg = ids[h:]
    rl = [set([]), set([])]
    for i, l in enumerate([pos, neg]):
        for pair in itertools.combinations(l, 2):
            rl[i].add(tuple(sorted(pair)))
    return rl


def read_as_dict(json):
    D = {}
    for ref in json:
        D[ref] = ref_list(json[ref]["seqids"])
    return D
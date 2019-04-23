import numpy as np
from sklearn.metrics import roc_curve, auc

def validate_dpair(dpair, seqid_list):
    ID_D = {id: 0 for id in seqid_list}
    ID_IDX = {id: i for i, id in enumerate(seqid_list)}
    ID_N = len(ID_D)
    arr = np.zeros([ID_N, ID_N])
    l = []
    for pair in dpair:
        id1 = pair[0]
        id2 = pair[1]
        if id1 not in ID_D:
            e = "Unrecognized sequence id: {}".format(id1)
            return e, False
        if id2 not in ID_D:
            e = "Unrecognized sequence id: {}".format(id2)
            return e, False
        ID_D[id1] += 1
        ID_D[id2] += 1
        arr[ID_IDX[id1]][ID_IDX[id2]] = dpair[pair]
        arr[ID_IDX[id2]][ID_IDX[id1]] = dpair[pair]
        l.append('{}\t{}\t{}'.format(id1, id2, dpair[pair]))

    counts = set(ID_D.values())
    if len(counts) > 1 or list(counts)[0] != (ID_N-1):
        e = "Input file does not contain all pairwise comparisons."
        return e, False
    else:
        return False, "\n".join(l)


def dpair_to_dscore(dpair):
    max_value = max(dpair.values())
    for pair in dpair:
        dpair[pair] = max_value - dpair[pair]
    return dpair



def benchmark(dpair, REF_SEQ_IDS, strlevel_queryset, PROTEIN_ROC_MAX_POINTS):
    size_arr = len(dpair)
    d = dpair_to_dscore(dpair)
    for i, strlevel in enumerate(strlevel_queryset, start=1):
        if i == 1:
            arr_prediction = np.zeros(size_arr)
            arr_actual = np.zeros(size_arr, dtype=np.int)
        else:
            arr_actual.fill(0)
        levels = {id: REF_SEQ_IDS[id][:i] for id in REF_SEQ_IDS}
        for n, pair in enumerate(d):
            levels1 = levels[pair[0]]
            levels2 = levels[pair[1]]
            if levels1 == levels2:
                arr_actual[n] = 1
            if i == 1:
                arr_prediction[n] = d[pair]
        fpr, tpr, thresholds = roc_curve(arr_actual, arr_prediction)
        roc_auc = auc(fpr, tpr)

        # Limit number of points for ROC drawing.
        if tpr.size > PROTEIN_ROC_MAX_POINTS:
            tpr1 = np.zeros(PROTEIN_ROC_MAX_POINTS)
            fpr1 = np.zeros(PROTEIN_ROC_MAX_POINTS)
            step = int(tpr.size/PROTEIN_ROC_MAX_POINTS)
            for i in range(0, PROTEIN_ROC_MAX_POINTS):
                tpr1[i] = tpr[i*step]
                fpr1[i] = fpr[i*step]
            tpr1[-1] = tpr[-1]
            fpr1[-1] = fpr[-1]
            tpr = tpr1
            fpr = fpr1

        #generate an array with strings
        tpr = np.char.mod('%f', tpr)
        #combine to a string
        tpr_str = ",".join(tpr)
        #generate an array with strings
        fpr = np.char.mod('%f', fpr)
        #combine to a string
        fpr_str = ",".join(fpr)
        yield strlevel, roc_auc, tpr_str, fpr_str
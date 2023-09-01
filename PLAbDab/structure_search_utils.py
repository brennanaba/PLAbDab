import numpy as np
import numba as nb
from PLAbDab.region_definitions import *


from ImmuneBuilder import ABodyBuilder2

predictor = ABodyBuilder2()


def predict_antibody(seqs):
    ab = predictor.predict(seqs)

    coords = ab.aligned_traces[ab.ranking.index(0)].cpu().numpy()
    H_numb = [x[0][0] for x in ab.numbered_sequences["H"]]
    L_numb = [x[0][0] + 128 for x in ab.numbered_sequences["L"]]
    numbers = np.array(H_numb + L_numb, dtype = int)

    return numbers, coords


def get_antibody(text):
    lines = [x for x in text.split("\n") if x[13:15] == "CA"]
    size = len(lines)
    numbers = np.empty(size, dtype=int)
    coords = np.empty((size, 3))

    for i in range(size):
        line = lines[i]
        assert (line[21] == "H") or (line[21] == "L"), "Chains must be labelled H for heavy and L for light" 
        chain_term = 128 if line[21] == "L" else 0
        number = int(line[22:26])
        if number <= 128:
            numbers[i] = number + chain_term
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            coords[i] = (x, y, z)
        else:
            numbers[i] = -1

    return numbers[numbers!=-1], coords[numbers!=-1]


def parse_antibody(file):
    with open(file) as f:
        txt = f.read()
    return get_antibody(txt)


@nb.njit
def get_residues(antibody, selection):
    numbers, coords = antibody

    ids = np.zeros(len(numbers), dtype=nb.int32)
    for i, n in enumerate(numbers):
        if n in selection:
            ids[i] = 1

    select = np.zeros((sum(ids), 3))
    count = 0
    for i, val in enumerate(ids):
        if val == 1:
            select[count] = coords[i]
            count = count + 1

    return select


@nb.njit
def remove_insertions(ab):
    nums = ab[0]
    l_ab = len(nums)
    mask = np.ones(l_ab, np.int64)
    for i in range(1,l_ab):
        if nums[i] == nums[i-1]:
            mask[i] = 0
            
    indices = np.empty(sum(mask), np.int64)
    count = 0
    
    for i in range(l_ab):
        if mask[i] == 1:
            indices[count] = i
            count = count + 1
    return nums[indices], ab[1][indices]


def get_CDR_lengths(antibody):
    len_h1 = str(len(get_residues(antibody, reg_def["CDRH1"])))
    len_h2 = str(len(get_residues(antibody, reg_def["CDRH2"])))
    len_h3 = str(len(get_residues(antibody, reg_def["CDRH3"])))
    len_l1 = str(len(get_residues(antibody, reg_def["CDRL1"])))
    len_l2 = str(len(get_residues(antibody, reg_def["CDRL2"])))
    len_l3 = str(len(get_residues(antibody, reg_def["CDRL3"])))
    return len_h1, len_h2, len_h3, len_l1, len_l2, len_l3


@nb.njit
def get_alignment_transform(fixed, moveable, anchors):
    fixed, moveable = remove_insertions(fixed), remove_insertions(moveable)
    anchors = np.intersect1d(np.intersect1d(anchors, fixed[0]), moveable[0])
  
    anch1 = get_residues(moveable, anchors)
    anch2 = get_residues(fixed, anchors)
    
    n_residues = anch1.shape[0]
    
    anch1_center = anch1.sum(0) / n_residues
    anch2_center = anch2.sum(0) / n_residues
    
    anch1 = anch1-anch1_center
    anch2 = anch2-anch2_center

    V, _, W = np.linalg.svd(anch1.T @ anch2)
    U = V @ W

    if np.linalg.det(U) < 0:
        U = (np.array([[1, 1, -1]]) * V) @ W

    return anch1_center, anch2_center, U


@nb.njit
def align(fixed, moveable, anchors):
    x, y, U = get_alignment_transform(fixed, moveable, anchors)

    return moveable[0], ((moveable[1] - x) @ U + y)


@nb.njit
def rmsd(ab1, ab2, selection=reg_def_CDR_all, anchors=reg_def_fw_all):
    residues1 = get_residues(ab1, selection)
    residues2 = get_residues(align(ab1, ab2, anchors), selection)
    l = len(residues1)

    total = 0
    for i in range(l):
        total = total + sum((residues1[i] - residues2[i]) ** 2)

    return np.sqrt(total / l)

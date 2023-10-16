import numpy as np
import sys
from common.Constant import CODE

def calculate_MHL(minK:int, maxK:int, cpg_info: np.ndarray):
    maxLength = max(np.count_nonzero(cpg_info != 0, axis=1))
    maxK = min(maxK, maxLength)
    if (minK > maxK):
        print("Error: maxK is too large.\n")
        sys.exit(-1)

    meth_count_list = np.zeros(maxK, dtype='int')
    total_count_list = np.zeros(maxK, dtype='int')

    for i in range(minK - 1, maxK):
        for j in range(np.shape(cpg_info)[1] - i):
            cpg_info_split = cpg_info[..., j:(j + i + 1)]
            meth_count_list[i] += np.count_nonzero(np.count_nonzero(cpg_info_split == CODE.METHYLATED.value, axis=1) == (i + 1))
            total_count_list[i] += np.count_nonzero(np.count_nonzero(cpg_info_split != 0, axis=1) == (i + 1))

    weight = 0.0
    temp = 0.0
    for i in range(minK - 1, maxK):
        if total_count_list[i] > 0:
            temp += (i + 1.0) * meth_count_list[i] / total_count_list[i]
            weight += i + 1.0
    MHL = temp / weight
    return MHL

def calculate_MBS(K:int, cpg_info: np.ndarray):
    kmerNum = 0.0
    temp = 0.0
    for cpg_read in cpg_info:
        cpg_count = np.count_nonzero(cpg_read != 0)
        if cpg_count > K:
            cpg_str = ''.join(map(str, cpg_read)) # convert to string
            cpg_str_list = list(filter(None, cpg_str.split('0'))) # split by "0" and filter None value
            temp += sum([len(x) ** 2 for x in cpg_str_list]) / cpg_count ** 2
            kmerNum += 1.0

    MBS = temp / kmerNum
    return MBS

def calculate_Entropy(K:int, cpg_info: np.ndarray):
    kmer_dict = dict()
    for cpg_read in cpg_info:
        cpg_count = np.count_nonzero(cpg_read != 0)
        if cpg_count > K:
            cpg_str = ''.join(map(str, cpg_read))  # convert to string
            for i in range(len(cpg_str) - K + 1):
                kmer = cpg_str[i: i + K]
                if "0" not in kmer:
                    try:
                        kmer_dict[kmer] += 1
                    except:
                        kmer_dict[kmer] = 1

    temp = 0.0
    sums = sum(kmer_dict.values())
    for i in kmer_dict.values():
        temp += i / sums * np.log2(i / sums)
    Entropy = abs(-1 / 4 * temp)
    return Entropy

def calculate_R2(cpg_info: np.ndarray):
    pass
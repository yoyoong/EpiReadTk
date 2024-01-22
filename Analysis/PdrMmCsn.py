from scipy.sparse import csr_matrix
from scipy.stats import norm
import numpy as np
import pandas as pd
import argparse, sys, os
import matplotlib.pyplot as plt

def upperlower(data, boxsize):
    (n1, n2) = data.shape  # n1 gene; n2 sample
    upper = np.zeros((n1, n2), dtype=np.float32)
    lower = np.zeros((n1, n2), dtype=np.float32)

    for i in range(0, n1):
        s1 = sorted(data[i, :])
        s2 = data[i, :].argsort()
        sum = int(np.sum(np.sign(s1)))
        n3 = n2 - sum
        h = round(boxsize / 2 * sum)
        k = 0
        while k < n2:
            s = 0
            while k + s + 1 < n2:
                if s1[k + s + 1] == s1[k]:
                    s = s + 1
                else:
                    break
            if s >= h:
                upper[i, s2[k:k + s + 1]] = data[i, s2[k]]
                lower[i, s2[k:k + s + 1]] = data[i, s2[k]]
            else:
                upper[i, s2[k:k + s + 1]] = data[i, s2[min(n2 - 1, k + s + h)]]
                lower[i, s2[k:k + s + 1]] = data[i, s2[max(n3 * (n3 > h) + 1, k - h)]]

            k = k + s + 1
    return (upper, lower)

def getCSNList(input_data, boxSize, alpha, csn_writer):
    data, metric, regions = input_data.values, input_data.index.values, input_data.columns.values
    (n1, n2) = data.shape
    eps = np.finfo(float).eps

    (upper, lower) = upperlower(data, boxsize=boxSize)

    for k in range(0, n2):
        region = regions[k]
        MM = data[0][k]
        PDR = data[1][k]

        B = np.zeros((n1, n2), dtype=np.float32)
        for j in range(0, n2):
            B[:, j] = (data[:, j] <= upper[:, k]) & (data[:, j] >= lower[:, k]) & (data[:, k] > 0)
        a = B.sum(axis=1)
        a = np.reshape(a, (n1, 1))
        # temp = (np.dot(B, B.T) * n2 - np.dot(a, a.T)) / np.sqrt(np.dot(a, a.T) * np.dot(n2 - a, (n2 - a).T) / (n2 - 1) + eps)
        temp = (B @ B.T * n2 - a @ a.T) / np.sqrt((a @ a.T) * ((n2 - a) @ (n2 - a).T) / (n2 - 1) + eps)

        np.fill_diagonal(temp, 0)
        matrix = csr_matrix(temp).tocoo()

        if len(matrix.data) >= 1:
            if matrix.data[0] >= norm.ppf(1 - alpha):
                csn_writer.writelines(str(region) + "\t" + str(MM) + "\t" + str(PDR) + "\t" + "red" + "\n")
            else:
                csn_writer.writelines(str(region) + "\t" + str(MM) + "\t" + str(PDR) + "\t" + "blue" + "\n")
        else:
            csn_writer.writelines(str(region) + "\t" + str(MM) + "\t" + str(PDR) + "\t" + "blue" + "\n")

        print(str(region) + "end!")

def main(args):

    input_data = pd.read_csv(args.inputFile, sep='\t', usecols=[0, 1, 2, 3, 11, 12])
    input_data = input_data[input_data["MM"].notna() & input_data["PDR"].notna()]
    input_data['duplicated'] = input_data.duplicated(['pos'], keep=False)
    input_data = input_data[input_data['duplicated'] == True]

    color_file = f"{args.tag}.color.txt"
    csn_writer = open(color_file, "w")
    csn_writer.writelines("region\tMM\tPDR\tcolor\n")
    getCSNList(input_data, 0.1, 0.01, csn_writer)
    csn_writer.close()

    # color_file = "MM_PDR.color.txt"
    csn_result = pd.read_csv(color_file, sep='\t', index_col=0)

    red_x = csn_result[csn_result["color"] == "red"]["MM"]
    red_y = csn_result[csn_result["color"] == "red"]["PDR"]


    f = np.poly1d(np.polyfit(red_x, red_y, 5))
    t = np.linspace(0, 1, 250)

    for index, row in csn_result.iterrows():
        if (row.color == "blue") and (row.PDR < f(row.MM)):
            csn_result.loc[index, 'color'] = "green"

    blue_x = csn_result[csn_result["color"] == "blue"]["MM"]
    blue_y = csn_result[csn_result["color"] == "blue"]["PDR"]
    green_x = csn_result[csn_result["color"] == "green"]["MM"]
    green_y = csn_result[csn_result["color"] == "green"]["PDR"]

    csn_result.to_csv(color_file)

    fig = plt.gcf()
    plt.scatter(red_x, red_y, c="red", s=2)
    plt.scatter(blue_x, blue_y, c="blue", s=2)
    plt.scatter(green_x, green_y, c="green", s=2)
    plt.plot(t, f(t), '-')
    plt.title("MM-PDR CSN plot")
    plt.xlabel("MM")
    plt.ylabel("PDR")
    plt.figure(figsize=(10, 10))
    tag = args.tag
    fig.savefig(tag + '.pdf')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='This script requires the following dependencies: numpy/pandas/scipy/matplotlib')
    parser.add_argument('--inputFile', type=str, required=False, help="inputdata metrics file, generate by summaryBySNP",
                        default='/sibcb2/bioinformatics2/hongyuyang/code/EpiReadTk/outputDir/Z000000PX.csv')
    parser.add_argument('--boxSize', type=float, required=False, help="Size of neighborhood, Default = 0.1 (nx(k) = ny(k) = 0.1*n)", default=0.1)
    parser.add_argument('--alpha', type=float, required=False, help="Significant level (eg. 0.001, 0.01, 0.05 ...), Default = 0.01", default=0.01)
    parser.add_argument('--tag', type=str, required=False, help="prefix of the output file(s)", default="test")
    parser.set_defaults(func='csn')
    args = parser.parse_args()

    try:
        args.func
    except:
        sys.exit()

    main(args)
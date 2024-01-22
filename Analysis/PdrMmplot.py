import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

data_path = '/sibcb2/bioinformatics2/hongyuyang/code/EpiReadTk/outputDir/Z000000PX.csv'
data = pd.read_csv(data_path, delimiter="\t", header=0)
data['order'] = [1, 2] * (len(data) // 2)
data_1 = data[data['order'] == 1]
data_2 = data[data['order'] == 2]

# metrics = ['MM', 'PDR', 'CHALM', 'MHL', 'MCR', 'MBS', 'Entropy']
metrics = ['MM', 'PDR', 'CHALM', 'MHL', 'MCR', 'MBS', 'Entropy']
for metric in metrics:
    data_1_metric_mask = data_1[metric].apply(pd.isna)
    data_2_metric_mask = data_2[metric].apply(pd.isna)
    mask = np.array(~data_1_metric_mask) * np.array(~data_2_metric_mask)
    mask_nonone_num = np.count_nonzero(mask)
    data_1_metric = data_1[mask][metric]
    data_2_metric = data_2[mask][metric]

    plt.scatter(data_1_metric, data_2_metric, c='blue', s=3)
    plt.title(metric)
    plt.xlabel('SNP1')
    plt.ylabel('SNP2')
    plt.show()
    print(f'{metric} nonone_num: {mask_nonone_num}')

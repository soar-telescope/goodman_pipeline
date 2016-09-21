import matplotlib.pyplot as plt
import numpy as np

source = '/data/simon/data/SOAR/work/goodman/test/extraction-tests/results_cuhear.txt'

ff = open(source)
ff = ff.readlines()
header = []
data = []
for i in range(len(ff)):
    if ff[i][0] == '#':
        header.append(ff[i][0])
    else:
        data.append(ff[i].split())
   
print(header)
non_linear = []
linear = []
for line in data:
    ref = float(line[0])
    non = float(line[1])
    lin = float(line[2])
    non_linear.append((ref - non) ** 2)
    linear.append((ref - lin) ** 2)
    print(line)
print('RMS error non-linear : %s '%np.sqrt(np.mean(non_linear)))
print('RMS error linear : %s ' % np.sqrt(np.mean(linear)))


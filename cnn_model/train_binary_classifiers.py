from binary_classifiers import *
from dataset import Dataset
import numpy as np
import sys

idx = int(sys.argv[1]) - 1


x = {'type': 'flat', 'downsample': False, 'hapstats': False}
y = {'type': 'multiclass', 'mut': 'keep', 'sel': 'remove', 'inbreeding': 'add'}

# count number of classes
data = Dataset(x, y)
data.apply_feature_target_configs()
n_classes = data.targets.shape[1]
del data

# make index list
ijs = []
for i in range(n_classes):
    for j in range(i+1, n_classes):
        ijs.append((i,j))

# initialize classifier
ij = ijs[idx]
data = BinaryDataset(ij, x, y)
classifier = BinaryClassifier(ij, data)

with open(classifier.logpath + 'dataset_info.txt', 'wt') as txt:
    txt.write(f'x = {x}\n')
    txt.write(f'y = {y}\n')
    txt.write(f'x path: {data.feature_path}\n')
    txt.write(f'y path: {data.target_path}\n')

max_epochs = 300
classifier.fit_model(max_epochs, verbose=1, top_k=1)

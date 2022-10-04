from multiclass_classifiers import *
from dataset import Dataset
import numpy as np


x = {'type': 'flat', 'downsample': False, 'hapstats': False, 'tajD': False}
y = {'type': 'multiclass', 'mut': 'keep', 'sel': 'remove', 'inbreeding': 'add'}
data = Dataset(x, y)
multiclass = InbreedingMultiClassClassifier(data)

with open(multiclass.logpath + 'dataset_info.txt', 'wt') as txt:
    txt.write(f'x = {x}\n')
    txt.write(f'y = {y}\n')
    txt.write(f'x path: {data.feature_path}\n')
    txt.write(f'y path: {data.target_path}\n')

max_epochs = 5000
multiclass.fit_model(max_epochs, verbose=1)

multiclass.plot_loss()
multiclass.plot_acc()
multiclass.plot_confusion_matrix()

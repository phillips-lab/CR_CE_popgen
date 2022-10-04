from base import *
import numpy as np

# model building
from tensorflow.keras import Model
from tensorflow.keras.optimizers import Adadelta, Adam, Adagrad, Nadam
from tensorflow.keras.losses import CategoricalCrossentropy
from tensorflow.keras.layers import *


##############
#   MODELS   #
##############
out_type = 'multiclass'
class MultiClassClassifier(BaseClassifier):
    def __init__(self, dataset, vnum=None):
        super().__init__('multiclass', dataset, vnum)
        self.kernel_ws = [[15, 20], [6, 7]]
        self.trim = bool(self.dataset.downsample)

    def compile_model(self):
        w1, w2 = self.kernel_ws[self.trim]
        h = 4 if self.dataset.hapstats else 3
        x = Input(self.dataset.x_train[0,:].shape)

        y = Conv2D(6, (h, w2), activation='relu')(x)
        y = Dropout(0.4)(y)
        y = Conv2D(12, (h, w2), activation='relu')(y)
        y = Dropout(0.4)(y)
        y = Conv2D(24, (h, w1), activation='relu')(y)
        y = Dropout(0.4)(y)
        y = Conv2D(48, (2, w1), activation='relu')(y)
        y = Flatten()(y)
        y = Dense(1000, activation='relu')(y)
        y = Dropout(0.4)(y)
        y = Dense(200, activation='relu')(y)
        y = Dropout(0.4)(y)
        y = Dense(20, activation='relu')(y)
        y = Dense(self.dataset.y_train.shape[1], name='multiclass_targets')(y)

        model = Model(x, y)
        loss = CategoricalCrossentropy(from_logits=True)
        optmzr = Adam(learning_rate=0.000015)
        model.compile(optmzr, loss, metrics=['accuracy'])

        return model


class InbreedingMultiClassClassifier(BaseClassifier):
    def __init__(self, dataset, vnum=None):
        super().__init__('inbreeding_multiclass', dataset, vnum)
        self.kernel_ws = [[15, 20], [6, 7]]
        self.trim = bool(self.dataset.downsample)

    def compile_model(self):
        w1, w2 = self.kernel_ws[self.trim]
        h = 4 if self.dataset.hapstats else 3
        x = Input(self.dataset.x_train[0,:].shape)

        y = Conv2D(6, (h, w2), activation='relu')(x)
        y = Dropout(0.4)(y)
        y = Conv2D(12, (h, w2), activation='relu')(y)
        y = Dropout(0.4)(y)
        y = Conv2D(24, (h, w1), activation='relu')(y)
        y = Dropout(0.4)(y)
        y = Conv2D(48, (2, w1), activation='relu')(y)
        y = Flatten()(y)
        y = Dense(1000, activation='relu')(y)
        y = Dropout(0.4)(y)
        y = Dense(200, activation='relu')(y)
        y = Dropout(0.4)(y)
        y = Dense(20, activation='relu')(y)
        y = Dense(self.dataset.y_train.shape[1], name='multiclass_targets')(y)

        model = Model(x, y)
        loss = CategoricalCrossentropy(from_logits=True)
        optmzr = Adam(learning_rate=0.000015)
        model.compile(optmzr, loss, metrics=['accuracy'])

        return model


if __name__ == '__main__':
    print('success')

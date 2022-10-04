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
out_type = 'multitarget'
class OneChannelClassifier(BaseClassifier):
    def __init__(self, in_type='flat', balancing='keep', mutation='keep'):
        super().__init__('one_channel', in_type, out_type, balancing, mutation)

    def compile_model(self):
        # build network
        x = Input(self.x_train[0,:].shape)

        y = Conv2D(6, (4, 15), activation='relu')(x)
        y = Dropout(0.4)(y)
        y = Conv2D(12, (4, 20), activation='relu')(y)
        y = Dropout(0.4)(y)
        y = Conv2D(24, (4, 15), activation='relu')(y)
        y = Dropout(0.4)(y)
        y = Conv2D(48, (2, 15), activation='relu')(y)
        y = Flatten()(y)
        y = Dense(1000, 'relu')(y)
        y = Dropout(0.4)(y)
        y = Dense(200, 'relu')(y)
        y = Dropout(0.4)(y)
        y = Dense(20, 'relu')(y)

        # concatenate last hidden layers
        mut_pred = Dense(self.mut_len, name='mutation')(y)
        slf_pred = Dense(self.slf_len, name='selfing')(y)
        sel_pred = Dense(self.sel_len, name='selection')(y)

        model = Model(x, [mut_pred, slf_pred, sel_pred])

        # losses
        losses = [CategoricalCrossentropy(from_logits=True), 
                  CategoricalCrossentropy(from_logits=True),
                  CategoricalCrossentropy(from_logits=True)]

        optmzr = Adam(learning_rate=0.00001)
        model.compile(optmzr, losses, loss_weights=[1.0, 1.0, 1.0],
                      metrics=['accuracy'])
        return model


class Conv3DClassifier(BaseClassifier):
    def __init__(self, in_type='stack', balancing='keep', mutation='keep'):
        super().__init__('conv3d', in_type, out_type, balancing, mutation)
        self.expand_fvec_dims()

    def expand_fvec_dims(self):
        """ Adds a dimension to fvecs for Conv3D layers """
        self.x_train = np.expand_dims(self.x_train, 4)
        self.x_test = np.expand_dims(self.x_test, 4)
        assert self.x_train.shape[1:] == (13, 25, 3, 1)

    def compile_model(self):
        # build network
        x = Input(self.x_train[0,:].shape)

        y = Conv3D(6, (4, 10, 2), activation='relu')(x)
        y = Dropout(0.4)(y)
        y = Conv3D(12, (4, 4, 1), activation='relu')(y)
        y = Dropout(0.4)(y)
        y = Conv3D(24, (4, 4, 1), activation='relu')(y)
        y = Dropout(0.4)(y)
        y = Conv3D(48, (4, 4, 1), activation='relu')(y)
        y = Flatten()(y)
        y = Dense(1000, 'relu')(y)
        y = Dropout(0.4)(y)
        y = Dense(200, 'relu')(y)
        y = Dropout(0.4)(y)
        y = Dense(20, 'relu')(y)

        # concatenate last hidden layers
        mut_pred = Dense(self.mut_len, name='mutation')(y)
        slf_pred = Dense(self.slf_len, name='selfing')(y)
        sel_pred = Dense(self.sel_len, name='selection')(y)

        model = Model(x, [mut_pred, slf_pred, sel_pred])

        # losses
        losses = [CategoricalCrossentropy(from_logits=True), 
                  CategoricalCrossentropy(from_logits=True),
                  CategoricalCrossentropy(from_logits=True)]

        optmzr = Adam(learning_rate=0.00001)
        model.compile(optmzr, losses, loss_weights=[1.0, 1.0, 1.0],
                      metrics=['accuracy'])
        return model

    
class VanillaClassifier(BaseClassifier):
    def __init__(self, in_type='stack', balancing='keep', mutation='keep'):
        super().__init__('vanilla', in_type, out_type, balancing, mutation)

    def compile_model(self):
        # build network
        x = Input(self.x_train[0,:].shape)

        y = Conv2D(6, (4, 5), activation='relu')(x)
        y = Dropout(0.4)(y)
        y = Conv2D(12, (4, 7), activation='relu')(y)
        y = Dropout(0.4)(y)
        y = Conv2D(24, (4, 5), activation='relu')(y)
        y = Dropout(0.4)(y)
        y = Conv2D(48, (2, 4), activation='relu')(y)
        y = Flatten()(y)
        y = Dense(1000, 'relu')(y)
        y = Dropout(0.4)(y)
        y = Dense(200, 'relu')(y)
        y = Dropout(0.4)(y)
        y = Dense(20, 'relu')(y)

        # concatenate last hidden layers
        mut_pred = Dense(self.mut_len, name='mutation')(y)
        slf_pred = Dense(self.slf_len, name='selfing')(y)
        sel_pred = Dense(self.sel_len, name='selection')(y)

        model = Model(x, [mut_pred, slf_pred, sel_pred])

        # losses
        losses = [CategoricalCrossentropy(from_logits=True), 
                  CategoricalCrossentropy(from_logits=True),
                  CategoricalCrossentropy(from_logits=True)]

        optmzr = Adam(learning_rate=0.000015)
        model.compile(optmzr, losses, loss_weights=[1.0, 1.0, 1.0],
                      metrics=['accuracy'])
        return model


class ParallelClassifier(BaseClassifier):
    def __init__(self, balancing='keep'):
        super().__init__('parallel', balancing)

    def compile_model(self):
        # build network
        x = Input(self.x_train[0,:].shape)

        def build_branch(x):
            y = Conv2D(6, (2, 5), activation='relu')(x)
            y = Dropout(0.2)(y)
            y = Conv2D(12, (3, 5), activation='relu')(y)
            y = Conv2D(24, (4, 5), activation='relu')(y)
            y = Dropout(0.2)(y)
            y = Conv2D(48, (5, 5), activation='relu')(y)
            y = Flatten()(y)
            y = Dense(100, 'relu')(y)
            y = Dropout(0.1)(y)
            y = Dense(20, 'relu')(y)
            return y

        # network branches
        m = build_branch(x) # mutation branch
        s = build_branch(x) # selfing branch
        sel = build_branch(x) # selection branch

        # concatenate last hidden layers
        mss = concatenate([m, s, sel])
        mut_pred = Dense(self.mut_len, 'softmax', name='mutation')(mss)
        slf_pred = Dense(self.slf_len, 'softmax', name='selfing')(mss)
        sel_pred = Dense(self.sel_len, 'softmax', name='selection')(mss)

        model = Model(x, [mut_pred, slf_pred, sel_pred])

        # losses
        losses = [CategoricalCrossentropy(from_logits=True), 
                  CategoricalCrossentropy(from_logits=True),
                  CategoricalCrossentropy(from_logits=True)]

        # loss weights for old data: [0.2, 1.0, 1.0]
        model.compile('adam', losses, loss_weights=[1.0, 1.0, 1.0],
                      metrics=['accuracy'])
        return model


class TowerClassifier(BaseClassifier):
    def __init__(self, balancing='keep', mutation='keep'):
        super().__init__('tower', balancing, mutation)

    def compile_model(self):
        # build network
        x = Input(self.x_train[0,:].shape)

        def build_tower(x):
            y = Conv2D(6, (2, 5), activation='relu')(x)
            y = Dropout(0.3)(y)
            y = Conv2D(12, (3, 5), activation='relu')(y)
            y = Dropout(0.3)(y)
            y = Conv2D(24, (4, 5), activation='relu')(y)
            y = Dropout(0.3)(y)
            y = Conv2D(48, (5, 5), activation='relu')(y)
            y = Flatten()(y)
            return y

        def build_branch(x):
            y = Dense(100, 'relu')(x)
            y = Dropout(0.3)(y)
            y = Dense(20, 'relu')(y)
            return y

        def build_network(x):
            # towers
            towers = [build_tower(x) for _ in range(3)]

            # tower branches
            branches = [[build_branch(t) for _ in range(3)] for t in towers]

            # connect to losses
            mut, slf, sel = [], [], []
            for b in branches: # loops through mut, slf, sel tower branches
                mut.append(b[0])
                slf.append(b[1])
                sel.append(b[2])

            # concatenate last layers
            muts = concatenate(mut)
            slfs = concatenate(slf)
            sels = concatenate(sel)
            mut_pred = Dense(self.mut_len,
                             'softmax',
                             name='mutation')(muts)
            slf_pred = Dense(self.slf_len,
                             'softmax',
                             name='selfing')(slfs)
            sel_pred = Dense(self.sel_len,
                             'softmax',
                             name='selection')(sels)

            return Model(x, [mut_pred, slf_pred, sel_pred])

        model = build_network(x)

        # losses
        losses = [CategoricalCrossentropy(from_logits=True), 
                  CategoricalCrossentropy(from_logits=True),
                  CategoricalCrossentropy(from_logits=True)]

        # loss weights for old data: [0.2, 1.0, 1.0]
        model.compile('adam', losses, loss_weights=[0.5, 0.8, 1.0],
                      metrics=['accuracy'])
        return model


class BounceClassifier(BaseClassifier):
    def __init__(self, balancing='keep'):
        super().__init__('bounce', balancing)

    def compile_model(self):
        # input template
        x = Input(self.x_train[0,:].shape)

        def conv_branch(x, channels, kernel_size):
            y = Conv2D(channels, kernel_size, activation='relu')(x)
            return y

        def dense_branch(x, units):
            y = Dense(units, activation='relu')(x)
            return y

        def conv_bounce(x, channels, kernel_size):
            y1 = conv_branch(x, channels, kernel_size)
            y2 = conv_branch(x, channels, kernel_size)
            y3 = conv_branch(x, channels, kernel_size)
            y = concatenate([y1, y2, y3])
            return y

        def dense_bounce(x, units, flatten=False):
            if flatten:
                x = Flatten()(x)
            y1 = dense_branch(x, units)
            y2 = dense_branch(x, units)
            y3 = dense_branch(x, units)
            y = concatenate([y1, y2, y3])
            return y

        # build network
        y = conv_bounce(x, 6, (2, 5))
        y = Dropout(0.2)(y)
        y = conv_bounce(y, 12, (3, 5))
        y = Dropout(0.2)(y)
        y = conv_bounce(y, 24, (4, 5))
        y = Dropout(0.2)(y)
        y = conv_bounce(y, 48, (5, 5))
        y = dense_bounce(y, 500, flatten=True)
        y = Dropout(0.1)(y)
        y = dense_bounce(y, 100)
        y = Dropout(0.1)(y)
        y = dense_bounce(y, 20)

        # final layer predic://kr-colab.slack.com/archives/DF99RAKU4tions
        mut_pred = Dense(self.mut_len, 'softmax', name='mutation')(y)
        slf_pred = Dense(self.slf_len, 'softmax', name='selfing')(y)
        sel_pred = Dense(self.sel_len, 'softmax', name='selection')(y)

        model = Model(x, [mut_pred, slf_pred, sel_pred])

        # losses
        losses = [CategoricalCrossentropy(from_logits=True), 
                  CategoricalCrossentropy(from_logits=True),
                  CategoricalCrossentropy(from_logits=True)]

        # loss weights for old data: [0.2, 1.0, 1.0]
        model.compile('adam', losses, loss_weights=[1.0, 1.0, 1.0],
                      metrics=['accuracy'])
        return model

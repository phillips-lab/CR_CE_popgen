from base import *
from dataset import Dataset
import numpy as np
from glob import glob

# keras cnn model 
from tensorflow.keras.layers import Input, Conv2D, Dropout, Dense, Flatten
from tensorflow.keras.losses import BinaryCrossentropy
from tensorflow.keras.optimizers import Adam
from tensorflow.keras import Model


class BinaryDataset(Dataset):
    def __init__(self, class_idxs, x, y):
        self._ij = class_idxs
        super().__init__(x, y)
        
    def prepare_for_training(self, test_size=0.1):
        self.apply_feature_target_configs()
        self.slice_ij_data()
        self.train_test_split(test_size)
        
    def slice_ij_data(self):
        class_idxs = self._ij
        targets = self.targets
        sim_idxs = [i for i, y in enumerate(targets) if np.argmax(y) in class_idxs]
        self.targets = targets[sim_idxs][:, class_idxs]
        self.features = self.features[sim_idxs]
        
# subclass from base and overload some methods
class BinaryClassifier(BaseClassifier):
    def __init__(self, class_idxs, dataset, vnum=None):
        self._ij = class_idxs
        self._i, self._j = class_idxs
        super().__init__('binary_ensemble', dataset, vnum=0)
        
    def compile_model(self):
        h, w1, w2 = 3, 15, 20 # kernel shapes
        in_shape = self.dataset.x_train[0].shape
        x = Input(in_shape)
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
        y = Dense(2, activation='softmax', name='targets')(y)
        
        model = Model(x, y)
        loss = BinaryCrossentropy(from_logits=False)
        optmzr = Adam(learning_rate=0.000015)
        model.compile(optmzr, loss, metrics=['accuracy'])

        return model
    
    #### OVERLOADED METHODS: slice data?, logpath stuff
    ### MAYBE UPDATE THIS: new_logpath = logpath + 'i,j_'
    def set_logpath(self, model_name, vnum=None):
        if vnum is None:
            self.set_logpath_(model_name)
        else:
            logpath = f'model_logs/{model_name}/version_{vnum}/'
            assert exists(logpath)
            self.vnum = vnum
            self.logpath = logpath + f'{self._i},{self._j}_'
            
    ### INCLUDE INDEX i,j IN MODEL CHECKPOINT NAME (and maybe also loss)
    def get_callbacks(self, max_epochs):
        """ configures checkpoint and lr scheduler callbacks """
        # checkpoint callback
        fname = self.logpath + 'epoch={epoch:02d}'
        if self.out_type == 'multitarget':
            fname += '-mut_acc={val_mutation_accuracy:.3f}'
            fname += '-selec_acc={val_selection_accuracy:.3f}'
            fname += '-self_acc={val_selfing_accuracy:.3f}.hdf5'
            monitor = 'mutation_accuracy'
        
        elif self.out_type == 'multiclass':
            fname += '-loss={loss:.3f}-val_acc={val_accuracy:.3f}.hdf5'
            monitor = 'loss'
        ckpt = ModelCheckpoint(filepath=fname,
                               monitor=monitor,
                               save_best_only=True,
                               mode='min')
        prog = PrintProgress(max_epochs)
        stop = EarlyStopping('val_loss', min_delta=1e-4, patience=30)

        return [ckpt, prog, stop]
    
    ### MIGHT NEED TO CORRECT SPLITS FOR i,j
    def get_ckpt_accs(self):
        ckpt_filepaths = glob(self.logpath + '*hdf5')
        ckpts = [c.split('/')[-1] for c in ckpt_filepaths]

        if self.out_type == 'multitarget':
            # get checkpoint selection accuracy
            accs = [c.split('-')[2] for c in ckpts]
            accs = [float(s.split('=')[-1]) for s in accs]
        elif self.out_type == 'multiclass':
            accs = [c.split('=')[3] for c in ckpts]
            accs = [float(f'{c.split(".")[0]}.{c.split(".")[1]}') for c in accs]
        return accs
        



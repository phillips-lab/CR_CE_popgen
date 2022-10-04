import numpy as np
import pandas as pd
import pickle
from sklearn.model_selection import train_test_split
from sklearn.utils import resample, shuffle
from os.path import exists


n_trim, n_rep = 3, 3
features_map = {'type': ('flat', 'stack'), # shape of feature vectors
                'downsample': (False, (n_trim, n_rep)), # either False or tuple
                'hapstats': (True, False)}

targets_map = {'type': ('multitarget', 'multiclass', 'multilabel'), 
               'mut': ('remove', 'keep'),
               'sel': ('remove', 'resample', 'keep'),
               'inbreeding': ('replace', 'add')}

# wrangled .pkl files --> train/test features/targets & empirical features
class Dataset:
    def __init__(self, features, targets):
        self.feature_config = features
        self.target_config = targets

        self.downsample = features['downsample']
        self.in_type = features['type']
        self.out_type = targets['type']

        # use default mut (keep) and sel (remove)
        if 'mut' not in targets.keys():
            self.mut = 'keep'
        else:
            self.mut = targets['mut']
        if 'sel' not in targets.keys():
            self.sel = 'remove'
        else:
            self.sel = targets['sel']
        
        assert self.in_type in ['stack', 'flat']
        assert self.out_type in ['multitarget', 'multiclass', 'multilabel']

        # remove haplotype stats nDiplo, H1, H12, H2/H1
        if 'hapstats' in features.keys():
            self.hapstats = features['hapstats']
        else:
            self.hapstats = True
        if 'tajD' in features.keys():
            self.tajd = features['tajD']
        else:
            self.tajd = True

        # use inbreeding for outcrossing
        if 'inbreeding' in targets.keys():
            self.inbreeding = targets['inbreeding']
        else:
            self.inbreeding = False

        self.data_dir = '/projects/haldane/mlukac/worm-landscapes/data/'
        self.init_target_names()
        if self.downsample:
            self.data_dir += 'trim/'
        self.init_feature_path()
        self.init_target_path()
        self.var_keys = ['mutation', 'selfing', 'selection']
        
    ##################
    ## INITIALIZING ##
    ##################
    def init_feature_path(self):
        if self.downsample:
            n_trim, n_rep = self.downsample
            self.feature_path = f'fvecs_{self.in_type}_trim={n_trim}_rep={n_rep}'
        else:
            self.feature_path = f'fvecs_{self.in_type}'
        if not exists(self.data_dir + f'{self.feature_path}.pkl'):
            print('Desired features do not exist. To make them, use Dataset.create().')
    
    def init_target_path(self):
        if self.downsample:
            n_rep = self.downsample[-1]
            self.target_path = f'{self.out_type}_targets_rep={n_rep}'
        else:
            self.target_path = f'{self.out_type}_targets'
        if not exists(self.data_dir + f'{self.target_path}.pkl'):
            print('Desired targets do not exist. To make them, use Dataset.create().')
            
    def add_inbreeding_target_names(self):
        assert hasattr(self, '_target_names')

        # no balancing selection in inbreeding sims
        names = self._target_names
        ii = [i for i, n in enumerate(names) if ', 0.0,' not in n or 'NDBa' in n]
        inbr_names = np.delete(names, ii)

    def init_target_names(self):
        names = self.load_data(f'{self.out_type}_names')
        self._target_names = names # unaltered target names
        self.target_names = names

 #       if self.inbreeding == 'add':
 #           self.add_inbreeding_target_names()
 #           for mut in ['1.0', '1.15', '1.5', '2.0']:
 #               for sel in ['N', 'ND', 'NDB', 'NDBa']:
 #                   self._target_names.append(f'({mut}, 0.0i, {sel})')

 #       # remove appropriate names if necessary
 #       if self.mut == 'remove':
 #           idxs = [i for i, n in enumerate(names) if '(2.0,' in n]
 #           names = np.delete(names, idxs)
 #       if self.sel == 'remove':
 #           idxs = [i for i, n in enumerate(names) if 'NDBa' in n]
 #           names = np.delete(names, idxs)

 #       self.target_names = names.tolist()
            
    def init_features(self):
        self.features = self.load_data(f'{self.feature_path}')
    
    def init_targets(self):
        self.targets = self.load_data(f'{self.target_path}')

    def create(self):
        '''Use Wrangler to downsample'''
        raise NotImplementedError()
    
    ##############
    ## PREPPING ##
    ##############
    def apply_feature_target_configs(self):
        # load data
        self.init_features()
        self.init_targets()
        
        # handle inbreeding first
        if self.inbreeding == 'replace':
            self._replace_outcrossing_with_inbreeding()
        elif self.inbreeding == 'add':
            self._add_inbreeding()

        # feature transforms
        if not self.hapstats:
            self._remove_hapstats()

        # remove some simulations
        if self.mut == 'remove':
            self.remove_2x_mutation_rate()
        if self.sel == 'remove':
            self._remove_balancing_selection()
            
    def prepare_for_training(self, test_size=0.3):
        self.apply_feature_target_configs()
        self.train_test_split(test_size)

    #############
    ## SETTING ##
    #############
    def set_features(self, features):
        self.features = features
        
    def set_targets(self, targets):
        self.targets = targets
            
    #############
    ## GETTING ##
    #############
    def load_data(self, fname):
        filename = f'{self.data_dir}{fname}.pkl'
        with open(filename, 'rb') as f:
            data = pickle.load(f)
        return data
    
    def get_empirical_features(self, worm, chrom, kb_indiv=None):
        # infer number of subwindows per domain
        if self.downsample:
            n_subw = 11
        else:
            n_subw = 25

        # load empirical features
        emp_path = f'data/empirical/{n_subw}subw/'
        if kb_indiv is None:
            emp_path += f'fvecs_{worm}_chr_{chrom}.pkl'
        else:
            kb, indiv = kb_indiv
            emp_path += f'fvecs_{worm}{kb}{indiv}_chr_{chrom}.pkl'
        with open(emp_path, 'rb') as f:
            emp_features = pickle.load(f)

        # flatten if necessary
        if self.in_type == 'flat':
            n1, n2, n3, n4 = emp_features.shape
            flat_emp_features = np.zeros((n1, n2, n3 * n4, 1))
            for i, x in enumerate(emp_features):
                flat_emp_features[i] = np.reshape(x, (n2, n3*n4, 1), 'F')
            emp_features = flat_emp_features

        if not self.hapstats:
            emp_features = self.remove_hapstats(emp_features)

        return emp_features

    def get_multitarget_path(self):
        tp = self.target_path.split('_')
        tp[0] = 'multitarget'
        multitarget_path = '_'.join(tp)
        return multitarget_path
    
    ##########################
    ## TRAIN/TEST SPLITTING ##
    ##########################
    def train_test_split(self, test_size):
        if self.out_type == 'multitarget':
            self.multitarget_train_test_split(test_size)
        if self.out_type == 'multiclass':
            self.multiclass_train_test_split(test_size)
            
    def multitarget_train_test_split(self, test_size):
        fvecs = self.features
        mutation, selfing, selection = self.targets.values()
        self.mut_len = mutation.shape[1]
        self.slf_len = selfing.shape[1]
        self.sel_len = selection.shape[1]
        split_data = train_test_split(fvecs,
                                      mutation,
                                      selfing,
                                      selection,
                                      test_size=test_size,
                                      random_state=23)
        if self.sel == 'resamp':
            split_data = self.multitarget_resamp_balancing_selection(split_data)
        # unpack split data
        (self.x_train, self.x_test,
         self.mut_train, self.mut_test,
         self.slf_train, self.slf_test,
         self.sel_train, self.sel_test) = split_data
        self.y_train = [self.mut_train, self.slf_train, self.sel_train]
        self.y_test = [self.mut_test, self.slf_test, self.sel_test]
            
    def multiclass_train_test_split(self, test_size):
        fvecs = self.features
        targets = self.targets
        self.target_len = targets.shape[1]
        split_data = train_test_split(fvecs, targets,
                                      test_size=test_size,
                                      random_state=23)
        if self.sel == 'resamp':
            self.multiclass_resamp_balancing_selection(split_data)
        self.x_train, self.x_test, self.y_train, self.y_test = split_data

    ##
    ## MERGE AND SHUFFLE
    ##
    def merge_data(self, fvecs, targs):
        '''
        Input:
            fvecs (tuple): contains feature arrays to merge
            targs (tuple): contains target arrays to merge
        Returns: merged features and targets
        '''
        fvecs = np.concatenate(fvecs)
        targs = np.concatenate(targs)
        fvecs, targs = shuffle(fvecs, targs, random_state=23)
        return fvecs, targs

    ##########################
    ##########################
    ## DATA TRANSFORMATIONS ##
    ##########################
    ##########################

    ##
    ## REMOVE HAPLOTYPE STATISTICS
    ##
    def remove_hapstats(self, features):
        if self.tajd:
            hapstat_idxs = [6,7,8,9]
        else:
            hapstat_idxs = [2,6,7,8,9]
        features = np.delete(features, hapstat_idxs, axis=1)
        return features

    def _remove_hapstats(self):
        fvecs = self.features
        fvecs = self.remove_hapstats(fvecs)
        self.set_features(fvecs)

    ##
    ## REMOVE BALANCING SELECTION
    ##
    def remove_balancing_selection(self):
        if self.out_type == 'multitarget':
            self.multitarget_remove_balancing_selection()
        elif self.out_type == 'multiclass':
            self.multiclass_remove_balancing_selection()
            
    def remove_balancing_selection_(self, features, targets):
        if self.out_type == 'multiclass':
            names = self.target_names
            assert len(names) == targets.shape[1]
            cols = [i for i, n in enumerate(names) if 'NDBa' in n]
            rows = [i for i, y in enumerate(targets) if np.argmax(y) in cols]
            self.target_names = np.delete(names, cols).tolist()
            targets = np.delete(targets, cols, axis=1)
            targets = np.delete(targets, rows, axis=0)
            features = np.delete(features, rows, axis=0)
        elif self.out_type == 'multitarget':
            raise NotImplementedError()
        else:
            raise ValueError('wrong out_type specified')
        return features, targets

    def _remove_balancing_selection(self):
        fvecs, targs = self.features, self.targets
        fvecs, targs = self.remove_balancing_selection_(fvecs, targs)
        self.set_features(fvecs)
        self.set_targets(targs)

    def multitarget_remove_balancing_selection(self):
        fvecs = self.features
        targets = self.targets
        self._target_names['selection'] = self._target_names['selection'][:-1]
        mutation, selfing, selection = targets.values()
        not_bal = selection[:,-1]==0 # rows without balancing
        # remove balancing rows (and column for sel)
        fvecs = fvecs[not_bal]
        mutation = mutation[not_bal]
        selfing = selfing[not_bal]
        selection = selection[not_bal,:-1]
        targets = dict(zip(self.var_keys, [mutation, selfing, selection]))
        
        self.set_features(fvecs)
        self.set_targets(targets)
    
    def multiclass_remove_balancing_selection(self):
        # load multitarget targets for slicing
        fvecs = self.features
        targets = self.targets
        names = self.target_names

        # balancing columns 
        bal_cols = [i for i in range(3, targets.shape[1]+1, 4)] 
        targets = np.delete(targets, bal_cols, axis=1)
    
        out_cols = [i for i, n in enumerate(names) if ', 0.0, ' in n]
        out_rows = [i for i, y in enumerate(targs) if np.argmax(y) in out_cols]

        self.features = np.delete(fvecs, out_rows, axis=0)
        self.targets = np.delete(targs, out_rows, axis=0)

    ##
    ## REPLACE OUTCROSSING WITH INBREEDING
    ##
    def replace_outcrossing_with_inbreeding(self, fvecs, targs):
        # remove outcrossing rows
        if self.out_type == 'multiclass':
            names = self._target_names
            assert len(names) == targs.shape[1]
            cols = [i for i, n in enumerate(names) if ', 0.0, ' in n]
            rows = [i for i, y in enumerate(targs) if np.argmax(y) in cols]
            fvecs = np.delete(fvecs, rows, axis=0)
            targs = np.delete(targs, rows, axis=0)
        elif self.out_type == 'multitarget':
            raise NotImplementedError()
        else:
            raise ValueError('wrong out_type specified')

        # merge inbreeding and shuffle 
        inb_fvecs = self.load_data(f'inbreeding/fvecs_{self.in_type}')
        inb_targs = self.load_data(f'inbreeding/{self.out_type}_targets')
        if not self.hapstats:
            inb_fvecs = self.remove_hapstats(inb_fvecs)
        fvecs, targs = self.merge_data((fvecs, inb_fvecs), (targs, inb_targs))

        return fvecs, targs

    def _replace_outcrossing_with_inbreeding(self):
        fvecs, targs = self.features, self.targets
        fvecs, targs = self.replace_outcrossing_with_inbreeding(fvecs, targs)
        self.set_features(fvecs)
        self.set_targets(targs)

    ##
    ## ADD INBREEDING
    ##
    def _add_inbreeding(self):
        if self.out_type == 'multitarget':
            self.multitarget_add_inbreeding()
        elif self.out_type == 'multiclass':
            self._multiclass_add_inbreeding()

    def _multiclass_add_inbreeding(self):
        names, fvecs, targs = self.multiclass_add_inbreeding()
        fvecs, targs = shuffle(fvecs, targs, random_state=23)
        self.target_names = names
        self.features = fvecs
        self.targets = targs

    def multiclass_add_inbreeding(self):
        # load features and targets
        fvecs = self.features
        targs = self.targets
        names = self.target_names
        inbr_fvecs = self.load_data(f'inbreeding/fvecs_{self.in_type}')
        inbr_targs = self.load_data(f'inbreeding/{self.out_type}_targets')
        assert len(self.target_names) == inbr_targs.shape[1]
        assert targs.shape[1] == inbr_targs.shape[1]

        # remove zero columns from inbreeding one-hots
        ii = [i for i, n in enumerate(names) if ', 0.0,' not in n or 'NDBa' in n]
        inbr_targs = np.delete(inbr_targs, ii, axis=1)

        # update inbreeding names
        inbr_names = np.delete(names, ii).tolist()
        for i, n in enumerate(inbr_names):
            n = n.split(', ')
            n[1] = '0.0i'
            n = ', '.join(n)
            inbr_names[i] = n

        # get zeros arrays to pad one-hots
        zs = np.zeros((len(targs), len(inbr_names)))
        inbr_zs = np.zeros((len(inbr_targs), len(names)))
        targs = np.hstack((targs, zs))
        inbr_targs = np.hstack((inbr_zs, inbr_targs))
        assert targs.shape[1] == inbr_targs.shape[1]

        # stack features, targets, names
        fvecs = np.vstack((fvecs, inbr_fvecs))
        targs = np.vstack((targs, inbr_targs))
        names = names + inbr_names

        return names, fvecs, targs
        
    def multitarget_add_inbreeding(self):
        pass
    
    def include_inbreeding(self):
        assert hasattr(self, 'features')
        assert hasattr(self, 'targets')

        # remove outcrossing data
        self.remove_outcrossing()

        # add inbreeding data and reshuffle
        self.add_inbreeding()
        fvecs, targs = shuffle(self.features, self.targets)
        self.set_features(fvecs)
        self.set_targets(targs)


    #################
    ## DEPRECATED? ##
    #################

    ##
    ## REMOVE OUTCROSSING
    ##
    def remove_outcrossing(self):
        if self.out_type == 'multitarget':
            self.multitarget_remove_outcrossing()
        elif self.out_type == 'multiclass':
            self.multiclass_remove_outcrossing()

    def multitarget_remove_outcrossing(self):
        pass

    def multiclass_remove_outcrossing(self):
        fvecs = self.features
        targs = self.targets
        names = self._target_names

        # get outcrossing indices
        out_cols = [i for i, n in enumerate(names) if ', 0.0, ' in n]
        out_rows = [i for i, y in enumerate(targs) if np.argmax(y) in out_cols]

        self.set_features(np.delete(fvecs, out_rows, axis=0))
        self.set_targets(np.delete(targs, out_rows, axis=0))

    ##
    ## RESAMPLING BALANCING SELECTION
    ##
    def resamp_balancing_selection(self):
        '''resample only for training set'''
        if self.out_type == 'multitarget':
            multitarget_resamp_balancing_selection()
        elif self.out_type == 'multiclass':
            multiclass_resamp_balancing_selection()
    
    def multiclass_resamp_balancing_selection(self):
        raise NotImplementedError()
    
    def multitarget_resamp_balancing_selection(self, split_data):
        (x_train, x_test,
         mut_train, mut_test,
         slf_train, slf_test,
         sel_train, sel_test) = split_data

        # get balancing indexes
        bal_idxs, = np.where(sel_train[:,-1] == 1) # NDBa indexes
        n_not_bal = int(np.min(np.sum(sel_train[:,:-1], axis=0)))
        n_bal = len(bal_idxs)
        n_samples = n_not_bal - n_bal
        boot_idxs = resample(bal_idxs, n_samples=n_samples)

        def boot(data):
            boot_data = data[boot_idxs]
            return np.concatenate([data, boot_data], axis=0)

        # resampling
        x_train = boot(x_train)
        mut_train = boot(mut_train)
        slf_train = boot(slf_train)
        sel_train = boot(sel_train)

        # pack data
        split_data = (x_train, x_test,
                      mut_train, mut_test,
                      slf_train, slf_test,
                      sel_train, sel_test
                     )

        return split_data

    ##
    ## REMOVE 2x MUTATION RATE
    ##
    def remove_2x_mutation_rate(self):
        if self.out_type == 'multitarget':
            self.multitarget_remove_2x_mutation_rate()
        elif self.out_type == 'multiclass':
            self.multiclass_remove_2x_mutation_rate()

    def multitarget_remove_2x_mutation_rate(self):
        fvecs = self.features
        targets = self.targets
        self._target_names['mutation'] = self._target_names['mutation'][:-1]
        mutation, selfing, selection = targets.values()
        not_2x = mutation[:,-1]==0
        fvecs = fvecs[not_2x]
        mutation = mutation[not_2x,:-1] # remove 2-1-2 rows/col
        selfing = selfing[not_2x]
        selection = selection[not_2x]
        targets = dict(zip(self.var_keys, [mutation, selfing, selection]))
        
        self.set_features(fvecs)
        self.set_targets(targets)
    
    def multiclass_remove_2x_mutation_rate(self):
        # load multitarget targets for slicing
        fvecs = self.features
        targets = self.targets
        multitargets = self.multitarget_counterpart()
        mutation = self.load_data(multitargets)['mutation']

        not_2x = mutation[:,-1]==0
        self.not_2x = not_2x
        # check for other removes
        if hasattr(self, 'lo_slf'):
            not_2x = np.logical_and(not_2x, self.lo_slf)
        if hasattr(self, 'not_bal'):
            not_2x = np.logical_and(not_2x, self.not_bal)
        fvecs = fvecs[not_2x]
        targets = targets[not_2x]
        mutation, selfing, selection = self.targets.values()
        self.mut_len = mutation.shape[1]
        self.slf_len = selfing.shape[1]
        self.sel_len = selection.shape[1]
        split_data = train_test_split(fvecs,
                                      mutation,
                                      selfing,
                                      selection,
                                      test_size=test_size,
                                      random_state=23)
        if self.sel == 'resamp':
            split_data = self.multitarget_resamp_balancing_selection(split_data)
        # unpack split data
        (self.x_train, self.x_test,
         self.mut_train, self.mut_test,
         self.slf_train, self.slf_test,
         self.sel_train, self.sel_test) = split_data
        self.y_train = [self.mut_train, self.slf_train, self.sel_train]
        self.y_test = [self.mut_test, self.slf_test, self.sel_test]
            
    ##
    ## REMOVE HIGH SELFING RATE
    ##
    def remove_high_selfing_rate(self):
        if self.out_type == 'multitarget':
            self.multitarget_remove_high_selfing_rate()
        elif self.out_type == 'multiclass':
            self.multiclass_remove_high_selfing_rate()

    def multitarget_remove_high_selfing_rate(self):
        pass

    def multiclass_remove_high_selfing_rate(self):
        pass


if __name__ == '__main__':
    features = {'type': 'flat', 'downsample': False}
    targets = {'type': 'multitarget', 'mut': 'keep', 'sel': 'remove'}

    dataset = Dataset(features, targets)
    dataset.prepare_for_training()


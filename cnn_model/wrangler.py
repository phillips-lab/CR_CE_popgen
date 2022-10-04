
""" Wrangles feature vectors and one-hot labels """
import numpy as np
import pickle
import pandas as pd
import glob
from collections import Counter
from os.path import exists



# feature vector .txt files --> wrangled .pkl files
# TODO: empirical data
class Wrangler:
    def __init__(self, sim_dir=None, dest_dir=None):
        # locations of feature vector files
        self.source_dir = '/projects/haldane/shared/worm_landscapes/'
        if sim_dir is None:
            self.source_sim = self.source_dir + '1mb_25subw_beta/SIM/*.txt'
        else:
            self.source_sim = self.source_dir + sim_dir + '*.txt'
        if dest_dir is None:
            self.dest_dir = '/projects/haldane/mlukac/worm-landscapes/data/'
        else:
            self.dest_dir = dest_dir
        
        # empirical features directories
        self.ce_25subw_dir = self.source_dir + '1mb_25subw_beta/CE_sorted/'
        self.cr_25subw_dir = self.source_dir + '1mb_25subw_beta/CR_sorted/'
        self.emp_11subw_dir = self.source_dir + '11wind/'
        self.emp_25subw_dir = self.source_dir + '1mb_25subw_beta/EMP_25sbw/'
        self.ce10 = self.emp_25subw_dir + 'CE10/'
        self.cr10 = self.emp_25subw_dir + 'CR10/'
        self.ce100 = self.emp_25subw_dir + 'CE100/'
        self.cr100 = self.emp_25subw_dir + 'CR100/'
        self.ce10_4ind = self.emp_25subw_dir + 'CE10_4ind/'
        self.cr10_4ind = self.emp_25subw_dir + 'CR10_4ind/'
        self.ce100_4ind = self.emp_25subw_dir + 'CE100_4ind/'
        self.cr100_4ind = self.emp_25subw_dir + 'CR100_4ind/'
        
        # print sim file names
        sim_paths = sorted(glob.glob(self.source_sim))
        sim_names = [sim.split('/')[-1].split('.1mb')[0] for sim in sim_paths]
        sim_param = [sim.split('Ne_5000_')[1] for sim in sim_names]
        print(f'read in {len(sim_names)} simulation names. Format example:')
        print(f'  sim_paths: {sim_paths[0]}')
        print(f'  sim_names: {sim_names[0]}')
        print(f'  sim_param: {sim_param[0]}')
        self.sim_paths = sim_paths
        self.sim_param = sim_param
        
        # feature shape (n_sims, n_stats, n_subw, n_dom)
        self.n_sims = len(sim_paths)
        self.n_stats = 13
        self.n_subw = 25
        self.n_dom = 3
    
    ## INIT FUNCTIONS 
    ## these are called when pulling fresh features
    def init_data_dir(self):
        print('initializing simulation info')
        self.init_sim_info()
        print('assigning target strings')
        self.assign_target_strings()
        print('one hot encoding targets')
        self.init_one_hot_encoding()
        print('initializing feature vectors')
        self.init_fvecs_stack_and_flat()
        print('trimming and downsampling')
        self.downsample_features_and_targets(3)
        
    def init_sim_info(self):
        # sim_info stores the simulation parameter keys/values
        # from this we can explore how balanced the sims are
        sim_param = self.sim_param
        n_sims = self.n_sims
        
        sim_keys = sim_param[0].split('_')[::2] + ['SBaA', 'SBaC']
        print('keys:', sim_keys)
        sim_vals = np.empty((len(sim_keys), n_sims), dtype=object)

        # last two keys are either ['SBA', 'SBC'] or ['SBaA', 'SBaC']
        for i, param in enumerate(sim_param):
            # get sim keys, vals
            keys_i = param.split('_')[::2]  # sim param names
            vals_i = param.split('_')[1::2] # sim param values
            vals_i[0:2] = [str(float(v)) for v in vals_i[0:2]] # Mut and Self
            sim_vals[:-4,i] = vals_i[:-2]

            # deal with selection 
            ben, bal = ['SBA', 'SBC'], ['SBaA', 'SBaC']
            sel_keys = keys_i[-2:]
            sel_vals = vals_i[-2:] 
            # beneficial, not balancing selection
            if sel_keys == ben:
                sim_vals[-4:-2,i] = sel_vals
                sim_vals[-2:,i] = '-1'
            # balancing, not beneficial selection
            elif sel_keys == bal:
                sim_vals[-2:,i] = sel_vals
                sim_vals[-4:-2,i] = '-1'
            else:
                raise Exception(f'selection keys are {sel_keys} but must be one of {ben}, {bal}')
                 
        # print unique sim params
 #       sim_info = dict(zip(sim_keys, sim_vals))
 #       sim_info_ = sim_info.copy()
 #       sim_info_['Mut'] = [str(float(m)) for m in sim_info_['Mut']]
 #               sim_vals[-2:,i] = sel_vals
 #               sim_vals[-4:-2,i] = '-1'
 #           else:
 #               raise Exception(f'selection keys are {sel_keys} but must be one of {ben}, {bal}')
                 
        # print unique sim params
        sim_info = dict(zip(sim_keys, sim_vals))
        sim_info_ = sim_info.copy()
        sim_info_['Mut'] = [str(float(m)) for m in sim_info_['Mut']]
        sim_info_['Self'] = [str(float(s)) for s in sim_info_['Self']]
        print('unique simulation parameters')
        for k, v in sim_info.items():
            print(k, np.unique(v))
            
        self.sim_info = sim_info
        
    def init_fvecs_stack_and_flat(self):        
        fname = self.dest_dir + 'fvecs_stack.pkl'
        if exists(fname):
            print('stacked fvecs already exist but were requested')
        else:
            self.init_fvecs_stack()
            
        fname = self.dest_dir + 'fvecs_flat.pkl'
        if exists(fname):
            print('flat fvecs already exist but were requested')
        else:
            self.init_fvecs_flat()    
            
    # stack (n_sims, n_stats, n_subw/3, 3) for Conv2D
    def init_fvecs_stack(self):
        sim_paths = self.sim_paths
        n_sims = self.n_sims
        n_stats = self.n_stats
        n_subw = self.n_subw
        n_dom = self.n_dom
        
        def stack_sim(sim):
            """ sim is pandas dataframe """
            fvec_stack = np.zeros((n_stats, n_subw, n_dom))
            for i in range(3):
                fvec_stack[:,:,i] = sim.iloc[i,4:].values.reshape(n_stats, n_subw)
            return fvec_stack
        
        fvecs_stack = np.zeros((n_sims, n_stats, n_subw, n_dom))
        for i, path in enumerate(sim_paths):
            sim = pd.read_table(path)
            fvecs_stack[i] = stack_sim(sim)
            n_nans = np.sum(np.isnan(fvecs_stack[i]))
            if n_nans != 0:
                print(f'nans detected in sim {i}')
            
        with open(self.dest_dir + 'fvecs_stack.pkl', 'wb') as f:
            pickle.dump(fvecs_stack, f)
    
    def init_emp_fvecs(self, n_subw):
        sort_key = lambda fname: int(fname.split('_')[-1].split('.')[0])
        worm_chr = {'CE': ['I', 'II', 'III', 'IV'], 'CR': ['I', 'II', 'III', 'V']}
        for worm, chrom in worm_chr.items():
            if n_subw == 11:
                source_dir = self.emp_11subw_dir
            if n_subw == 25 and worm == 'CE':
                source_dir = self.ce_25subw_dir
            if n_subw == 25 and worm == 'CR':
                source_dir = self.cr_25subw_dir
            for c in chrom:
                fname = source_dir + f'{worm}*_chr_{c}_*sorted*txt'
                paths = sorted(glob.glob(fname), key=sort_key)
                self.save_emp_fvecs(n_subw, paths, worm, c)

    def save_emp_fvecs(self, n_subw, paths, worm, chrom):
        n_sims = len(paths)
        n_stats = 13
        n_dom = 3
        dest_dir = self.dest_dir + f'empirical/{n_subw}subw/'
        
        def stack_table(fvec):
            fvec_stack = np.zeros((n_stats, n_subw, n_dom))
            for i in range(n_dom):
                dom_stack = fvec.iloc[i,4:].values.reshape(n_stats, n_subw)
                fvec_stack[:,:,i] = dom_stack
            return fvec_stack

        fvecs_stack = np.zeros((n_sims, n_stats, n_subw, n_dom))
        for i, path in enumerate(paths):
            fvec = pd.read_table(path)
            fvecs_stack[i] = stack_table(fvec)
            n_nans = np.sum(np.isnan(fvecs_stack[i]))
            if n_nans > 0:
                print(f'NaNs detected in {path}')

        with open(dest_dir + f'fvecs_{worm}_chr_{chrom}.pkl', 'wb') as f:
            pickle.dump(fvecs_stack, f)


    ## ONE-HOT ENCODING
    def assign_target_strings(self):
        """use simulation parameters to bin the sims into (mutation, selfing, selection) labels
        mutation and selection don't need to be altered, but selection:
        
          Neutral (N):
          SDA_0_SDC_0_SBA_0_SBC_0

          Neutral & Deleterious (ND):
          SDA_[3|15|7.5]_SDC_[3|15]_SBA_0_SBC_0

          Neutral & Deleterious & Beneficial (NDB):
          SDA_[3|15|7.5]_SDC_[3|15]_SBA_[7.5|15]_SBC_[15]

          Neutral & Deleterious & Balancing (NDBa):
          SDA_[15]_SDC_[15]_SBaA_[15]_SBaC_[15]
          """
        sim_info = self.sim_info
        n_sims = self.n_sims
        
        def classify_selection(i):
            sda, sdc = sim_info['SDA'][i], sim_info['SDC'][i]
            sba, sbc = sim_info['SBA'][i], sim_info['SBC'][i]
            sbaa, sbac = sim_info['SBaA'][i], sim_info['SBaC'][i]

            sel_type = None
            if sda == '0' and sdc == '0':
                sel_type = 'N'
            elif sba == '0' and sbc == '0':
                sel_type = 'ND'
            elif sbaa == '-1' and sbac == '-1':
                sel_type = 'NDB'
            elif sba == '-1' and sbc == '-1':
                sel_type = 'NDBa'
            else:
                raise Exception(f'{sim_names[i]} selection type not classified!')
            return sel_type

        mutation = sim_info['Mut']
        selfing = sim_info['Self']
        selection = [classify_selection(i) for i in range(n_sims)]
        self.target_strings = [mutation, selfing, selection]
        
    def init_one_hot_encoding(self):
        self.multitarget_one_hot_encoding()
        self.multiclass_one_hot_encoding()
        self.multilabel_one_hot_encoding()
        
    def multitarget_one_hot_encoding(self):
        n_sims = self.n_sims
        # multitarget is length mut one-hots, length slf one-hots, length sel one-hots
        # Example:  (1, 0, N) = [(1,0,0), (1,0,0,0,0), (1,0,0,0)]
        names = ['mutation', 'selfing', 'selection']
        values = [['1.0', '1.15', '1.5', '2.0'],          # mutation
                  ['0.0', '0.9', '0.98', '0.999', '1.0'], # selfing
                  ['N', 'ND', 'NDB', 'NDBa']]             # selection
        values_oh = [np.eye(len(v)) for v in values]
        target_strings_oh = [np.zeros((n_sims, len(v))) for v in values]
        
        multitarget_values = dict(zip(names, values))
        multitarget_values_oh = dict(zip(names, values_oh))
        target_strings = dict(zip(names, self.target_strings))
        target_strings_oh = dict(zip(names, target_strings_oh))
        
        # fill in target_strings_oh
        for name, oh in target_strings_oh.items():
            values = multitarget_values[name]
            values_oh = multitarget_values_oh[name]
            strings = target_strings[name]
            for i in range(n_sims):
                j = np.argwhere(np.array(values) == strings[i]).item()
                oh[i] = values_oh[j]
            target_strings_oh[name] = oh
                
        with open(self.dest_dir + 'multitarget_names.pkl', 'wb') as f:
            pickle.dump(multitarget_values, f)
        with open(self.dest_dir + 'multitarget_targets.pkl', 'wb') as f:
            pickle.dump(target_strings_oh, f)
    
    def multiclass_one_hot_encoding(self):
        n_sims = self.n_sims
        names = ['mutation', 'selfing', 'selection']
        values = [['1.0', '1.15', '1.5', '2.0'],          # mutation
                  ['0.0', '0.9', '0.98', '0.999', '1.0'], # selfing
                  ['N', 'ND', 'NDB', 'NDBa']]             # selection
        targets = dict(zip(names, values))
        # multiclass is length mut * slf * sel one-hots
        # Example:  (1, 0, N) = (1,0,...,0)
        multiclass_names = []
        for mut in targets['mutation']:
            for slf in targets['selfing']:
                for sel in targets['selection']:
                    multiclass_names.append(f'({mut}, {slf}, {sel})')
                    
        multiclasses = dict(zip(multiclass_names, np.eye(len(multiclass_names))))
        multiclass_oh = np.zeros((n_sims, len(multiclass_names)))
        
        ts = self.target_strings
        for i in range(n_sims):
            target = f'({ts[0][i]}, {ts[1][i]}, {ts[2][i]})'
            multiclass_oh[i] = multiclasses[target]
        
        with open(self.dest_dir + 'multiclass_targets.pkl', 'wb') as f:
            pickle.dump(multiclass_oh, f)
        with open(self.dest_dir + 'multiclass_names.pkl', 'wb') as f:
            pickle.dump(multiclass_names, f)
        
    def multilabel_one_hot_encoding(self): # (deprecated?)
        # multilabel is length mut+slf+sel three-hots 
        #   (1, 0, N) = (1,0,0,1,0,0,0,0,1,0,0,0)
        with open(self.dest_dir + 'multitarget_targets.pkl', 'rb') as f:
            oh = pickle.load(f)
        mut_oh = oh['mutation']
        slf_oh = oh['selfing']
        sel_oh = oh['selection']
        multilabel_oh = np.hstack((mut_oh, slf_oh, sel_oh))
        
        with open(self.dest_dir + 'multilabel_targets.pkl', 'wb') as f:
            pickle.dump(multilabel_oh, f)

            
    ## DERIVED DATA
    ## these are called when we need to transform primary pickled data
    def multiclass_target_frequency(self):
        with open(self.dest_dir + 'multiclass_targets.pkl', 'rb') as f:
            targets = pickle.load(f)
        with open(self.dest_dir + 'multiclass_names.pkl', 'rb') as f:
            names = pickle.load(f)
        target_strings = [names[np.argmax(y)] for y in targets]
        target_hist = Counter(target_strings)
        for name in names:
            print(name, '\t', target_hist[name])
                        
    ## TRIM/DOWNSAMPLE SUBWINDOWS
    def repeat_targets(self, n_rep):
        # check if already exist
        from os.path import exists
        
        out_type = 'multitarget'
        fname = self.dest_dir + f'trim/{out_type}_targets_rep={n_rep}.pkl'
        if not exists(fname):
            self.repeat_multitarget_targets(n_rep)
        
        out_type = 'multiclass'
        fname = self.dest_dir + f'trim/{out_type}_targets_rep={n_rep}.pkl'
        if not exists(fname):
            self.repeat_multiclass_targets(n_rep)
        
    def repeat_multitarget_targets(self, n_rep):
        with open(self.dest_dir + f'multitarget_targets.pkl', 'rb') as f:
            targets = pickle.load(f)
        for name in ['mutation', 'selfing', 'selection']:
            targets[name] = np.repeat(targets[name], n_rep, axis=0)
        with open(self.dest_dir + f'trim/multitarget_targets_rep={n_rep}.pkl', 'wb') as f:
            pickle.dump(targets, f)
    
    def repeat_multiclass_targets(self, n_rep):
        with open(self.dest_dir + f'multiclass_targets.pkl', 'rb') as f:
            targets = pickle.load(f)
        targets = np.repeat(targets, n_rep, axis=0)
        with open(self.dest_dir + f'trim/multiclass_targets_rep={n_rep}.pkl', 'wb') as f:
            pickle.dump(targets, f)
    
    # Main method for downsampling data
    def downsample_features_and_targets(self, n_rep):
        # targets first
        self.repeat_targets(n_rep)
        
        # now feature vectors
        n_trim = 3
        n_subw = 11
        """Arguments:
            fvecs (np.array): feature vectors fixina be trimmed
            n_trim (int): number of terminal windows that need trimmin
            n_subw (int): number of subwindows per domain to retain
            n_rep (int): number of replicate samples from each feature vec
        """
        # infer desired shape from fvecs
        # (n_sims, n_stats, n_subw, 3 if stack else 1 if flat)
        fv_shape = fvecs.shape
        shape_type = fv_shape[3]
        fvec_type = 'stack' if shape_type==3 else 'flat' if shape_type==1 else ''
        assert len(fvec_type) > 0

        # trimmed fvecs array
        if fvec_type == 'stack':
            n_subw_per_domain = fv_shape[2]
            trim_fv_shape = (n_rep * fv_shape[0], fv_shape[1], n_subw, fv_shape[3])
        elif fvec_type == 'flat':
            assert fv_shape[2] % 3 == 0
            n_subw_per_domain = int(fv_shape[2] / 3)
            trim_fv_shape = (n_rep * fv_shape[0], fv_shape[1], n_subw * 3, fv_shape[3])

        trim_fv = np.zeros(trim_fv_shape)

        rng = np.random.default_rng()
        for n, fv in enumerate(fvecs):
            # separate fvecs to domains left, center, right
            # 1. remove n_trim terminal windows
            ## case a: fvec stack
            if fvec_type == 'stack':
                left_arm, center, right_arm = [fv[:,:,i] for i in range(3)]
                left_trim = left_arm[:, n_trim:]
                right_trim = right_arm[:, :-n_trim]

            ## case b: fvec flat
            elif fvec_type == 'flat':
                left = n_subw_per_domain
                right = 2 * n_subw_per_domain
                left_trim = fv[:, n_trim:left]
                center = fv[:, left:right]
                right_trim = fv[:, right:-n_trim]

            # 2. sample n_subw subwindows
            ## step a: get sample indices for each domain
            idxs_arm = np.arange(n_subw_per_domain - n_trim)
            idxs_center = np.arange(n_subw_per_domain)

            for rep_i in range(n * n_rep, (n+1) * n_rep):
                left_i = sorted(rng.choice(idxs_arm, size=n_subw))
                center_i = sorted(rng.choice(idxs_center, size=n_subw))
                right_i = sorted(rng.choice(idxs_arm, size=n_subw))

                ## step b: slice each domain according to sample
                left = left_trim[:, left_i]
                centr = center[:, center_i]
                right = right_trim[:, right_i]

                # 3. reassemble fvecs
                ## case a: stacked fvecs
                if fvec_type == 'stack':
                    trim_fv[rep_i,:,:,0] = left
                    trim_fv[rep_i,:,:,1] = centr
                    trim_fv[rep_i,:,:,2] = right
                ## case b: flat fvecs
                elif fvec_type == 'flat':
                    trim_fv[rep_i,:, :n_subw] = left
                    trim_fv[rep_i,:, n_subw:2*n_subw] = centr
                    trim_fv[rep_i,:, 2*n_subw:] = right

        # save as trim_fvec_stack/flat
        return trim_fv
    
    ## RESHAPING
    # flat (n_sims, n_stats, n_subw, 1) for Conv2D
    def stack_to_flat(self):
        n_sims = self.n_sims
        n_stats = self.n_stats
        n_subw = self.n_subw
        n_dom = self.n_dom
        fvecs_flat = np.zeros((n_sims, n_stats, n_dom * n_subw, 1))
        
        with open(self.dest_dir + 'fvecs_stack.pkl', 'rb') as f:
            fvecs_stack = pickle.load(f)
        for i, fv in enumerate(fvecs_stack):
            fvecs_flat[i] = np.reshape(fv, (n_stats, 3*n_subw, 1), 'F')
        
        with open(self.dest_dir + 'fvecs_flat.pkl', 'wb') as f:
            pickle.dump(fvecs_flat, f)

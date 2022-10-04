from abc import ABC, abstractmethod
import numpy as np

# plotting
import matplotlib.pyplot as plt
plt.rcParams['figure.figsize'] = (10, 10)
plt.rcParams['figure.dpi'] = 200
import itertools

# tensorflow
from tensorflow.keras.callbacks import *
from tensorflow.keras.models import load_model
from tensorflow.io.gfile import glob, exists, mkdir, remove
from tensorflow.math import confusion_matrix
from tensorflow.random import set_seed

# data wrangling
import pickle
from sklearn.model_selection import train_test_split
from sklearn.utils import resample


class PrintProgress(Callback):
    def __init__(self, max_epochs):
        self.max_epochs = max_epochs
        self.epoch_step = int(max_epochs * 0.1)
        self.percent = 0.
        super().__init__()

    def on_epoch_end(self, epoch, logs=None):
        # print epoch and logs every 5 or 10 percent
        if epoch % self.epoch_step == 0:
            print(f'epoch {epoch}/{self.max_epochs}')


##################
#   BASE MODEL   #
##################
class BaseClassifier(ABC):
    def __init__(self, 
                 model_name,  # for model_log dir 
                 dataset,     # Dataset object
                 vnum=None):  # optional for loading pretrained
        self.dataset = dataset
        self.in_type = dataset.in_type
        self.out_type = dataset.out_type
        self.model = None
        self.is_trained = False
        self.set_logpath(model_name, vnum)

    #################
    #   MODEL FIT   #
    #################
    @abstractmethod
    def compile_model(self):
        """ construct network layers """
        raise NotImplementedError()
    
    def fit_model(self, max_epochs, verbose=0, summarize=False, top_k=5):
        """ fit to training data """
        self.dataset.prepare_for_training()
        if self.model is None:
            self.model = self.compile_model()
        else:
            self.model = self.load_best_model()

        if summarize: self.model.summary()

        # callbacks
        callbacks = self.get_callbacks(max_epochs)

        # train network
        x_train = self.dataset.x_train
        y_train = self.dataset.y_train
        valid_data = (self.dataset.x_test, self.dataset.y_test)
        #set_seed(23)
        self.history = self.model.fit(x_train, y_train, 
                                      validation_data=valid_data, 
                                      epochs=max_epochs, 
                                      callbacks=callbacks,
                                      verbose=verbose)
        self.is_trained = True

        # clean checkpoints
        self.keep_top_k_models(k=top_k)

        # pickle training history
        with open(self.logpath + 'training_history.pkl', 'wb') as f:
            pickle.dump(self.history.history, f)
    
    # TODO handle multilabel case maybe?
    def keep_top_k_models(self, k, metric='acc'):
        if metric == 'acc':
            self.keep_top_k_models_from_acc(k)
        elif metric == 'loss':
            self.keep_top_k_models_from_loss(k)

    def keep_top_k_models_from_acc(self, k):
        # get checkpoint filenames and accuracies
        ckpt_filepaths = glob(self.logpath + '*hdf5')
        accs = self.get_ckpt_accs()

        # get best and worst model indexes
        ordered_idxs = np.argsort(accs)
        best_idxs = ordered_idxs[-k:]
        worst_idxs = ordered_idxs[:-k]

        # save best path, remove not top k checkpoints
        self.best_model_path = ckpt_filepaths[best_idxs[-1]]
        for i in worst_idxs:
            remove(ckpt_filepaths[i])

    def keep_top_k_models_from_loss(self, k):
        # get checkpoint filenames and accuracies
        ckpt_filepaths = glob(self.logpath + '*hdf5')
        losses = self.get_ckpt_losses()

        # get best and worst model indexes
        ordered_idxs = np.argsort(losses)
        best_idxs = ordered_idxs[:k]
        worst_idxs = ordered_idxs[k:]

        # save best path, remove not top k checkpoints
        self.best_model_path = ckpt_filepaths[best_idxs[0]]
        for i in worst_idxs:
            remove(ckpt_filepaths[i])

    def load_best_model(self):
        if not hasattr(self, 'best_model_path'):
            self.set_best_model_path()
        print(f'loading best model {self.best_model_path.split("/")[-1]}')
        return load_model(self.best_model_path)

    def load_data(self, name):
        filename = f'data/{name}.pkl'
        with open(filename, 'rb') as f:
            data = pickle.load(f)
        return data

    def predict(self, filename):
        """ Predicts on empirical data pickled in filename """
        assert self.is_trained 

        # load empirical data
        data = self.load_data(filename)
        preds = self.model.predict(data)
        return preds

    #####################
    #   GETS AND SETS   #
    #####################
    def set_logpath(self, model_name, vnum=None):
        if vnum is None:
            self.set_logpath_(model_name)
        else:
            logpath = f'model_logs/{model_name}/version_{vnum}/'
            assert exists(logpath)
            self.vnum = vnum
            self.logpath = logpath

    def set_logpath_(self, model_name):
        """ Sets model version number and log path """
        # make model name directory if not exist
        if not exists(f'model_logs/{model_name}'):
            mkdir(f'model_logs/{model_name}')

        # make version directory
        if not exists(f'model_logs/{model_name}/version_0'):
            vnum = 0
        else:
            # vnum is 1+ most recent vnum
            versions = glob(f'model_logs/{model_name}/version_*')
            current_vnum = max([int(v.split('_')[-1]) for v in versions])
            vnum = current_vnum + 1

        # save version number and logpath
        self.vnum = vnum
        self.logpath = f'model_logs/{model_name}/version_{vnum}/'
        mkdir(self.logpath)

    def set_best_model_path(self):
        # get checkpoint accuracy metrics
        ckpt_filepaths = glob(self.logpath + '*hdf5')
        accs = self.get_ckpt_accs()
        self.best_model_path = ckpt_filepaths[np.argmax(accs)]
        
    def get_ckpt_losses(self):
        ckpt_filepaths = glob(self.logpath + '*hdf5')
        ckpts = [c.split('/')[-1] for c in ckpt_filepaths]
        
        if self.out_type == 'multitarget':
            # get checkpoint selection accuracy
            accs = [c.split('-')[2] for c in ckpts]
            accs = [float(s.split('=')[-1]) for s in accs]
        elif self.out_type == 'multiclass':
            loss = [c.split('=')[2] for c in ckpts]
            loss = [float(c.split("-")[0]) for c in loss]
        return loss

    def get_ckpt_accs(self):
        ckpt_filepaths = glob(self.logpath + '*hdf5')
        ckpts = [c.split('/')[-1] for c in ckpt_filepaths]
        
        if self.out_type == 'multitarget':
            # get checkpoint selection accuracy
            accs = [c.split('-')[2] for c in ckpts]
            accs = [float(s.split('=')[-1]) for s in accs]
        elif self.out_type == 'multiclass':
            accs = [c.split('=')[2] for c in ckpts]
            accs = [float(f'{c.split(".")[0]}.{c.split(".")[1]}') for c in accs]
        return accs

    # TODO handle multilabel maybe?
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
            fname += '-val_acc={val_accuracy:.3f}.hdf5'
            monitor = 'accuracy'
        ckpt = ModelCheckpoint(filepath=fname, 
                               monitor=monitor,
                               save_best_only=True,
                               mode='max')
        prog = PrintProgress(max_epochs)
        stop = EarlyStopping('val_loss', min_delta=1e-7, patience=100)

        return [ckpt, prog, stop]

    def get_history(self, *keys):
        """ Unpacks keys from history """
        assert self.model is not None and self.history is not None
        history = self.history.history
        key_histories = []
        for key in keys:
            key_histories.append(history[key])
        return tuple(key_histories)

    def get_cm(self, truth, preds, softmax=True, num_classes=None):
        if softmax:
            truth = np.argmax(truth, axis=1)
            preds = np.argmax(preds, axis=1)
        cm = confusion_matrix(truth, preds, num_classes=num_classes)
        return cm.numpy()

    def get_int_from_softmax(self, softmax_dict):
        """
        Given softmax dict with keys=(mutation, selfing, selection)
        returns the string representation of each prediction
        """
        # convert softmax predictions to argmax then strings
        int_dict = softmax_dict.copy()
        for var, str_labs in self.dataset.target_names.items():
            int_dict[var] = np.argmax(softmax_dict[var], axis=1)
        return int_dict
        
    def get_string_from_softmax(self, softmax_dict):
        """
        Given softmax dict with keys=(mutation, selfing, selection)
        returns the string representation of each prediction
        """
        # convert softmax predictions to argmax then strings
        str_dict = softmax_dict.copy()
        
        # multitarget case
        if self.out_type == 'multitarget':
            for var, str_labs in self.dataset.target_names.items():
                argmax = np.argmax(softmax_dict[var], axis=1)
                str_dict[var] = np.array([str_labs[i] for i in argmax])
        else:
            raise Exception('this function only used for multitarget')
        
        return str_dict

    def get_prediction_dict(self):
        """ 
        Predicts on test data
        Returns dictionary with var strings as keys
        and string label predictions as values
        """
        assert self.is_trained
        model = self.load_best_model()

        if self.out_type == 'multitarget':
            # predict with model, convert to dict, get string labels
            softmax_preds = model.predict(self.dataset.x_test)
            softmax_preds_dict = dict(zip(self.var_keys, softmax_preds))
            preds = self.get_string_from_softmax(softmax_preds_dict)

            # get ground truth string labels
            softmax_truth = [self.mut_test, self.slf_test, self.sel_test]
            softmax_truth_dict = dict(zip(self.var_keys, softmax_truth))
            truth = self.get_string_from_softmax(softmax_truth_dict)

            # merge to single dictionary
            prediction_results = dict()
            for var in self.var_keys:
                var_truth = f'{var}_truth'
                var_pred = f'{var}_pred'
                prediction_results[var_truth] = truth[var]
                prediction_results[var_pred] = preds[var]

        elif self.out_type == 'multiclass':
            # predict with model, get string predictions
            softmax_preds = model.predict(self.dataset.x_test)
            argmax_preds = np.argmax(softmax_preds, axis=1)
            preds = [self.dataset.target_names[a] for a in argmax_preds]

            # get ground truth string predictions
            softmax_truth = self.dataset.y_test
            argmax_truth = np.argmax(softmax_truth, axis=1)
            truth = [self.dataset.target_names[a] for a in argmax_truth]

            # merge prediction results into single dictionary
            prediction_results = dict(zip(['truth', 'preds'], [truth, preds]))

        return prediction_results

    ################
    #   PLOTTING   #
    ################
    def str2ints(self, varname, str_truth, str_preds):
        int_truth = str_truth.copy()
        int_preds = str_preds.copy()
        for idx, vlab in enumerate(self.dataset.target_names[varname]):
            ii, = np.where(str_truth==vlab)
            jj, = np.where(str_preds==vlab)
            int_truth[ii] = idx
            int_preds[jj] = idx
        int_truth = int_truth.astype('int')
        int_preds = int_preds.astype('int')
        return int_truth, int_preds

    def get_cmcm(self):
        """ 
        Plots mutation, selfing confusion matrices sliced by entries in
        selection confusion matrix. A confusion matrix of confusion matrices.
        """
        # get prediction dictionary and selection cm
        pred_results = self.get_prediction_dict()
        sel_truth = pred_results['selection_truth']
        sel_preds = pred_results['selection_pred']
        sel_ints = self.str2ints('selection', sel_truth, sel_preds)
        sel_int_truth, sel_int_preds = sel_ints
        sel_cm = self.get_cm(sel_int_truth, sel_int_preds, 
                             softmax=False, num_classes=self.sel_len)

        # get mutation and selfing sliced confusion matrices
        sel_len = self.sel_len
        mut_len = self.mut_len
        slf_len = self.slf_len
        mut_cms = np.zeros((sel_len, sel_len, mut_len, mut_len), dtype='int')
        slf_cms = np.zeros((sel_len, sel_len, slf_len, slf_len), dtype='int')

        for i in range(sel_len):
            for j in range(sel_len):
                # slice into selection cm pairs
                true_lab = self.dataset.target_names['selection'][i]
                pred_lab = self.dataset.target_names['selection'][j]
                truth = pred_results['selection_truth']
                preds = pred_results['selection_pred']
                i_slice = truth==true_lab
                j_slice = preds==pred_lab
                ij_slice = np.logical_and(i_slice, j_slice)

                assert sum(ij_slice) == sel_cm[i,j] #sanity check

                # inner variable cms
                def inner_cm(varname):
                    n_classes = len(set(pred_results[f'{varname}_truth']))
                    ij_truth = pred_results[f'{varname}_truth'][ij_slice]
                    ij_pred = pred_results[f'{varname}_pred'][ij_slice]
                    ij_truth, ij_pred = self.str2ints(varname, ij_truth, ij_pred)
                    ij_cm = self.get_cm(ij_truth, ij_pred, 
                                softmax=False, num_classes=n_classes)
                    return ij_cm
                mut_cms[i,j] = inner_cm('mutation')
                slf_cms[i,j] = inner_cm('selfing')
                    
        return mut_cms, slf_cms

    def _plot_cmcm(self, cms, varname):
        out_len, _, in_len, _ = cms.shape
        cmap = plt.cm.Blues
        fig, ax = plt.subplots(out_len, out_len, dpi=200, 
                               sharex=True, sharey=True)
        plt.suptitle(f'selection/{varname} heirarchical confusion matrix',
                    fontsize=20)
        for i in range(out_len):
            for j in range(out_len):
                ij_cm = cms[i,j]
                ax[i,j].imshow(ij_cm, interpolation='nearest', cmap=cmap)
                tick_marks = np.arange(in_len)
                if j==0:
                    ax[i,j].set_ylabel(self.dataset.target_names['selection'][i], fontsize=15)
                if i==out_len-1:
                    ax[i,j].set_xlabel(self.dataset.target_names['selection'][j], fontsize=15)
                ax[i,j].set_xticks(tick_marks)
                ax[i,j].set_xticklabels(self.dataset.target_names[varname])
                ax[i,j].set_yticks(tick_marks)
                ax[i,j].set_yticklabels(self.dataset.target_names[varname])
                thresh = ij_cm.max() / 2.
                for r, c in itertools.product(range(in_len), range(in_len)):
                    clr = 'white' if ij_cm[r,c] > thresh else 'black'
                    ax[i,j].text(c, r, format(ij_cm[r,c], 'd'),
                                 horizontalalignment='center',
                                 color=clr)
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        plt.savefig(self.logpath + f'{varname}_cmcm.png')
        plt.close()

    def plot_cmcm(self):
        mut_cms, slf_cms = self.get_cmcm()
        self._plot_cmcm(mut_cms, 'mutation')
        self._plot_cmcm(slf_cms, 'selfing')

    def plot_cm(self, cm, classes, filename, 
                normalize=False, 
                cmap=plt.cm.Blues):
        cm_norm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
        print(f'{filename} confusion matrix')
        print(cm_norm)
        if normalize: cm = cm_norm

        plt.imshow(cm, interpolation='nearest', cmap=cmap)
        if self.out_type == 'multitarget':
            plt.title('%s confusion matrix' % filename, fontsize=15)
            plt.colorbar()
            tick_marks = np.arange(len(classes))
            plt.xticks(tick_marks, classes, rotation=45)
            plt.yticks(tick_marks, classes)
        elif self.out_type == 'multiclass':
            plt.title('confusion matrix', fontsize=15)
            plt.colorbar()
            tick_marks = np.arange(len(classes))
            plt.xticks(tick_marks, classes, rotation=90, fontsize=7)
            plt.yticks(tick_marks, classes, fontsize=7)

        if self.out_type == 'multitarget':
            fmt = '.2f' if normalize else 'd'
            thresh = cm.max() / 2.
            rows, cols = cm.shape
            for i, j in itertools.product(range(rows), range(cols)):
                plt.text(j, i, format(cm[i, j], fmt),
                         horizontalalignment='center',
                         color='white' if cm[i, j] > thresh else 'black')

        plt.ylabel('True label', fontsize=12)
        plt.xlabel('Predicted label', fontsize=12)
        plt.tight_layout()
        plt.savefig(self.logpath + f'{filename}.png')
        plt.close()

    def plot_confusion_matrix(self):
        assert self.is_trained
        # load best checkpoint
        self.model = self.load_best_model()

        preds = self.model.predict(self.dataset.x_test)
        names = self.load_data(f'{self.out_type}_names')

        if self.out_type == 'multitarget':
            # plot mutation confusion matrix
            mut_cm = self.get_cm(self.mut_test, preds[0])
            mut_labs = names['mutation']
            self.plot_cm(mut_cm, mut_labs, 'mutation')

            # plot selfing confusion matrix
            slf_cm = self.get_cm(self.slf_test, preds[1])
            slf_labs = names['selfing']
            self.plot_cm(slf_cm, slf_labs, 'selfing')

            # plot selection confusion matrix
            sel_cm = self.get_cm(self.sel_test, preds[2])
            sel_labs = names['selection']
            self.plot_cm(sel_cm, sel_labs, 'selection')
        elif self.out_type == 'multiclass':
            cm = self.get_cm(self.dataset.y_test, preds)
            self.plot_cm(cm, self.dataset.target_names, 'confusion')

    def plot_history(self, key):
        """ plot training and validation history unless learning rate """
        assert key in self.history.history.keys()
        fig, ax = plt.subplots(dpi=200)
        if key == 'lr':
            lr, = self.get_history(key)
            ax.plot(lr, label='learning rate')
            plt.xlabel('epoch')
            plt.ylabel('learning rate')
        else:
            train_key, valid_key = key, f'val_{key}'
            training, validation = self.get_history(train_key, valid_key)
            ax.plot(training, label='training')
            ax.plot(validation, label='validation')
            plt.xlabel('epoch')
            plt.ylabel(key)
            plt.legend()

        # save in logpath
        plt.savefig(self.logpath + f'{key}.png')
        plt.close()

    def plot_loss(self):
        """ Plot and save training losses """
        if self.out_type == 'multitarget':
            for k in ['loss', 'mutation_loss', 'selfing_loss', 'selection_loss']:
                self.plot_history(k)
        else:
            self.plot_history('loss')

    def plot_acc(self):
        """ Plot and save training accuracies """
        if self.out_type == 'multitarget':
            for k in ['mutation_accuracy', 'selfing_accuracy', 'selection_accuracy']:
                self.plot_history(k)
        else:
            self.plot_history('accuracy')

    def plot_lr(self):
        """ Plot and save learning rate vs epochs """
        self.plot_history('lr')
    
    def plot_var_counts(self, cond_on='selection', plot_truth=False):
        """
        Visualization of 4D lattice where 
        for each (cond_on_truth, cond_on_pred) pair of classes
        there is a cross section matrix with var1_classes rows
        and var2_classes cols.
        The color represents relative frequencies within each cross section
        The numbers represent frequency of that (var1, var2) pair
        among cond_on_truth
        """
        # get condition on truth, predictions, misclassified boolean, labels
        pred_results = self.get_prediction_dict()
        cond_on_truth = pred_results[f'{cond_on}_truth'] # obs labels
        cond_on_pred = pred_results[f'{cond_on}_pred'] # obs pred labels
        cond_on_misclass = cond_on_truth!=cond_on_pred # boolean
        cond_on_labs = self.dataset.target_names[cond_on] # unique labels
        n_labs = len(cond_on_labs) # number of unique labels
        
        # get not cond on variable names and unique labels
        im_vars = [x for x in self.var_keys if x!=cond_on]
        im_labs = [self.dataset.target_names[var] for var in im_vars]
            
        # plot data
        image_rows, image_cols = [len(im_lab) for im_lab in im_labs]
        tensor_shape = (n_labs, n_labs, image_rows, image_cols)
        image_tensor_counts = np.zeros(tensor_shape)
        image_tensor_freqs = np.zeros_like(image_tensor_counts)
        
        # plot 2d hists
        fig, ax = plt.subplots(n_labs, n_labs, 
                                dpi=200, 
                                sharex=True, sharey=True)
        for i in range(n_labs): # truth index
            for j in range(n_labs): # prediction index
                # row and column labels and counts 
                true_class = cond_on_labs[i]
                pred_class = cond_on_labs[j]
                if i==n_labs-1: ax[i,j].set_xlabel(pred_class, fontsize=12)
                if j==0: ax[i,j].set_ylabel(true_class, fontsize=12)
                row_bool = cond_on_truth==true_class
                col_bool = cond_on_pred==pred_class
                plot_bool = np.logical_and(row_bool, col_bool)
                
                # set ticks
                ax[i,j].set_xticks(range(len(im_labs[1])))
                ax[i,j].set_xticklabels(im_labs[1])
                ax[i,j].set_yticks(range(len(im_labs[0])))
                ax[i,j].set_yticklabels(im_labs[0])
                
                # make image of count data
                for ii in range(image_rows): # var1 index
                    for jj in range(image_cols): # var2 index
                        if plot_truth:
                            var1 = pred_results[f'{im_vars[0]}_truth']
                            var2 = pred_results[f'{im_vars[1]}_truth']
                            fig_title = f'{im_vars[0]} and {im_vars[1]} joint ground truth'
                        else:
                            var1 = pred_results[f'{im_vars[0]}_pred']
                            var2 = pred_results[f'{im_vars[1]}_pred']
                            fig_title = f'{im_vars[0]} and {im_vars[1]} joint predictions'
                        var1_bool = var1==im_labs[0][ii]
                        var2_bool = var2==im_labs[1][jj]
                        var1_and_var2 = np.logical_and(var1_bool, var2_bool)
                        marginal_class = np.logical_and(plot_bool, var1_and_var2)
                        class_count = int(sum(marginal_class))
                        image_tensor_counts[i,j,ii,jj] = class_count
                
        # make frequency data
        sums = np.sum(image_tensor_counts, axis=1, keepdims=True) # keepdims for broadcasting
        image_tensor_freqs = np.divide(image_tensor_counts, sums, 
                                       out=image_tensor_freqs, 
                                       where=sums!=0)
        for i in range(n_labs):
            for j in range(n_labs):
                for ii in range(image_rows):
                    for jj in range(image_cols):
                        number = round(image_tensor_freqs[i,j,ii,jj], 2)
                        thresh = image_tensor_counts[i,j,:,:].max() / 2
                        color = 'white' if int(thresh)==0 else 'white'
                        color = 'white' if image_tensor_counts[i,j,ii,jj]<=thresh else 'black'
                        ax[i,j].text(jj, ii, number, 
                                     ha='center', 
                                     va='center', 
                                     color=color)
                ax[i,j].imshow(image_tensor_counts[i,j,:,:])
        fig.tight_layout()
        
        # padding dicts
        title = {'mutation': 1., 'selfing': 0.9, 'selection': 0.82}
        xtitle = {'mutation': -0.02, 'selfing': 0.09, 'selection': 0.13}
        ytitle = {'mutation': 0.02, 'selfing': -0.02, 'selection': -0.03}
        hspace = {'mutation': 0.05, 'selfing': -0.6, 'selection': -0.7}
        wspace = {'mutation': -0.4, 'selfing': 0.1, 'selection': 0.1}
        
        # axis titles
        fig.text(0.5, title[cond_on], fig_title, 
                    ha='center', 
                    fontsize=20)
        fig.text(0.5, xtitle[cond_on], f'predicted {cond_on}', 
                    ha='center', 
                    fontsize=15)
        fig.text(ytitle[cond_on], 0.5, f'true {cond_on}', 
                    va='center', 
                    rotation='vertical', 
                    fontsize=15)
        
        plt.subplots_adjust(hspace=hspace[cond_on], wspace=wspace[cond_on])
        plt.savefig(self.logpath + f'{cond_on}_counts.png')
        plt.close()


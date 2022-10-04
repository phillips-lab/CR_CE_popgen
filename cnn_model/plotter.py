import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['figure.figsize'] = (10, 10)
plt.rcParams['figure.dpi'] = 200
import itertools
import pickle



class Plotter:
    def __init__(self, dataset):
        self.dataset = dataset


class ModelPlotter(Plotter):
    def __init__(self, model):
        self.super().__init__(model.dataset)
        self.model = model
        
    def get_cm(self, truth, preds, softmax=True, num_classes=None):
        if softmax:
            truth = np.argmax(truth, axis=1)
            preds = np.argmax(preds, axis=1)
        cm = confusion_matrix(truth, preds, num_classes=num_classes)
        return cm.numpy()

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
                    ylabs = self.dataset.target_names['selection'][i]
                    ax[i,j].set_ylabel(ylabs, fontsize=15)
                if i==out_len-1:
                    xlabs = self.dataset.target_names['selection'][j]
                    ax[i,j].set_xlabel(xlabs, fontsize=15)
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
                            fig_title = f'{im_vars[0]} and {im_vars[1]} joint truth'
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
        sums = np.sum(image_tensor_counts, axis=1, keepdims=True)
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
                        if image_tensor_counts[i,j,ii,jj]<=thresh:
                            color = 'white'  
                        else: 
                            color = 'black'
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

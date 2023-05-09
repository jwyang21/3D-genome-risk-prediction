# used the codes of https://github.com/jaredleekatzman/DeepSurv/blob/41eed003e5b892c81e7855e400861fa7a2d9da4f/deepsurv/deep_surv.py, after proper modification
# referred to https://github.com/JoonHyeongPark/ExplainableMLKGNN/blob/main/model.py

from __future__ import print_function, absolute_import

import lasagne
import numpy
import time
import json
import h5py

import theano
import theano.tensor as T

from lifelines.utils import concordance_index

from deepsurv_logger import DeepSurvLogger

from lasagne.regularization import regularize_layer_params, l1, l2
from lasagne.nonlinearities import rectify,selu

import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.rcParams['figure.dpi'] = 300 #150
plt.rc('font', family='FreeSans', size=7)
plt.rc('figure', figsize=(1.5, 1.5))

plt.rc('xtick', labelsize=7)
plt.rc('ytick', labelsize=7)

class DeepSurv:
    def __init__(self, n_in,
    learning_rate, hidden_layers_sizes = None,
    lr_decay = 0.0, momentum = 0.9,
    L2_reg = 0.0, L1_reg = 0.0,
    activation = "rectify",
    dropout = None,
    batch_norm = False,
    standardize = False,
    ):
        self.X = T.fmatrix('x')  
        self.E = T.ivector('e') 

        self.offset = numpy.zeros(shape = n_in, dtype=numpy.float32)
        self.scale = numpy.ones(shape = n_in, dtype=numpy.float32)

        network = lasagne.layers.InputLayer(shape=(None,n_in),
            input_var = self.X)

        self.standardize = standardize

        if activation == 'rectify':
            activation_fn = rectify
        elif activation == 'selu':
            activation_fn = selu
        else:
            raise IllegalArgumentException("Unknown activation function: %s" % activation)

        for n_layer in (hidden_layers_sizes or []):
            if activation_fn == lasagne.nonlinearities.rectify:
                W_init = lasagne.init.GlorotUniform()
            else:
                W_init = lasagne.init.GlorotUniform()

            network = lasagne.layers.DenseLayer(
                network, num_units = n_layer,
                nonlinearity = activation_fn,
                W = W_init
            )

            if batch_norm:
                network = lasagne.layers.batch_norm(network)

            if not dropout is None:
                network = lasagne.layers.DropoutLayer(network, p = dropout)

        network = lasagne.layers.DenseLayer(
            network, num_units = 1,
            nonlinearity = lasagne.nonlinearities.linear,
            W = lasagne.init.GlorotUniform()
        )

        self.network = network
        self.params = lasagne.layers.get_all_params(self.network,
                                                    trainable = True)
        self.hidden_layers = lasagne.layers.get_all_layers(self.network)[1:]

        self.partial_hazard = T.exp(self.risk(deterministic = True)) 

        self.hyperparams = {
            'n_in': n_in,
            'learning_rate': learning_rate,
            'hidden_layers_sizes': hidden_layers_sizes,
            'lr_decay': lr_decay,
            'momentum': momentum,
            'L2_reg': L2_reg,
            'L1_reg': L1_reg,
            'activation': activation,
            'dropout': dropout,
            'batch_norm': batch_norm,
            'standardize': standardize
        }

        self.n_in = n_in
        self.learning_rate = learning_rate
        self.lr_decay = lr_decay
        self.L2_reg = L2_reg
        self.L1_reg = L1_reg
        self.momentum = momentum
        self.restored_update_params = None

    def _negative_log_likelihood(self, E, deterministic = False):
        risk = self.risk(deterministic)
        hazard_ratio = T.exp(risk)
        log_risk = T.log(T.extra_ops.cumsum(hazard_ratio))
        uncensored_likelihood = risk.T - log_risk
        censored_likelihood = uncensored_likelihood * E
        num_observed_events = T.sum(E)
        neg_likelihood = -T.sum(censored_likelihood) / num_observed_events
        return neg_likelihood

    def _get_loss_updates(self,
    L1_reg = 0.0, L2_reg = 0.001,
    update_fn = lasagne.updates.nesterov_momentum,
    max_norm = None, deterministic = False,
    momentum = 0.9,
    **kwargs):
        loss = (
            self._negative_log_likelihood(self.E, deterministic)
            + regularize_layer_params(self.network,l1) * L1_reg
            + regularize_layer_params(self.network, l2) * L2_reg
        )

        if max_norm:
            grads = T.grad(loss,self.params)
            scaled_grads = lasagne.updates.total_norm_constraint(grads, max_norm)
            updates = update_fn(
                scaled_grads, self.params, **kwargs
            )
        else:
            updates = update_fn(
                loss, self.params, **kwargs
            )

        if momentum:
            updates = lasagne.updates.apply_nesterov_momentum(updates, 
                self.params, self.learning_rate, momentum=momentum)

        if self.restored_update_params:
            for p, value in zip(updates.keys(), self.restored_update_params):
                p.set_value(value)
            self.restored_update_params = None
            
        self.updates = updates

        return loss, updates

    def _get_train_valid_fn(self,
    L1_reg, L2_reg, learning_rate,
    **kwargs):

        loss, updates = self._get_loss_updates(
            L1_reg, L2_reg, deterministic = False,
            learning_rate=learning_rate, **kwargs
        )
        train_fn = theano.function(
            inputs = [self.X, self.E],
            outputs = loss,
            updates = updates,
            name = 'train'
        )

        valid_loss, _ = self._get_loss_updates(
            L1_reg, L2_reg, deterministic = True,
            learning_rate=learning_rate, **kwargs
        )

        valid_fn = theano.function(
            inputs = [self.X, self.E],
            outputs = valid_loss,
            name = 'valid'
        )
        return train_fn, valid_fn

    def get_concordance_index(self, x, t, e, **kwargs):
        compute_hazards = theano.function(
            inputs = [self.X],
            outputs = -self.partial_hazard
        )
        partial_hazards = compute_hazards(x)
        ### EDIT
        if numpy.isnan(partial_hazards).sum().sum() == partial_hazards.shape[0]:
            partial_hazards2 = partial_hazards.copy()
            identical_value = 0 #set partial_hazards of all samples as a constant, since all partial_hazard values are NaN which makes computing c-index impossible
            for i in range(partial_hazards.shape[0]):
                partial_hazards2[i] = [identical_value]
            return concordance_index(t, partial_hazards2, e)
        else:
            return concordance_index(t, partial_hazards, e)

    def _standardize_x(self, x):
        return (x - self.offset) / self.scale

    def prepare_data(self,dataset):
        if isinstance(dataset, dict):
            x, e, t = dataset['x'], dataset['e'], dataset['t']

        if self.standardize:
            x = self._standardize_x(x)

        sort_idx = numpy.argsort(t)[::-1]
        x = x[sort_idx]
        e = e[sort_idx]
        t = t[sort_idx]

        return (x, e, t)

    def train(self,
    train_data, valid_data= None,
    n_epochs = 500,
    logger = None,
    update_fn = lasagne.updates.nesterov_momentum,
    verbose = True,
    **kwargs):
        if logger is None:
            logger = DeepSurvLogger('DeepSurv')

        if self.standardize:
            self.offset = train_data['x'].mean(axis = 0)
            self.scale = train_data['x'].std(axis = 0)

        x_train, e_train, t_train = self.prepare_data(train_data)

        if valid_data:
            x_valid, e_valid, t_valid = self.prepare_data(valid_data)

        best_validation_loss = numpy.inf
        best_params = None
        best_params_idx = -1

        lr = theano.shared(numpy.array(self.learning_rate,
                                    dtype = numpy.float32))
        momentum = numpy.array(0, dtype= numpy.float32)

        train_fn, valid_fn = self._get_train_valid_fn(
            L1_reg=self.L1_reg, L2_reg=self.L2_reg,
            learning_rate=lr,
            momentum = momentum,
            update_fn = update_fn, **kwargs
        )

        start = time.time()
        early_stopping_flag = 0
        early_stopping_num_epochs = 10 #acts as PATIENCE which is typically used in early stopping
        min_train_epoch = 25
        
        for epoch in range(n_epochs):
            # Power-Learning Rate Decay
            lr = self.learning_rate / (1 + epoch * self.lr_decay)
            logger.logValue('lr', lr, epoch)

            if self.momentum and epoch >= 10:
                momentum = self.momentum

            loss = train_fn(x_train, e_train)
            
            #### EDIT
            # if loss is NaN, replace loss with a very large number to prevent error
            if str(loss) == 'nan':
                print("replace NaN loss with a very large value (99999)")
                loss = 99999

            logger.logValue('loss', loss, epoch)

            ci_train = self.get_concordance_index(
                x_train,
                t_train,
                e_train,
            )
            logger.logValue('c-index',ci_train, epoch)
            
            if verbose and (epoch % 100 == 0):
                logger.print_progress_bar(epoch, n_epochs, loss, ci_train)
                
            # validation
            if valid_data: #record validation performance for every epoch. 
                validation_loss = valid_fn(x_valid, e_valid)
                logger.logValue('valid_loss', validation_loss, epoch)

                ci_valid = self.get_concordance_index(
                    x_valid,
                    t_valid,
                    e_valid
                )
                logger.logValue('valid_c-index', ci_valid, epoch)

                if validation_loss < best_validation_loss:
                    best_params = [param.copy().eval() for param in self.params]
                    best_params_idx = epoch
                    best_validation_loss = validation_loss
                    early_stopping_flag = 0
                else:
                    if epoch > min_train_epoch:
                        early_stopping_flag += 1

            if early_stopping_flag >= early_stopping_num_epochs:
                if epoch > min_train_epoch:
                    break

        if verbose:
            logger.logMessage('Finished Training with %d iterations in %0.2fs' % (
                epoch + 1, time.time() - start
            ))
        logger.shutdown()
        
        logger.history['best_valid_loss'] = best_validation_loss
        logger.history['best_params'] = best_params
        logger.history['best_params_idx'] = best_params_idx

        return logger.history

    def to_json(self):
        return json.dumps(self.hyperparams)

    def save_model(self, filename, weights_file = None):
        with open(filename, 'w') as fp:
            fp.write(self.to_json())

        if weights_file:
            self.save_weights(weights_file)

    def save_weights(self,filename):
        def save_list_by_idx(group, lst):
            for (idx, param) in enumerate(lst):
                group.create_dataset(str(idx), data=param)

        weights_out = lasagne.layers.get_all_param_values(self.network, trainable=False)
        if self.updates:
            updates_out = [p.get_value() for p in self.updates.keys()]
        else:
            raise Exception("Model has not been trained: no params to save!")

        with h5py.File(filename, 'w') as f_out:
            weights_grp = f_out.create_group('weights')
            save_list_by_idx(weights_grp, weights_out)

            updates_grp = f_out.create_group('updates')
            save_list_by_idx(updates_grp, updates_out)

    def load_weights(self, filename):
        def load_all_keys(fp):
            results = []
            for key in fp:
                dataset = fp[key][:]
                results.append((int(key), dataset))
            return results

        def sort_params_by_idx(params):
            return [param for (idx, param) in sorted(params, 
            key=lambda param: param[0])]

        with h5py.File(filename, 'r') as f_in:
            weights_in = load_all_keys(f_in['weights'])
            updates_in = load_all_keys(f_in['updates'])

        sorted_weights_in = sort_params_by_idx(weights_in)
        lasagne.layers.set_all_param_values(self.network, sorted_weights_in, 
            trainable=False)

        sorted_updates_in = sort_params_by_idx(updates_in)
        self.restored_update_params = sorted_updates_in

    def risk(self,deterministic = False):
        return lasagne.layers.get_output(self.network, deterministic = deterministic)

    def predict_risk(self, x):
        risk_fxn = theano.function(
            inputs = [self.X],
            outputs = self.risk(deterministic= True),
            name = 'predicted risk'
        )
        return risk_fxn(x)

    def recommend_treatment(self, x, trt_i, trt_j, trt_idx = -1):
        x_trt = numpy.copy(x)
        x_trt[:,trt_idx] = trt_i
        h_i = self.predict_risk(x_trt)
        x_trt[:,trt_idx] = trt_j;
        h_j = self.predict_risk(x_trt)

        rec_ij = h_i - h_j
        return rec_ij

    def plot_risk_surface(self, data, i = 0, j = 1,
        figsize = (6,4), x_lims = None, y_lims = None, c_lims = None):
        fig = plt.figure(figsize=figsize)
        X = data[:,i]
        Y = data[:,j]
        Z = self.predict_risk(data)

        if not x_lims is None:
            x_lims = [np.round(np.min(X)), np.round(np.max(X))]
        if not y_lims is None:
            y_lims = [np.round(np.min(Y)), np.round(np.max(Y))]
        if not c_lims is None:
            c_lims = [np.round(np.min(Z)), np.round(np.max(Z))]

        ax = plt.scatter(X,Y, c = Z, edgecolors = 'none', marker = '.')
        ax.set_clim(*c_lims)
        plt.colorbar()
        plt.xlim(*x_lims)
        plt.ylim(*y_lims)
        plt.xlabel('$x_{%d}$' % i, fontsize=18)
        plt.ylabel('$x_{%d}$' % j, fontsize=18)

        return fig

def load_model_from_json(model_fp, weights_fp = None):
    with open(model_fp, 'r') as fp:
        json_model = fp.read()
    print('Loading json model:',json_model)
    hyperparams = json.loads(json_model)

    model = DeepSurv(**hyperparams)

    if weights_fp:
        model.load_weights(weights_fp)

    return model

# -*- coding: utf-8 -*-
"""
Created on Fri Aug 30 12:02:42 2019

@author: Andrew Wentzel
"""
from numpy.random import seed
seed(1)
from tensorflow.compat.v1 import set_random_seed
set_random_seed(2)

from PatientSet import PatientSet
from Constants import Constants
import Metrics
from analysis import *
from Models import *
from Toxicity import ClusterStats
from copy import copy
import numpy as np
import pandas as pd
from Boruta import BorutaPy
from sklearn.feature_selection import mutual_info_classif, f_classif, SelectPercentile, SelectKBest
from sklearn.model_selection import cross_validate, cross_val_predict, LeaveOneOut
from sklearn.metrics import accuracy_score, recall_score, roc_auc_score, roc_curve, f1_score
from sklearn.naive_bayes import BernoulliNB, ComplementNB, GaussianNB
from sklearn.linear_model import LogisticRegression
from sklearn.neighbors import KNeighborsClassifier
from sklearn.preprocessing import OneHotEncoder, QuantileTransformer, PowerTransformer
from sklearn.cluster import AgglomerativeClustering
from sklearn.tree import DecisionTreeClassifier
from sklearn.discriminant_analysis import QuadraticDiscriminantAnalysis
from sklearn.ensemble import VotingClassifier, ExtraTreesClassifier, RandomForestClassifier, AdaBoostClassifier
from sklearn.svm import SVC
from NCA import NeighborhoodComponentsAnalysis
from sklearn.base import BaseEstimator, ClassifierMixin
from collections import namedtuple, OrderedDict
from imblearn import under_sampling, over_sampling, combine
from scipy.special import softmax
from scipy.stats import kruskal

import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=UserWarning)
#warnings.filterwarnings("ignore", category=RuntimeWarning)

class MetricLearningClassifier(BaseEstimator, ClassifierMixin):

    def __init__(self, n_components = None,
                 random_state = 1,
                 resampler = None,
                 use_softmax = True,
                 metric_learner = None):
        self.n_components = n_components
        self.transformer = metric_learner
        if metric_learner is None:
            self.transformer = NeighborhoodComponentsAnalysis(n_components = n_components,
                                                      max_iter = 1000,
                                                      random_state = random_state)
        self.group_parameters = namedtuple('group_parameters', ['means', 'inv_covariance', 'max_dist'])
        self.resampler = resampler
        self.use_softmax = use_softmax

    def fit(self, x, y):
        self.transformer.fit(x, y)
        self.groups = OrderedDict()
        if self.resampler is not None:
            xtemp, ytemp = self.resampler.fit_resample(x,y)
            if len(np.unique(ytemp)) == len(np.unique(y)):
                x = xtemp
                y = ytemp
        for group in np.unique(y):
            self.groups[group] = self.group_params(x, y, group)

    def group_params(self, x, y, group):
        targets = np.argwhere(y == group).ravel()
        x_target = self.transformer.transform(x[targets])
        fmeans = x_target.mean(axis = 0)
        inv_cov = np.linalg.pinv(np.cov(x_target.T))
        train_dists = self.mahalanobis_distances(x, self.group_parameters(fmeans, inv_cov, 0))
        parameters = self.group_parameters(fmeans, inv_cov, train_dists.max())
        return parameters

    def mahalanobis_distances(self, x, group):
        x_offset = self.transformer.transform(x) - group.means
        left_term = np.dot(x_offset, group.inv_covariance)
        mahalanobis = np.dot(left_term, x_offset.T).diagonal()
        return mahalanobis

    def predict_proba(self, x):
        all_distances = []
        for group_id, group_params in self.groups.items():
            distances = self.mahalanobis_distances(x, group_params)
            proximity = np.clip(1 - (distances/group_params.max_dist), 0.00001, 1)
            all_distances.append(proximity)
        output = np.hstack(all_distances).reshape(-1, len(self.groups.keys()))
        if self.use_softmax:
            output = softmax(output)
        else:
            output = output/output.sum(axis = 1).reshape(-1,1)
        return output

    def predict(self, x):
        labels = list(self.groups.keys())
        probs = self.predict_proba(self, x)
        max_probs =  np.argmax(probs, axis = 1).ravel()
        ypred = np.zeros(max_probs.shape).astype(np.dtype(labels[0]))
        for i in range(max_probs.shape[0]):
            ypred[i] = labels[max_probs[i]]
        return ypred[i]

    def fit_predict(self, x, y):
        self.fit(x,y)
        return self.predict(x)



class RecallBasedModel:

    def __init__(self, model, recall_threshold = None, feature_selection = False):
        self.model = copy(model)
        self.recall_threshold = recall_threshold
        self.probability_threshold = None
        self.feature_selection = feature_selection

    def fit(self, x, y):
        self.get_feature_args(x,y)
        xfit = x[:, self.features_to_use]
        self.model.fit(xfit, y.ravel())
        if self.recall_threshold is not None:
            self.tune_threshold(xfit, y.ravel())

    def predict(self, x):
        xsubset = x[:, self.features_to_use]
        if self.probability_threshold is not None:
            probs = self.model.predict_proba(xsubset)[:,1].ravel()
            prediction = probs >= self.probability_threshold
        else:
            prediction = self.model.predict(xsubset)
        return prediction.astype('bool')

    def fit_predict(self, x, y):
        self.fit(x,y)
        return self.predict(x,y);

    def tune_threshold(self, x, y):
        ypred = self.model.predict_proba(x)
        sorted_scores = sorted(ypred[:,1], key = lambda x: -x)
        threshold_i = 0
        ythresh = ypred[:,1] >= sorted_scores[threshold_i]
        while recall_score(y, ythresh) < self.recall_threshold and threshold_i < len(sorted_scores) - 1:
            threshold_i += 1
            ythresh = ypred[:,1] >= sorted_scores[threshold_i]
        self.probability_threshold = sorted_scores[threshold_i]
        self.all_thresholds = sorted_scores

    def increment_threshold(self, increment = 1):
        current_index = self.all_thresholds.index(self.probability_threshold)
        new_index = np.clip(0, current_index + increment, len(self.all_thresholds) -1)
        self.probability_threshold = self.all_thresholds[new_index]

    def get_feature_args(self, x, y, percentile = 80, k = 40):
        if self.feature_selection == 'info':
            info_score = mutual_info_classif(x, y)
            self.features_to_use = np.argwhere(info_score > 0).ravel()
            if len(self.features_to_use) <= 1:
                self.features_to_use = np.argwhere(x.std(axis = 0) > 0).ravel()
        elif self.feature_selection == 'percentile':
            selector = SelectPercentile(percentile = percentile)
            selector.fit(x,y)
            self.features_to_use = np.argwhere(selector.get_support()).ravel()
        elif self.feature_selection == 'kbest':
            k = np.min([int(np.ceil(percentile*x.shape[1]/100)), k])
            selector = SelectKBest(k = k).fit(x,y)
            self.features_to_use = np.argwhere(selector.get_support()).ravel()
        else:
            self.features_to_use = np.argwhere(x.std(axis = 0) > 0).ravel()

    def __str__(self):
        string = str(self.model)
        string += '\n num features ' + str(len(self.features_to_use))
        string += '\n threshold ' + str(self.probability_threshold)
        string += '\n num_in_threshold: ' + str(self.all_thresholds.index(self.probability_threshold))
        return string + '\n'

class StackedClassifier:

    def __init__(self, default_models,
                 min_misses = 0,
                 recall_threshold = None,
                 feature_selection = False,
                 num_feature_splits = 2):
        self.gen_models = lambda: [RecallBasedModel(copy(m), recall_threshold, feature_selection) for m in default_models]
        self.min_misses = min_misses
        self.num_models = len(default_models)
        self.num_feature_splits = num_feature_splits

    def fit(self, x, y):
        models = []
        x_groups = self.split_features(x, y)
        for x_set in x_groups:
            model_set = self.gen_models()
            model_group = []
            for model in model_set:
                model.fit(x_set, y)
                model_group.append(model)
            models.append(model_group)
        self.models = models

    def predict(self, x, min_votes = None):
        min_votes = self.num_models if min_votes is None else min_votes
        x_groups= self.split_features(x)
        assert(len(x_groups) == len(self.models))
        predictions = []
        for group in range(len(x_groups)):
            model_set = self.models[group]
            x_set = x_groups[group]
            for model in model_set:
                ypred = model.predict(x_set)
                predictions.append(ypred.reshape(-1,1))
        predictions = np.hstack(predictions)
        ypred = predictions.sum(axis = 1) >= min_votes
        return ypred

    def fit_predict(self, x,y, min_votes = None):
        self.fit(x,y)
        return self.predict(x, min_votes)

    def split_features(self, x,y = None):
        #feature clusteirng or whatever here
        #just split if train is false, find grups if true
        if self.num_feature_splits <= 1:
            return [x]
        if y is not None:
            info_score = np.nan_to_num(mutual_info_classif(x, y))
            scores = np.argsort(-info_score)
            self.feature_group_args = [[] for g in range(self.num_feature_splits)]
            current_feature = 0
            while current_feature < x.shape[1]:
                for group in range(self.num_feature_splits):
                    self.feature_group_args[group].append(scores[current_feature])
                    current_feature += 1
                    if current_feature >= x.shape[1]:
                        break
#            self.feature_group_args = []
#            clusterer = AgglomerativeClustering(n_clusters = self.num_feature_splits)
#            features = x.transpose() #should I regularize here?
#            groups = clusterer.fit_predict(features)
#            for g in np.unique(groups):
#                args = np.argwhere(groups == g).ravel()
#                self.feature_group_args.append(args)
        x_groups = []
        for group_args in self.feature_group_args:
            x_groups.append(x[:, group_args])
        return x_groups

    def __repr__(self):
        string = ""
        for modelset in self.models:
            for model in modelset:
                string += str(model)
        return string

class IterativeClassifier:

    def __init__(self, default_models,
                 min_misses = 0,
                 recall_threshold = None,
                 feature_selection = False,
                 num_feature_splits = 1):
        self.gen_ensemble = lambda: StackedClassifier(default_models,
                                                      min_misses,recall_threshold,
                                                      feature_selection,
                                                      num_feature_splits)
        self.min_misses = min_misses
        self.num_models = len(default_models)
        self.num_feature_splits = num_feature_splits

    def fit(self, x, y):
        current_model = self.gen_ensemble()
        models = [current_model]
        y_pred = current_model.fit_predict(x,y)
        while not self.is_done(y, y_pred) and len(models) < 7:
            args = np.argwhere(y_pred == 0).ravel()
#            print(len(args))
            if len(args) < 5:
                break
            xsubset = x[args,:]
            ysubset = y[args]
            new_model = self.gen_ensemble()
            new_predict = new_model.fit_predict(xsubset, ysubset)
            new_args = np.argwhere(new_predict > 0)
            if len(new_args) <= 0:
                break
            y_pred[args[new_args]] = True
            models.append(new_model)
        self.models = models
#        print()

    def predict(self, x):
        y_pred = np.zeros((x.shape[0],))
        for model in self.models:
            valid = np.argwhere(model.predict(x) > 0).ravel()
            y_pred[valid] = 1
        return y_pred

    def predict_proba(self, x):
        y_pred = np.zeros((x.shape[0],))
        model_pos = 1
        for model in self.models:
            valid = np.argwhere(model.predict(x) > 0).ravel()
            y_pred[valid] = 1/model_pos
            model_pos += 1
        return y_pred

    def fit_predict(self, x, y):
        self.fit(x,y)
        return self.predict(x)

    def is_done(self, y, y_pred):
        false_predictions = np.argwhere(y_pred == 0).ravel()
        false_negatives = y[false_predictions].sum()
        return (false_negatives <= self.min_misses)

    def test(self, x, y):
        ypred = self.fit_predict(x, y)
        yproba = self.predict_proba(x)
        print('AUC score', roc_auc_score(y, yproba))
        print('recall ', recall_score(y, ypred))
        print('f1 ', f1_score(y, ypred))
        return yproba


def feature_matrix(db):
    discrete_dists = Metrics.discretize(-db.tumor_distances, n_bins = 15, strategy='uniform')
    t_volumes = np.array([np.sum([g.volume for g in gtvs]) for gtvs in db.gtvs]).reshape(-1,1)
    discrete_volumes = Metrics.discretize(t_volumes, n_bins = 15, strategy='uniform')

    discrete_sim = Metrics.augmented_sim(discrete_dists, Metrics.jaccard_distance)
#    predicted_doses = TreeKnnEstimator().predict_doses([discrete_sim], db)

    x = np.hstack([
        TreeKnnEstimator().predict_doses([discrete_sim], db),
        discrete_dists,
        discrete_volumes,
        db.prescribed_doses.reshape(-1,1),
        db.dose_fractions.reshape(-1,1),
        db.has_gtvp.reshape(-1,1),
        OneHotEncoder(sparse = False).fit_transform(db.lateralities.reshape(-1,1)),
        OneHotEncoder(sparse = False).fit_transform(db.subsites.reshape(-1,1)),
        OneHotEncoder(sparse = False).fit_transform(db.t_categories.reshape(-1,1)),
               ])
    return x

def cv_ensemble(features,
                toxicity,
                classifiers = None,
                num_feature_splits = 1,
                recall_threshold = .99,
                feature_selection = False,
                regularizer = None):
    if isinstance(features, PatientSet):
        features = feature_matrix(features)
    if classifiers is None:
        classifiers = [
                BernoulliNB(),
                LogisticRegression(solver = 'lbfgs',max_iter=3000)
                       ]
    ensemble = IterativeClassifier(classifiers,
                                   num_feature_splits = num_feature_splits,
                                   recall_threshold = recall_threshold,
                                   feature_selection = feature_selection)
    copy(ensemble).test(features, toxicity)
    loo = LeaveOneOut()
    ypred = []
    for train_index, test_index in loo.split(features):
        if regularizer is None:
            xtrain, xtest = Metrics.normalize_and_drop(features[train_index] , features[test_index])
        else:
            xtrain = regularizer.fit_transform(features[train_index])
            xtest = regularizer.transform(features[test_index])
        ensemble.fit(xtrain, toxicity[train_index])
        ypred.append(ensemble.predict_proba(xtest)[0])
    ypred = np.array(ypred)
    ypred /= ypred.max()
    print(roc_auc_score(toxicity, ypred), recall_score(toxicity, ypred > 0))
    print(ClusterStats().get_contingency_table(ypred, toxicity))

def roc_cv(classifier, x, y, feature_selection_method = None, regularizer = None, resampler = None):
    x = feature_selection(x, y, feature_selection_method)
    loo = LeaveOneOut()
    ypred = np.zeros(y.shape)
    regularizer = QuantileTransformer() if regularizer is None else regularizer
    for train_idx, test_idx in loo.split(x):
        xtrain, ytrain = x[train_idx], y[train_idx]
        if resampler is not None:
            xtrain, ytrain = resampler.fit_resample(xtrain, ytrain)
        xtrain = regularizer.fit_transform(xtrain)
        xtest = regularizer.transform(x[test_idx])
        classifier.fit(xtrain, ytrain)
        ypred[test_idx] = classifier.predict_proba(xtest)[0,1]
    results = {}
    results['AUC'] = roc_auc_score(y, ypred)
    results['f1'] = f1_score(y, ypred > .5)
    results['output'] = ypred
    results['resampler'] = resampler
    results['regularizer'] = regularizer
    results['classifier'] = classifier.fit(x,y)
    return results

def organ_data_to_df(arr,
                     ids = None,
                     suffix = '_distance',
                     organ_list = None,
                     to_merge = None):
    organ_list = organ_list if organ_list is not None else Constants.organ_list
    assert(arr.shape[1] == len(organ_list))
    columns = [o+suffix for o in organ_list]
    if to_merge is not None and ids is None:
        ids = to_merge.index.values
    df = pd.DataFrame(arr, index = ids, columns = columns)
    if to_merge is not None:
        df = df.join(to_merge, how = 'inner')
    return df

def feature_selection(x,y, method):
    return x

def generate_features(db = None,
                      file = 'data/ds_topFeatures.csv',
                      db_features = ['hpv'],
                      db_organ_list = None):
    if db is None:
        db = PatientSet(root = 'data\\patients_v*\\',
                        use_distances = False)
    df = pd.read_csv(file, index_col = 1).drop('Unnamed: 0', axis = 1)
    df['T.category'] = df['T.category'].apply(lambda x: int(x[1]))
    ft = df.FT.values
    ar = df.AR.values
    tox = df.TOX.values
    df = db.to_dataframe(db_features, df, organ_list = db_organ_list)
    df = df.drop(['FT', 'AR', 'TOX'], axis = 1)
    return  df, ft, ar, tox

def test_classifiers(db = None, log = False,
                     db_features = ['hpv'],
                     predicted_doses = None,
                     pdose_organ_list = None,
                     db_features_organ_list = None,
                     regularizer = QuantileTransformer(),
                     additional_features = None):

    if log:
        from time import time
        from datetime import datetime
        timestamp = datetime.fromtimestamp(time()).strftime('%Y_%m_%d_%H%M%S')
        f = open(Constants.toxicity_log_file_root + timestamp +'.txt', 'w', buffering = 1)
        def write(string):
            print(string)
            f.write(str(string)+'\n')
    else:
        write = lambda string: print(string)
    df, ft, ar, tox = generate_features(db,
                                        db_features = db_features,
                                        db_organ_list = db_features_organ_list)

    if predicted_doses is not None:
        if not isinstance(predicted_doses, np.ndarray):
            predicted_doses = default_rt_prediction(db)
        o_args = np.array([Constants.organ_list.index(o) for o in pdose_organ_list if o in Constants.organ_list])
        df = organ_data_to_df(predicted_doses[:,o_args],
                          suffix = '_pdoses',
                          to_merge = df,
                          organ_list = pdose_organ_list)

    write('features: ' + ', '.join([str(c) for c in df.columns]) + '\n')
    from xgboost import XGBClassifier
    classifiers = [
                    DecisionTreeClassifier(),
                    XGBClassifier(booster = 'gblinear'),
                    XGBClassifier(10, booster = 'gblinear'),
                    XGBClassifier(14, booster = 'gblinear'),
                    XGBClassifier(20),
                    LogisticRegression(C = 1, solver = 'lbfgs', max_iter = 3000),
                    MetricLearningClassifier(use_softmax = True),
                    MetricLearningClassifier(
                            resampler = under_sampling.OneSidedSelection()),
                    MetricLearningClassifier(
                            resampler = under_sampling.CondensedNearestNeighbour()),
                    ExtraTreesClassifier(n_estimators = 200),
                   ]
    results = []
    resamplers = [None,
#                  under_sampling.InstanceHardnessThreshold(
#                          estimator = MetricLearningClassifier(),
#                          cv = 18),
#                  under_sampling.InstanceHardnessThreshold(cv = 18),
                  over_sampling.SMOTE(),
                  combine.SMOTEENN(),
                  under_sampling.InstanceHardnessThreshold(),
                  under_sampling.RepeatedEditedNearestNeighbours(),
                  under_sampling.EditedNearestNeighbours(),
                  under_sampling.CondensedNearestNeighbour(),
                  ]

    for classifier in classifiers:
        write(classifier)
        for resampler in resamplers:
            write(resampler)
            for outcome in [(ft, 'feeding_tube'), (ar, 'aspiration')]:
                try:
                    roc = roc_cv(classifier, df.values, outcome[0],
                                 regularizer = regularizer,
                                 resampler = resampler)
                    roc['outcome'] = outcome[1]
                    write(outcome[1])
                    write(roc['AUC'])
                    results.append(roc)
                except Exception as e:
                    print(e)
            write('\n')
    if log:
        f.close()

def plot_correlations(db,
                      use_predicted_dose = True,
                      use_distances = True,
                      use_volumes = True,
                      max_p = .15,
                      tox_name = 'toxicity'):
    db_features = ['hpv', 'smoking', 'ages',
                   'packs_per_year', 'dose_fractions',
                   'prescribed_doses', 'has_gtvp']
    data, ft, ar, tox = generate_features(db, db_features =db_features)
    if tox_name == 'toxicity':
        y = tox
    elif tox_name in ['feeding_tube', 'ft']:
        y = ft
    elif tox_name in ['aspiration', 'ar', 'aspiration_rate']:
        y = ar
    if use_predicted_dose:
        pred_doses = default_rt_prediction(db)
        data = organ_data_to_df(pred_doses,
                                ids = db.ids, to_merge = data,
                                suffix = '_pred_dose')
    if use_distances:
        data = organ_data_to_df(db.tumor_distances,
                                ids = db.ids, to_merge = data,
                                suffix = '_tumor_distance')
    if use_volumes:
        data = organ_data_to_df(db.volumes,
                                ids = db.ids, to_merge = data,
                                suffix = '_volumes')
    clusters = data.hc_ward.values
    high_cluster = np.argwhere(clusters == 2).ravel()
    toxicity = np.argwhere(y > 0).ravel()
    outliers = set(high_cluster) - set(toxicity)
    inliers = set(high_cluster) - outliers
    data = data.assign(classes = pd.Series(-np.ones(y.shape), index = data.index).values)
    data.classes.iloc[sorted(inliers)] = 0
    data.classes.iloc[sorted(outliers)] = 1
    inlier_data = (data[data.classes == 0])
    outlier_data = (data[data.classes == 1])
    pvals = {}
    for col in sorted(data.drop(['hc_ward', 'classes'], axis = 1).columns):
        v1 = inlier_data[col].values
        v2 = outlier_data[col].values
        pval = kruskal(v1, v2).pvalue
        if pval < max_p:
            pvals[col] = pval
    sorted_pvals = sorted(pvals.items(), key = lambda x: x[1])
    print(sorted_pvals)
    labels, vals = list(zip(*sorted_pvals))
    plt.barh(np.arange(len(vals)), max_p-np.array(vals), tick_label = labels, left = 1-max_p)
    plt.xlabel('1 - pvalue for kruskal-wallis test')
    plt.title('1 - pvalue between cluster 2 with and without ' + tox_name + ' per feature')

db = PatientSet(root = 'data\\patients_v*\\',
                    use_distances = False)
plot_correlations(db, tox_name = 'feeding_tube', max_p = .25)
#p_doses = default_rt_prediction(db)
pdose_organs = ['Soft_Palate', 'SPC', 'Extended_Oral_Cavity', 'Hard_Palate', 'Mandible', 'Brainstem', 'Lower_Lip']
feature_organs = ['Lt_Masseter_M', 'Rt_Masseter_M']
test_classifiers(db, log = True,
                 db_features = ['hpv', 'volumes', 'smoking'],
                 predicted_doses = p_doses,
                 pdose_organ_list = pdose_organs,
                 db_features_organ_list = feature_organs)

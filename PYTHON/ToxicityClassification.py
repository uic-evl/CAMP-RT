# -*- coding: utf-8 -*-
"""
Created on Fri Aug 30 12:02:42 2019

@author: Andrew Wentzel
"""

from PatientSet import PatientSet
from Constants import Constants
import Metrics
from analysis import *
from Models import *
from copy import copy
import numpy as np
from sklearn.feature_selection import mutual_info_classif, f_classif, SelectPercentile
from sklearn.model_selection import cross_validate, cross_val_predict, LeaveOneOut
from sklearn.metrics import accuracy_score, recall_score, roc_auc_score, roc_curve, f1_score
from sklearn.naive_bayes import BernoulliNB, ComplementNB, GaussianNB
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import OneHotEncoder
from sklearn.cluster import AgglomerativeClustering

class RecallBasedModel:

    def __init__(self, model, recall_threshold = None, feature_selection = False):
        self.model = copy(model)
        self.recall_threshold = recall_threshold
        self.probability_threshold = None
        self.feature_selection = feature_selection

    def fit(self, x, y):
        if self.feature_selection:
            self.get_feature_args(x,y)
        else:
            self.features_to_use = np.arange(x.shape[1]).ravel()
        xfit = x[:, self.features_to_use]
        self.model.fit(xfit, y.ravel())
        if self.recall_threshold is not None:
            self.tune_threshold(xfit, y.ravel())

    def predict(self, x):
        xsubset = x[:, self.features_to_use]
        if self.probability_threshold is not None:
            probs = self.model.predict_proba(xsubset)[:,1].ravel()
            prediction = probs > self.probability_threshold
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
        ythresh = ypred[:,1] > sorted_scores[threshold_i]
        while recall_score(y, ythresh) < self.recall_threshold and threshold_i < len(sorted_scores) - 1:
            threshold_i += 1
            ythresh = ypred[:,1] > sorted_scores[threshold_i]
        self.probability_threshold = sorted_scores[threshold_i]
        self.all_thresholds = sorted_scores

    def increment_threshold(self, increment = 1):
        current_index = self.all_thresholds.index(self.probability_threshold)
        new_index = np.clip(0, current_index + increment, len(self.all_thresholds) -1)
        self.probability_threshold = self.all_thresholds[new_index]

    def get_feature_args(self, x, y):
        if self.feature_selection == 'info':
            info_score = mutual_info_classif(x, y)
            self.features_to_use = np.argwhere(info_score > 0).ravel()
        else:
            selector = SelectPercentile(percentile = 80)
            selector.fit(x,y)
            self.features_to_use = np.argwhere(selector.get_support()).ravel()

class StackedClassifier:

    def __init__(self, default_models,
                 min_misses = 0,
                 recall_threshold = None,
                 feature_selection = False,
                 num_feature_splits = 2):
        self.gen_models = lambda: [RecallBasedModel(m, recall_threshold, feature_selection) for m in default_models]
        self.min_misses = min_misses
        self.num_models = len(default_models)
        self.num_feature_splits = num_feature_splits

    def fit(self, x, y):
        models = []
        x_groups = self.split_features(x, train = True)
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
        x_groups= self.split_features(x, train = False)
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

    def split_features(self, x, train = False):
        #feature clusteirng or whatever here
        #just split if train is false, find grups if true
        if self.num_feature_splits <= 1:
            return [x]
        if train:
            self.feature_group_args = []
            clusterer = AgglomerativeClustering(n_clusters = self.num_feature_splits)
            features = x.transpose() #should I regularize here?
            groups = clusterer.fit_predict(features)
            for g in np.unique(groups):
                args = np.argwhere(groups == g).ravel()
                self.feature_group_args.append(args)
        x_groups = []
        for group_args in self.feature_group_args:
            x_groups.append(x[:, group_args])
        return x_groups

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
        while not self.is_done(y, y_pred) and len(models) < 5:
            args = np.argwhere(y_pred == 0).ravel()
            xsubset = x[args,:]
            ysubset = y[args]
            new_model = self.gen_ensemble()
            new_predict = new_model.fit_predict(xsubset, ysubset)
            new_args = np.argwhere(new_predict > 0)
            y_pred[args[new_args]] = True
            models.append(new_model)
        self.models = models

    def predict(self, x):
        y_pred = np.zeros((x.shape[0],))
        for model in self.models:
            valid = np.argwhere(model.predict(x) > 0).ravel()
            y_pred[valid] = 1
        return y_pred

    def fit_predict(self, x, y):
        self.fit(x,y)
        return self.predict(x)

    def is_done(self, y, y_pred):
        false_predictions = np.argwhere(y_pred == 0).ravel()
        false_negatives = y[false_predictions].sum()
        return (false_negatives <= self.min_misses)


def feature_matrix(db):
    discrete_dists = Metrics.discretize(-db.tumor_distances, n_bins = 15, strategy='uniform')
    t_volumes = np.array([np.sum([g.volume for g in gtvs]) for gtvs in db.gtvs]).reshape(-1,1)
    discrete_volumes = Metrics.discretize(t_volumes, n_bins = 15, strategy='uniform')

    discrete_sim = Metrics.augmented_sim(discrete_dists, Metrics.jaccard_distance)
    predicted_doses = TreeKnnEstimator().predict_doses([discrete_sim], db)

    x = np.hstack([
        discrete_dists,
        discrete_volumes,
        db.prescribed_doses.reshape(-1,1),
        db.dose_fractions.reshape(-1,1),
        db.has_gtvp.reshape(-1,1),
        OneHotEncoder(sparse = False).fit_transform(db.lateralities.reshape(-1,1)),
        OneHotEncoder(sparse = False).fit_transform(db.subsites.reshape(-1,1)),
        predicted_doses
        #OneHotEncoder(sparse = False).fit_transform(db.t_categories.reshape(-1,1)),
               ])
    return x

#db = PatientSet(root = 'data\\patients_v*\\',
#                use_distances = False)

#toxicity = (db.feeding_tubes + db.aspiration) > 0
#features = feature_matrix(db)
classifiers = [BernoulliNB(), LogisticRegression(solver = 'lbfgs',max_iter=6000,class_weight='balanced')]

ensemble = IterativeClassifier(classifiers, recall_threshold = .6)
loo = LeaveOneOut()
ypred = []
for train_index, test_index in loo.split(features):
    ensemble.fit(features[train_index], toxicity[train_index])
    ypred.append(ensemble.predict(features[test_index])[0])
ypred = np.array(ypred)
print(roc_auc_score(toxicity, ypred), recall_score(toxicity, ypred))
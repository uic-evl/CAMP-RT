# -*- coding: utf-8 -*-
"""
Created on Thu Aug  8 13:11:50 2019

@author: Andrew Wentzel
"""
from numpy.random import seed
seed(1)
import numpy as np
import matplotlib.pyplot as plt
from analysis import *
from sklearn.cluster import KMeans, AgglomerativeClustering
from collections import namedtuple
import Metrics
from PatientSet import PatientSet
from Constants import Constants
import rpy2.robjects.numpy2ri
from rpy2.robjects.packages import importr
rpy2.robjects.numpy2ri.activate()
from sklearn.manifold import MDS
from NCA import NeighborhoodComponentsAnalysis
from Boruta import BorutaPy

from sklearn.naive_bayes import BernoulliNB, ComplementNB, GaussianNB
from sklearn.preprocessing import OneHotEncoder, quantile_transform
from sklearn.model_selection import cross_validate, cross_val_predict, LeaveOneOut
from sklearn.metrics import accuracy_score, recall_score, roc_auc_score, roc_curve, f1_score
from sklearn.ensemble import ExtraTreesClassifier, VotingClassifier
from sklearn.feature_selection import mutual_info_classif, f_classif, SelectKBest, mutual_info_regression
from sklearn.linear_model import LogisticRegression
from sklearn.manifold import MDS
from sklearn.utils import resample

from imblearn.combine import SMOTEENN, SMOTETomek

cluster_result = namedtuple('cluster_result', ['method', 'cluster', 'correlation', 'model'])

class ClusterStats():

    def __init__(self, clusterers = None, similarity_func = None, upsample = False):
        self.clusterers = clusterers
        self.upsample = upsample
        self.similarity_func = similarity_func
        self.mds = MDS(n_components = 30, dissimilarity = 'precomputed')
        if clusterers is None:
            c_range = range(2,5)
            self.clusterers = {}
            self.clusterers['Kmeans'] = [KMeans(n_clusters = i) for i in c_range]
            self.clusterers['Agglomerative_ward'] = [AgglomerativeClustering(n_clusters = i) for i in c_range]
            self.clusterers['Agglomerative_complete'] = [AgglomerativeClustering(n_clusters = i, linkage = 'complete') for i in c_range]

    def fisher_exact_test(self, c_labels, y):
        if len(set(y)) == 1:
            print('fisher test run with no positive class')
            return 0
#        assert(len(set(y)) == 2)
        #call fishers test from r
        contingency = self.get_contingency_table(c_labels, y)
        stats = importr('stats')
        pval = stats.fisher_test(contingency)[0][0]
        return pval

    def get_contingency_table(self, x, y):
        #assumes x and y are two equal length vectors, creates a mxn contigency table from them
        cols = sorted(list(np.unique(y)))
        rows = sorted(list(np.unique(x)))
        tabel = np.zeros((len(rows), len(cols)))
        for row_index in range(len(rows)):
            row_var = rows[row_index]
            for col_index in range(len(cols)):
                rowset = set(np.argwhere(x == row_var).ravel())
                colset = set(np.argwhere(y == cols[col_index]).ravel())
                tabel[row_index, col_index] = len(rowset & colset)
        return tabel

    def analyze_clusters(self, target_var, name, clusterer, doses, subset):
        result = []
        distance = self.get_dose_embedding(doses, target_var, subset)
        clusters = clusterer.fit_predict(distance).ravel()
        n_clusters = len(set(clusters))
        method = name + str(n_clusters)

        overall_correlation = self.fisher_exact_test(clusters, target_var)
        result.append( cluster_result(method, 'all',
                                      overall_correlation,
                                      clusterer))
#        print(method, overall_correlation)

        for c in np.unique(clusters):
            correlation = self.fisher_exact_test(clusters == c, target_var)
            result.append( cluster_result(method, str(c+1),
                                          correlation, clusterer))
        return result

    def get_dose_embedding(self, features, outcome, subset = True):
        if subset:
            features = self.subset_features(features, outcome)
        if self.similarity_func is not None:
            similarity = Metrics.dist_to_sim(Metrics.reduced_augmented_sim(features, self.similarity_func))
            features = self.mds.fit_transform(similarity)
        return features

    def cluster_by_dose(self, target_var, doses, args = None, subset = True):
        if args is not None:
            assert( isinstance(args, list) )
            doses = doses[:, args]
        results = []
        for cname, clusterers in self.clusterers.items():
            for clusterer in clusterers:
                results.extend(self.analyze_clusters(target_var, cname, clusterer, doses, subset))
        results = sorted(results, key = lambda x: x.correlation)
        return results

    def get_optimal_clustering(self, doses, target_var, args = None,
                               subset = False, patient_subset = None):
        clusters = np.zeros(target_var.shape)
        if patient_subset is not None:
            target = target_var[patient_subset]
            doses = doses[patient_subset,:]
        else:
            target = target_var
        result = self.cluster_by_dose(target, doses,
                                      args, subset)
        result = [r for r in result if r.cluster is 'all']
        if args is not None:
            doses = doses[:, args]
        clusters[patient_subset] = result[0].model.fit_predict(doses).ravel() + 1
        pval = self.fisher_exact_test(clusters, target_var)
        clusterer_data = cluster_result(method = result[0].method,
                                        cluster = result[0].cluster,
                                        correlation = pval,
                                        model = result[0].model)
        optimal = (clusters, clusterer_data)
        return optimal

    def subset_features(self, x, y):
        mutual_info = mutual_info_classif(x, y)
        good_features = np.argwhere(mutual_info > 0).ravel()
        return x[:, good_features]

def recall_based_model(x, y, model, score_threshold = .99, selector = None):
    loo = LeaveOneOut()
    loo.get_n_splits(x)
    y_out = np.zeros(y.shape)
    resampler = SMOTETomek()
    if len(np.argwhere(y > 0)) <= 1:
        return y_out
    for train_index, test_index in loo.split(x):
        xtrain, ytrain = x[train_index], y[train_index]
        if len(set(ytrain)) <= 1:
            y_out[test_index] = ytrain[0]
            continue
        xtest = x[test_index]
        if selector is not None:
            if selector == 'info':
                xtrain, xtest = select_by_info(xtrain, ytrain, xtest)
            elif selector == 'correlation':
                xtrain, xtest = select_by_fscore(xtrain, ytrain, xtest)
            else:
                fitted = selector.fit_transform(xtrain, ytrain)
                if fitted.ndim > 1 and fitted.shape[1] > 1:
                    xtrain = fitted
                    xtest = selector.transform(xtest)
        xfit, yfit = resampler.fit_resample(xtrain, ytrain.ravel())
        model.fit(xfit, yfit)
        yfit_pred = model.predict_proba(xtrain)
        sorted_scores = sorted(yfit_pred[:, 1], key = lambda x: -x)
        threshold_i = 0
        y_pred = yfit_pred[:,1] >= sorted_scores[threshold_i]
        while recall_score(y[train_index], y_pred) < score_threshold:
            threshold_i = threshold_i + 1
            y_pred = yfit_pred[:,1] >= sorted_scores[threshold_i]
        y_out[test_index] = model.predict_proba(xtest)[:,1] >= sorted_scores[threshold_i]
    return y_out

def ensemble_recall_based_model(x, y, model_list, score_threshold = 1, selector = None):
    if not isinstance(model_list, list):
        return recall_based_model(x, y, model_list, score_threshold, selector)
    labels = np.zeros(y.shape)
    for model in model_list:
        labels = labels + recall_based_model(x, y, model,
                                             score_threshold,
                                             selector = selector)
    return labels == len(model_list)

def select_by_info(xtrain, ytrain, xtest, threshold = 0):
    mutual_info = mutual_info_classif(xtrain, ytrain)
    good_features = np.argwhere(mutual_info > threshold).ravel()
    if xtest.ndim >1:
        xout = xtest[:,good_features]
    else:
        xout = xtest[good_features]
    return xtrain[:,good_features], xout

def select_by_fscore(xtrain, ytrain, xtest, threshold = .33):
    if threshold >= 1:
        return xtrain, xtest
    pvals = f_classif(xtrain, ytrain)[1]
    good_features = np.argwhere(pvals < threshold).ravel()
    if len(good_features) <= 1:
        return select_by_fscore(xtrain, ytrain, xtest, threshold + .1)
    if xtest.ndim >1:
        xout = xtest[:,good_features]
    else:
        xout = xtest[good_features]
    return xtrain[:,good_features], xout

def cv_feature_selection(x, y, method = 'info', **kwargs):
    if method == 'boruta':
        return cv_boruta_feature_selection(x, y, **kwargs)
    else:
        return cv_info_feature_selection(x, y, **kwargs)

def cv_info_feature_selection(x, y, percentile = 80, k = 40):
    n_patients = len(y)
    x_out = []
    k = np.min([int(np.ceil(percentile*x.shape[1]/100)), k])
    selector = SelectKBest(k = k)
    for p in range(n_patients):
        xtrain = np.delete(x, p, axis = 0)
        ytrain = np.delete(y, p, axis = 0)
        selector.fit(xtrain, ytrain)
        xfit = selector.transform(x[p,:].reshape(1,-1)).ravel()
        x_out.append( xfit )
    return np.vstack(x_out)

def cv_boruta_feature_selection(x, y, use_weak = True):
    n_patients = len(y)
    x_out = np.zeros(x.shape)
    boruta = BorutaPy( ExtraTreesClassifier(500, max_depth = 7),
                      n_estimators = 'auto')
    for p in range(n_patients):
        xtrain = np.delete(x, p, axis = 0)
        ytrain = np.delete(y, p, axis = 0)
        xtest = x[p]
        boruta.fit(xtrain, ytrain)
        mask = boruta.support_
        if use_weak:
            mask = mask | boruta.support_weak_
        x_out[p] = xtest*mask
    nonzero = np.argwhere( x_out.sum(axis = 0) > 0).ravel()
    x_out = x_out[:, nonzero]
    return x_out

def get_rankings(x, y, tiers, models):
    if isinstance(models, list):
        estimators = [('estimator' + str(k), models[k]) for k in range(len(models))]
        models = VotingClassifier(estimators, voting = 'soft')
    classes = np.unique(tiers)
    loo = LeaveOneOut()
    x = cv_feature_selection(x, y)
    scores = cross_val_predict(models, x, y, cv = loo, method='predict_proba')[:,1]
    for c in classes:
        class_args = np.argwhere(tiers == c).ravel()
        scores[class_args] += c
    return Metrics.minmax_scale(scores)

def get_labels(x, y, models,depth,
               selector = None,
               min_size = 10,
               min_ratio = 0,
               feature_selection = 'info'):
    x_subset = cv_feature_selection(x, y, feature_selection)
    labels = ensemble_recall_based_model(x_subset, y,
                                         models,
                                         selector = selector).astype('int32')
    positives = np.argwhere(labels > 0).ravel()
    negatives = np.argwhere(labels < .0001).ravel()
    n_pos, n_neg = len(positives), len(negatives)
    if depth < 4 and min([n_pos, n_neg]) > min_ratio * len(y):
        if n_pos > min_size:
            labels[positives] += get_labels(x[positives], y[positives],
                  models, depth+1, selector, min_size).astype('int32')
        if n_neg > min_size:
            labels[negatives] += get_labels(x[negatives], y[negatives],
                  models, depth+1, selector, min_size).astype('int32')
    return labels

def tiered_model(db, toxicity, models, feature_selection, selector = None):
    cluster_stats = ClusterStats()
    features = feature_matrix(db)
    new_classes = get_labels(features, toxicity, models, 0,
                             selector = selector,
                             feature_selection = feature_selection)

    new_classes = Metrics.discretize(new_classes.reshape(-1,1),
                             n_bins = len(set(new_classes)),
                             strategy = 'uniform').ravel()
    rankings = get_rankings(features, toxicity, new_classes, models)
    print(cluster_stats.get_contingency_table(new_classes, toxicity))
    print(cluster_stats.fisher_exact_test(new_classes, toxicity))
    print(roc_auc_score(toxicity, rankings))
    return new_classes, rankings

def evaluate_tiered_model(db, feature_selection = 'info'):
    toxicity = db.feeding_tubes + 2*db.aspiration
    labelers = [
            BernoulliNB(),
            ComplementNB(),
            GaussianNB(),
            LogisticRegression(solver = 'lbfgs', max_iter = 2000),
            ]

    ft_classes, ft_rankings = tiered_model(db, db.feeding_tubes, labelers, feature_selection)
    a_classes, a_rankings = tiered_model(db, db.aspiration, labelers, feature_selection)

    new_classes = 2*(a_classes > a_classes.min()) + (ft_classes > ft_classes.min())
    print(ClusterStats().get_contingency_table(new_classes, toxicity))
    print(ClusterStats().fisher_exact_test(new_classes, toxicity))
    print(f1_score(toxicity, new_classes, average = 'micro'), f1_score(toxicity, new_classes, average = 'macro'))


def feature_matrix(db):
    discrete_dists = Metrics.discretize(-db.tumor_distances, n_bins = 15, strategy='uniform')
    t_volumes = np.array([np.sum([g.volume for g in gtvs]) for gtvs in db.gtvs]).reshape(-1,1)
    discrete_volumes = Metrics.discretize(t_volumes, n_bins = 15, strategy='uniform')
    x = np.hstack([
        discrete_dists,
        #discrete_volumes,
        db.prescribed_doses.reshape(-1,1),
        db.dose_fractions.reshape(-1,1),
        #db.has_gtvp.reshape(-1,1),
        OneHotEncoder(sparse = False).fit_transform(db.lateralities.reshape(-1,1)),
        OneHotEncoder(sparse = False).fit_transform(db.subsites.reshape(-1,1)),
        OneHotEncoder(sparse = False).fit_transform(db.t_categories.reshape(-1,1)),
        pca(db.lymph_nodes[:, np.argwhere(db.lymph_nodes.std(axis = 0) > 0).ravel()],4)
               ])
    return x

def get_model_auc(x, y, model):
    ypred = cross_val_predict(model, x, y, cv = LeaveOneOut(), method = 'predict_proba')
    ypred = ypred[:,1]
    roc_score = roc_auc_score(y, ypred)
    fpr, tpr, thresholds = roc_curve(y, ypred)
    plt.plot(fpr, tpr)
    return fpr, tpr, thresholds, roc_score

def cluster_with_model(db, similarity, toxicity, model = None, selector = None):
    cluster_stats = ClusterStats()
    features = feature_matrix(db)
    features = cv_feature_selection(features, toxicity)
    predicted_doses = TreeKnnEstimator().predict_doses([similarity], db)

    if model is not None:
        labels = ensemble_recall_based_model(features, toxicity,
                                             model, selector = selector)
        if isinstance(model, list):
            print(recall_score(toxicity, labels))
        else:
            print(get_model_auc(features, toxicity, model)[-1])
        print(cluster_stats.get_contingency_table(labels, toxicity))
        print(cluster_stats.fisher_exact_test(labels, toxicity))
    else:
        labels = np.ones(toxicity.shape)
    cluster_results = cluster_stats.get_optimal_clustering(predicted_doses, toxicity,
                              patient_subset = np.argwhere(labels).ravel())
    print(cluster_results[1])
    print(KnnEstimator().get_error(predicted_doses, db.doses).mean())
    print(cluster_stats.get_contingency_table(cluster_results[0], toxicity))

    return cluster_results[0], labels

def rescale(x1, x2 = None):
    scale = lambda x: (x - x1.min(axis = 0))/(x1.max(axis = 0) - x1.min(axis = 0))
    if x2 is not None:
        return scale(x1), scale(x2)
    return scale(x1)

def nca_bootstrap_fit(xtrain, ytrain, xtest, n_samples = 100):
    nca = NeighborhoodComponentsAnalysis(5, max_iter=1000)
    transform = np.zeros((nca.n_components, xtrain.shape[1]))
    if n_samples == 0:
        nca.fit(xtrain, ytrain)
        return nca.transform(xtest).ravel()
    for n in range(n_samples):
        xnew, ynew = resample(xtrain, ytrain, stratify = ytrain)
        nca.fit(xnew, ynew)
        transform += nca.components_/n_samples
    return np.dot(xtest, transform.T).ravel()

def mutual_info_fit(xtrain, ytrain, xtest, n_samples = 0, pow_scale= 1/2):
    if n_samples == 0:
        info = mutual_info_classif(xtrain, ytrain)
        return xtest*(info**pow_scale)
    info = np.zeros(xtest.shape)
    for n in range(n_samples):
        xnew, ynew = resample(xtrain, ytrain, stratify = ytrain)
        info += mutual_info_classif(xnew, ynew)/n_samples
    return xtest*(info**pow_scale)

def info_regression_fit(xtrain, ytrain, xtest, n_samples = 0, pow_scale= 1/2):
    if n_samples == 0:
        info = mutual_info_regression(xtrain, ytrain)
        return xtest*(info**pow_scale)
    info = np.zeros(xtest.shape)
    for n in range(n_samples):
        xnew, ynew = resample(xtrain, ytrain, stratify = ytrain)
        info += mutual_info_regression(xnew, ynew)/n_samples
    return xtest*(info**pow_scale)

def svm_fit(xtrain, ytrain, xtest, n_samples = 0):
    from sklearn.svm import SVC
    svm = SVC(kernel = 'linear')
    if n_samples == 0:
        svm.fit(xtrain, ytrain)
        coef = svm.coef_/np.linalg.norm(svm.coef_)
        return xtest*coef
    coefs = np.zeros(xtest.shape)
    for n in range(n_samples):
        xnew, ynew = resample(xtrain, ytrain, stratify = ytrain)
        svm.fit(xnew, ynew)
        new_coefs = svm.coef_/np.linalg.norm(svm.coef_)
        coefs += new_coefs/n_samples
    return xtest*coefs

def tree_importance_fit(xtrain, ytrain, xtest, n_estimators = 100):
    if n_estimators == 0:
        n_estimators = 100 #this automatically does the bagging I think
    tree = ExtraTreesClassifier(n_estimators = n_estimators)
    tree.fit(xtrain, ytrain)
    return xtest*tree.feature_importances_

def boruta_fit(xtrain, ytrain, xtest, n_estimators = 0):
    boruta = BorutaPy(ExtraTreesClassifier(200, max_depth=5), n_estimators = 'auto')
    boruta.fit(xtrain, ytrain)
    if n_estimators == 0:
        return boruta.support_*xtest
    xout = []
    for n in range(n_estimators):
        xnew, ynew = resample(xtrain, ytrain, stratify = ytrain)
        boruta.fit(xtrain, ytrain)
        xout.append( boruta.support_*xtest/n_estimators )
    return np.sum(xout, axis = 0)



def select_supervised_features(known, estimated, toxicity, weight_func, n_samples = 0, outliers = None):
    transformed_predicted = []
    for p in range(toxicity.shape[0]):
        to_delete = p
        if outliers is not None:
            to_delete = outliers[:]
            to_delete.append(p)
        xtrain = np.delete(known, to_delete, axis = 0)
        ytrain = np.delete(toxicity, to_delete, axis = 0)
        xtest = estimated[p].reshape(1,-1)
        xtrain, xtest = rescale(xtrain, xtest)
        transformed = weight_func(xtrain, ytrain, xtest, n_samples)
        transformed_predicted.append(transformed)
    transformed_predicted = np.vstack(transformed_predicted)
    return transformed_predicted

from Models import *
from ErrorChecker import ErrorChecker

def save_prediction(transformed_predicted, name):
    transformed_predicted = transformed_predicted[:, np.argwhere(transformed_predicted.sum(axis = 0) > 0).ravel()]
    clustering = ClusterStats().get_optimal_clustering(transformed_predicted, toxicity)
    print(clustering[1])
    print(ClusterStats().get_contingency_table(clustering[0], toxicity))

    def get_class_db(db, classes):
        new_db = copy(db)
        new_db.classes = classes
        return new_db

    fancy_sim = Metrics.dist_to_sim(Metrics.reduced_augmented_sim(transformed_predicted, Metrics.mse))

    export(get_class_db(db, clustering[0]), similarity = fancy_sim, score_file = 'data\\' + name + '_similarity.csv')

#db = PatientSet(root = 'data\\patients_v*\\',
#                use_distances = False)


#toxicity = (db.feeding_tubes + db.aspiration) > 0
#
#discrete_dists = discretize(-db.tumor_distances)
#discrete_sim = Metrics.augmented_sim(discrete_dists, Metrics.jaccard_distance)
#
#predicted_doses = TreeKnnEstimator().predict_doses([discrete_sim], db)
#
#outliers = list(ErrorChecker().get_data_outliers(db.doses))
#
#true_doses = db.doses
#est_doses = predicted_doses

known = np.hstack([true_doses, feature_matrix(db)])
guessed = np.hstack([est_doses, feature_matrix(db)])


boruta = BorutaPy(ExtraTreesClassifier(300, max_depth=7), n_estimators = 300)
n_samples = 200
boruta.fit(rescale(known), toxicity)
knownfit = rescale(known)
xest = rescale(guessed)
support = np.zeros((knownfit.shape[1],))
weak_support = np.zeros(support.shape)
save_prediction(rescale(guessed)*boruta.support_, 'boruta_nobootstrap')
print(np.argwhere(boruta.support_ > 0).ravel())
for n in range(n_samples):
    xfit, y = resample(knownfit, toxicity, stratify = toxicity)
    boruta.fit(xfit, y)
    support += (boruta.support_)/(n_samples)
    weak_support += (boruta.support_weak_)/(n_samples)
    print(n, support*n)
supported = support > .66 #2 standard deviations?

save_prediction(rescale(guessed)*supported, 'boruta_bootstrap')





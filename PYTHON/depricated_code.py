# -*- coding: utf-8 -*-
"""
Created on Wed Jun 19 16:57:03 2019

@author: Andrew
"""

class NodeSimilarityModel():

    def __init__(self):
        pass

    def get_similarity(self, db):
        node_matrix = db.lymph_nodes
        subsites = db.subsites
        lateralities = db.lateralities
        num_patients = node_matrix.shape[0]
        similarity = np.zeros((num_patients, num_patients))
        for i1 in range(num_patients):
            for i2 in range(num_patients):
                if i1 == i2:
                    continue
                p1 = node_matrix[i1, :]
                p2 = node_matrix[i2, :]
                subsite1 = subsites[i1]
                subsite2 = subsites[i2]
                laterality1 = lateralities[i1]
                laterality2 = lateralities[i2]
                same_laterality = 1 if (laterality1 == laterality2) else 0
                same_subsite = 1 if subsite1 == subsite2 else 0
                similarity[i1, i2] = self.similarity(p1, p2, same_laterality, same_subsite)
        return similarity

    def similarity(self, x, y, j, k):
        numerator = x.dot(y)
        denominator = x.dot(x) + y.dot(y) - x.dot(y)
        if numerator == 0 or denominator == 0:
            return 0
        return numerator/denominator
    
class TreeEstimator():
    subsite_map = {'BOT': 0, 'GPS': 1, 'NOS': 2, 'Soft palate': 3, 'Tonsil': 4}
    laterality_map = {'B': 0, 'L': 1, 'R': 2}
    def __init__(self, num_pca_components = 10, n_estimators =10, min_samples_split = 4, max_depth = None):
        from sklearn.ensemble import RandomForestRegressor
        self.model = RandomForestRegressor(min_samples_split = min_samples_split, 
                                           n_estimators = n_estimators,
                                           max_depth = max_depth)
        self.num_pca_components = num_pca_components

    def evaluate(self, data):
        x = self.get_input_features(data)
        y = data.doses
        i = 0
        errors = []
        while i < (data.get_num_patients() - 10):
            values = np.arange(i, i+10, 1)
            x_train = np.delete(x, values, axis = 0)
            y_train = np.delete(y, values, axis = 0)
            x_test = x[values, :]
            y_test = y[values, :]
            i += 10
            y_predict = self.predict_doses( x_train, y_train, x_test, y_test)
            error = self.get_error(y_predict, y_test)
            errors.append(error)
        return np.mean(errors, axis = 0)
        
    def predict_doses(self, x_train, y_train, x_test, y_test):
        self.model.fit(x_train, y_train)
        y_predict = self.model.predict(x_test)
        return y_predict
        
    def get_error(self, predicted_doses, dose_matrix):
        differences = np.abs(predicted_doses - dose_matrix)
        percent_error = np.sum(differences, axis = 1)/np.sum(dose_matrix, axis = 1)
        return percent_error
    
    def get_input_features(self, data):
        num_patients = data.get_num_patients()
        pca = lambda x: Metrics.pca(x, self.num_pca_components)
        distances = pca(data.tumor_distances)
        lymph_nodes = pca(data.lymph_nodes)
        tumor_volumes = np.zeros((num_patients, 2))
        for i in range(num_patients):
            gtvs = data.gtvs[i]
            gtvp_volume = gtvs[0].volume
            gtvn_volume = 0
            for gtvn in gtvs[1:]:
                gtvn_volume += gtvn.volume
            tumor_volumes[i, :] = (gtvp_volume, gtvn_volume)
        laterality = data.lateralities.reshape(num_patients, 1)
        laterality = np.vectorize(TreeEstimator.laterality_map.__getitem__)(laterality)
        subsites = data.subsites.reshape(num_patients, 1)
        subsites = np.vectorize(TreeEstimator.subsite_map.__getitem__)(subsites)
        total_doses = data.prescribed_doses.reshape(num_patients, 1)
        clusters = data.classes.reshape(num_patients, 1)
        features = np.hstack([distances, lymph_nodes, tumor_volumes, total_doses, subsites, laterality])
        return features
            
class MLPEstimator(TreeEstimator):
    
    def __init__(self, num_pca_components = 10, n_estimators =10, min_samples_split = 4, max_depth = None):
        from sklearn.neural_network import MLPRegressor
        self.model = MLPRegressor(hidden_layer_sizes=(100,100,100), solver = 'lbfgs', alpha = .01)
        self.num_pca_components = num_pca_components
        
class KnnTreeEstimator(KnnEstimator):
    
    def __init__(self, num_pca_components = 10, n_estimators = 50, min_samples_split = 4, max_depth = None):
        from sklearn.ensemble import RandomForestRegressor
        self.model = RandomForestRegressor(n_estimators = n_estimators)
        self.num_pca_components = num_pca_components
        
    
    def get_prediction(self, data, scores, args, p):
        dose_matrix = data.doses
#        matched_scores = scores[args].reshape(len(args), 1
        matched_doses = dose_matrix[args, :]
        features = self.get_features(data)
        self.model.fit(features[args], matched_doses)
        return self.model.predict(features[p, :].reshape(1, features.shape[1]))
        
    def get_features(self, data):
        num_patients = data.get_num_patients()
        pca = lambda x: Metrics.pca(x, self.num_pca_components)
        lymph_nodes = pca(data.lymph_nodes)
        distances = pca(data.tumor_distances)
        tumor_volumes = np.zeros((num_patients, 2))
        for i in range(num_patients):
            gtvs = data.gtvs[i]
            gtvp_volume = gtvs[0].volume
            gtvn_volume = 0
            for gtvn in gtvs[1:]:
                gtvn_volume += gtvn.volume
            tumor_volumes[i, :] = (gtvp_volume, gtvn_volume)
        laterality = data.lateralities.reshape(num_patients, 1)
        laterality = np.vectorize(TreeEstimator.laterality_map.__getitem__)(laterality)
        subsites = data.subsites.reshape(num_patients, 1)
        subsites = np.vectorize(TreeEstimator.subsite_map.__getitem__)(subsites)
        total_doses = data.prescribed_doses.reshape(num_patients, 1)
        features = np.hstack([distances, lymph_nodes, tumor_volumes, total_doses, laterality, subsites])
        return features
    
    def get_num_matches(self, p, similarity, clusters):
        #for later better use probs
        num_cluster_values = len(np.where(clusters == clusters[p])[0])
        num_matches = np.max(np.sqrt([2*num_cluster_values, 3]))
        return int(num_matches)
        
class TreeSimilarity():
    
    subsite_map = {'BOT': 0, 'GPS': 1, 'NOS': 2, 'Soft palate': 3, 'Tonsil': 4}
    laterality_map = {'Bilateral': 0, 'L': 1, 'R': 2}
    
    def __init__(self, num_pca_components = 10, n_estimators =10, min_samples_split = 4, max_depth = None):
        from sklearn.ensemble import RandomForestRegressor
        self.model = RandomForestRegressor(min_samples_split = min_samples_split, 
                                           n_estimators = n_estimators,
                                           max_depth = max_depth)
        self.num_pca_components = num_pca_components
        
    def get_similarity(self, data):
        true_similarity = self.get_true_similarity(data)
        similarity = np.zeros(true_similarity.shape)
        x = self.get_input_features(data)
        for p in range(data.get_num_patients()):
            y_train = np.delete(true_similarity, p, axis = 0)
            x_train = np.delete(x, p, axis = 0)
            self.model.fit(x_train, y_train)
            similarity[p, :] = self.model.predict(x[p])
            similarity[p, p] = 0
        return similarity
            
    
    def get_true_similarity(self, data):
        n_patients = data.get_num_patients()
        doses = data.doses
        error_matrix = np.zeros((n_patients, n_patients))
        for p1 in range(n_patients):
            for p2 in range(p1 + 1, n_patients):
                dose_difference = np.abs(doses[p1,:] - doses[p2, :])
                error_matrix[p1, p2] = np.mean(dose_difference)
        similarity_matrix = 1 - (error_matrix - error_matrix.max())/(error_matrix.max() - error_matrix.min())
        similarity_matrix += similarity_matrix.transpose()
        return similarity_matrix
    
    def get_input_features(self, data):
        num_patients = data.get_num_patients()
        pca = lambda x: Metrics.pca(x, self.num_pca_components)
        distances = pca(data.tumor_distances)
        lymph_nodes = pca(data.lymph_nodes)
        tumor_volumes = np.zeros((num_patients, 2))
        for i in range(num_patients):
            gtvs = data.gtvs[i]
            gtvp_volume = gtvs[0].volume
            gtvn_volume = 0
            for gtvn in gtvs[1:]:
                gtvn_volume += gtvn.volume
            tumor_volumes[i, :] = (gtvp_volume, gtvn_volume)
        laterality = data.lateralities.reshape(num_patients, 1)
        laterality = np.vectorize(TreeSimilarity.laterality_map.__getitem__)(laterality)
        subsites = data.subsites.reshape(num_patients, 1)
        subsites = np.vectorize(TreeSimilarity.subsite_map.__getitem__)(subsites)
        total_doses = data.prescribed_doses.reshape(num_patients, 1)
        clusters = data.classes.reshape(num_patients, 1)
        features = np.hstack([distances, lymph_nodes, tumor_volumes, total_doses, laterality, subsites, clusters])
        return features
    
    
class DiscreteClassifierSimilarity(ClassifierSimilarity):
    
    def __init__(self, model = None, num_pca_components = 6, max_error = 3, min_matches = 3):
        from sklearn.multiclass import OneVsRestClassifier
        if model is None:
            from sklearn.naive_bayes import MultinomialNB
            model = MultinomialNB()
        self.model = OneVsRestClassifier(model)
        self.num_pca_components = num_pca_components
        self.max_error = max_error
        self.min_matches = min_matches
        
    def get_similarity(self, data, similarity = None):
        features = self.get_input_features(data)
        true_matches = self.get_true_matches(data)
        predicted_matches = np.zeros(true_matches.shape)
        for p in range(true_matches.shape[0]):
            train_features = np.delete(features, p, axis = 0)
            train_matches = np.delete(true_matches, p, axis = 0)
            predict_features = features[p,:].reshape(1, len(features[p,:]))
            self.model.fit(train_features , train_matches)
            predicted_matches[p, :] = self.model.predict_proba(predict_features)
            predicted_matches[p, p] = 0
            print(len(np.where(predicted_matches[p,:] > .4)[0]))
        return predicted_matches
    
    def get_input_features(self, data):
        num_patients = data.get_num_patients()
        pca = lambda x: Metrics.pca(x, self.num_pca_components)
        distances = np.round( pca(data.tumor_distances), 1)
        lymph_nodes = np.round( pca(data.lymph_nodes), 1)
        tumor_volumes = np.zeros((num_patients, 2))
        for i in range(num_patients):
            gtvs = data.gtvs[i]
            gtvp_volume = np.round( gtvs[0].volume, 1)
            gtvn_volume = 0
            for gtvn in gtvs[1:]:
                gtvn_volume += np.round( gtvn.volume, 1)
            tumor_volumes[i, :] = (gtvp_volume, gtvn_volume)
        laterality = data.lateralities.reshape(num_patients, 1)
        subsites = data.subsites.reshape(num_patients, 1)
        total_doses = data.prescribed_doses.reshape(num_patients, 1)
        clusters = data.classes.reshape(num_patients, 1)
        features = np.hstack([distances, lymph_nodes, tumor_volumes, total_doses, laterality, subsites])
        return features

class ClassifierSimilarity():
    
    def __init__(self, model = None, num_pca_components = 6, max_error = 3, min_matches = 3):
        from sklearn.multiclass import OneVsRestClassifier
        if model is None:
            from sklearn.linear_model import LogisticRegression
            model = LogisticRegression()
        self.model = OneVsRestClassifier(model)
        self.num_pca_components = num_pca_components
        self.max_error = max_error
        self.min_matches = min_matches
        
    def get_similarity(self, data, similarity = None):
        features = self.get_input_features(data)
        features = (features - features.mean(axis=0))/features.std(axis=0)
        features[:, 0] = 0
        true_matches = self.get_true_matches(data)
        predicted_matches = np.zeros(true_matches.shape)
        tsim_scores = TsimModel().get_similarity(data)
        for p in range(true_matches.shape[0]):
            train_features = np.delete(features, p, axis = 0)
            train_matches = np.delete(true_matches, p, axis = 0)
            predict_features = features[p,:]
            for p2 in range(true_matches.shape[1]):
                if p == p2:
                    continue
                train_features[: , 0] = np.delete(tsim_scores[:, p2], p, axis = 0)
                predict_features[0] = tsim_scores[p,p2]
                y = train_matches[:, p2].reshape(-1,1)
                print(y.shape, ' ', train_features.shape)
                self.model.fit(train_features , y)
                print(self.model.predict_proba(predict_features).shape)
                predicted_matches[p, p2] = self.model.predict_proba(predict_features)
#            predicted_matches[p, p] = 0
            print(len(np.where(predicted_matches[p,:] > .4)[0]))
        return predicted_matches
        
    def get_true_matches(self, data):
        dose_error = self.get_match_error(data)
        match_matrix = np.zeros(dose_error.shape)
        n_patients = data.get_num_patients()
        for p in range(n_patients):
            errors = dose_error[p, :]
            matches = []
            max_error = self.max_error
            while len(matches) < self.min_matches:
                matches = np.where(errors < max_error)[0]
                max_error = max_error + .2
            match_matrix[p, matches] = 1
        return match_matrix
        
    def get_match_error(self, data):
        n_patients = data.get_num_patients()
        doses = data.doses
        error_matrix = np.zeros((n_patients, n_patients))
        for p1 in range(n_patients):
            for p2 in range(p1 + 1, n_patients):
                dose_difference = np.abs(doses[p1,:] - doses[p2, :])
                error_matrix[p1, p2] = np.mean(dose_difference)
        error_matrix += error_matrix.transpose()
        return error_matrix
    
    def get_input_features(self, data):
        num_patients = data.get_num_patients()
        pca = lambda x: Metrics.pca(x, self.num_pca_components)
        distances = pca(data.tumor_distances)
        lymph_nodes = pca(data.lymph_nodes)
        tumor_volumes = np.zeros((num_patients, 2))
        for i in range(num_patients):
            gtvs = data.gtvs[i]
            gtvp_volume = gtvs[0].volume
            gtvn_volume = 0
            for gtvn in gtvs[1:]:
                gtvn_volume += gtvn.volume
            tumor_volumes[i, :] = (gtvp_volume, gtvn_volume)
        laterality = data.lateralities.reshape(num_patients, 1)
        laterality = np.vectorize(TreeSimilarity.laterality_map.__getitem__)(laterality)
        subsites = data.subsites.reshape(num_patients, 1)
        subsites = np.vectorize(TreeSimilarity.subsite_map.__getitem__)(subsites)
        total_doses = data.prescribed_doses.reshape(num_patients, 1)
        clusters = data.classes.reshape(num_patients, 1)
        features = np.hstack([distances, lymph_nodes, tumor_volumes, total_doses, subsites])
        return features

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

##old code from testing neural nets
        def get_autoencoderish_model(features):
    input_x = Input(shape=(features.shape[1],))
    encoder = Sequential([
            Dense(45, input_dim=features.shape[1], activation = 'relu'),
            Dense(100, activation = 'relu'),
            Dense(100, activation = 'relu'),
            Dense(100, activation = 'relu'),
            ])(input_x)

    decoder = Sequential([
            Dense(100,input_dim = 4, activation = 'relu',
                  activity_regularizer = regularizers.l2(.01)),
            Dense(45, activation = 'relu'),
            ])(encoder)
    model = Model(input_x, decoder)
    encoder_model= Model(input_x, encoder)
#    optimizer = optimizers.SGD(lr = .01, decay = 1e-12, momentum = .1)
    optimizer = optimizers.Adam()
    model.compile(loss = losses.mean_absolute_error,
                  optimizer = optimizer)
    return(model, encoder_model)

def get_regression_model(features, activation = 'relu', lr = .01):
    model = Sequential([
            Dense(45, input_dim=features.shape[1], activation = activation),
            Dense(100, activation = activation),
            Dense(200, activation = activation),
            Dense(45, activation = activation)
            ])
    optimizer = optimizers.SGD(lr = lr, decay = 1e-4, momentum = 0.05)
    model.compile(loss = losses.mean_absolute_error,
                  optimizer = optimizer)
    return(model)

def run_autoencoder(db):
    from keras.models import Sequential, Model
    from keras.layers import Dense, Activation, Input
    from keras import losses, optimizers,regularizers
    from sklearn.model_selection import LeaveOneOut
    features = get_input_distance_features(db)
    features = (features - features.mean(axis = 0))/features.std(axis = 0)
    clusters = db.classes.astype('int32')
    doses = db.doses

    loo = LeaveOneOut()
    loo.get_n_splits(features)
    regression_errors = []
    nn_sim = np.zeros((db.get_num_patients(),db.get_num_patients()))
    p1 = 0
    for train,test in loo.split(features, doses):
        model, encoder_model = get_autoencoderish_model(features)
        x_train = features[train]
        y_train = doses[train]
        x_test = features[test]
        y_test = doses[test]
        model.fit(x_train, y_train, epochs = 3, batch_size = 4, shuffle = True, verbose = 0)
        regression_error = model.evaluate(x_test, y_test)
        regression_errors.append(regression_error)
        print(regression_error)

        x_embedding = encoder_model.predict(features)
        for p2 in range(db.get_num_patients()):
            if p1 == p2:
                continue
            nn_sim[p1, p2] = 1/np.linalg.norm(x_embedding[p1] - x_embedding[p2])
        p1 += 1
    print(np.mean(regression_errors))
    nn_sim = (nn_sim - nn_sim.min(axis = 0))/(nn_sim.max(axis = 0) - nn_sim.min(axis = 0))

    threshold_grid_search(db, nn_sim)

def get_features(db, holdout = set([]) ):
    n_patients = db.get_num_patients()
    x = get_input_distance_features(db)
    normalizer = Normalizer()
    normalizer.fit( x[[ i for i in np.arange(n_patients) if i not in holdout] ] )
    x = normalizer.transform(x)
    a = []
    b = []
    y = []
    a_validate = []
    b_validate = []
    y_validate = []
    val_pairs = []
    true_error = SimilarityFuser(min_matches = 7, max_error = .20).get_true_matches(db)
    for p1 in range(n_patients):
        for p2 in range(p1 + 1, n_patients):
            loss = 0 if true_error[p1,p2] > 0 else 1
            if p1 in holdout or p2 in holdout:
                a_validate.append(x[p1])
                b_validate.append(x[p2])
                y_validate.append( loss )
                val_pairs.append((p1,p2))
            else:
                a.append(x[p1])
                b.append(x[p2])
                y.append( loss )
    a = np.array(a)
    b = np.array(b)
    y = np.array(y).ravel()
    a_validate = np.array(a_validate)
    b_validate = np.array(b_validate)
    y_validate = np.array(y_validate)

    args = np.arange(len(y))
    np.random.shuffle(args)

    return (a[args], b[args], y[args], a_validate, b_validate, y_validate, val_pairs)

def get_similarity_model(n_features, encoding_size = 25, reg = .000001):
    patient_a = layers.Input(shape = (n_features,))
    patient_b = layers.Input(shape = (n_features,))
    activation = 'selu'
    encoder = Sequential([
                Dense(50, input_dim = n_features, activation = activation,
                      activity_regularizer = regularizers.l2( reg )),
                Dense(100, activation = activation,
                      activity_regularizer = regularizers.l2( reg )),
                layers.Dropout(.5, seed = 0),
                Dense(encoding_size, activation = 'relu'),
                ])
    encoded_a = encoder(patient_a)
    encoded_b = encoder(patient_b)
    distance_layer = layers.dot([encoded_a, encoded_b], axes = 1, normalize = True)
    model = Model([patient_a, patient_b], distance_layer)
#    optimizer = optimizers.SGD(lr = .001, decay = 1e-8, momentum = .01)
    optimizer = optimizers.Adam()
    model.compile(optimizer = optimizer, loss = losses.mean_absolute_error)
    return(model)

def get_distance_model(n_features, encoding_size = 25, reg = 0.000001):
    patient_a = layers.Input(shape = (n_features,))
    patient_b = layers.Input(shape = (n_features,))
    activation = 'relu'
    encoder = Sequential([
                Dense(500, input_dim = n_features, activation = activation,
                      activity_regularizer = regularizers.l2( reg )),
                layers.Dropout(.1, seed = 0),
                Dense(encoding_size, activation = 'relu'),
                layers.BatchNormalization(),
                ])
    encoded_a = encoder(patient_a)
    encoded_b = encoder(patient_b)
    distance_layer = layers.Lambda(lambda x: K.expand_dims(K.mean(K.square(x[0] - x[1]),axis=-1),1),
                                   output_shape=(1,))([encoded_a, encoded_b])
    distance_activation = Activation('sigmoid')(distance_layer)
    model = Model(inputs=[patient_a, patient_b], outputs = distance_activation)
    distance_model = Model(inputs=[patient_a, patient_b], outputs = distance_layer)
#    optimizer = optimizers.SGD(lr = .01, decay = 1e-8, momentum = .01, nesterov = True)
    optimizer = optimizers.Adam(lr = .001, decay = 1e-8)
#    optimizer = optimizers.Adadelta()
#    optimizer = optimizers.RMSprop()
    model.compile(optimizer = optimizer,
                  loss = losses.mean_squared_error)
    return(model, distance_model)

### stuff with meetric learning?
def get_all_features(data, num_pca_components = 10):
    num_patients = data.get_num_patients()
    tumor_volumes = np.zeros((num_patients, 2))
    tumor_count = np.array([len(gtv) for gtv in data.gtvs]).reshape(-1,1)
    for i in range(num_patients):
        gtvs = data.gtvs[i]
        gtvp_volume = gtvs[0].volume
        gtvn_volume = 0
        for gtvn in gtvs[1:]:
            gtvn_volume += gtvn.volume
        tumor_volumes[i, :] = (gtvp_volume, gtvn_volume)
    laterality = data.lateralities.reshape(num_patients, 1)
    laterality = np.vectorize(Constants.laterality_map.__getitem__)(laterality)
    subsites = data.subsites.reshape(num_patients, 1)
    subsites = np.vectorize(Constants.subsite_map.__getitem__)(subsites)
    total_doses = data.prescribed_doses.reshape(num_patients, 1)
    ages = data.ages.reshape(-1,1)
    features = np.hstack([tumor_volumes, total_doses, subsites,
                          laterality, tumor_count,
                          data.tumor_distances, data.volumes])
    pca_features = pca(features, features.shape[1])
    return (pca_features - pca_features.mean(axis = 0))/pca_features.std(axis = 0)

def get_input_tumor_features(data, num_pca_components = 10):
    num_patients = data.get_num_patients()
    tumor_volumes = np.zeros((num_patients, 2))
    tumor_count = np.array([len(gtv) for gtv in data.gtvs]).reshape(-1,1)
    for i in range(num_patients):
        gtvs = data.gtvs[i]
        gtvp_volume = gtvs[0].volume
        gtvn_volume = 0
        for gtvn in gtvs[1:]:
            gtvn_volume += gtvn.volume
        tumor_volumes[i, :] = (gtvp_volume, gtvn_volume)
    laterality = data.lateralities.reshape(num_patients, 1)
    laterality = np.vectorize(Constants.laterality_map.__getitem__)(laterality)
    subsites = data.subsites.reshape(num_patients, 1)
    subsites = np.vectorize(Constants.subsite_map.__getitem__)(subsites)
    total_doses = data.prescribed_doses.reshape(num_patients, 1)
    features = np.hstack([tumor_volumes, total_doses,
                          subsites, tumor_count,
                          data.ajcc8.reshape(-1,1)])
    return copy.copy(features)

def get_input_organ_features(data):
    return copy.copy(data.volumes)

def get_input_distance_features(data, num_pca_components = 10):
    num_patients = data.get_num_patients()
    tumor_volumes = np.zeros((num_patients, 2))
    tumor_count = np.array([len(gtv) for gtv in data.gtvs]).reshape(-1,1)
    for i in range(num_patients):
        gtvs = data.gtvs[i]
        gtvp_volume = gtvs[0].volume
        gtvn_volume = 0
        for gtvn in gtvs[1:]:
            gtvn_volume += gtvn.volume
        tumor_volumes[i, :] = (gtvp_volume, gtvn_volume)
    laterality = data.lateralities.reshape(num_patients, 1)
    laterality = np.vectorize(Constants.laterality_map.__getitem__)(laterality)
    subsites = data.subsites.reshape(num_patients, 1)
    subsites = np.vectorize(Constants.subsite_map.__getitem__)(subsites)
    total_doses = data.prescribed_doses.reshape(num_patients, 1)
    features = np.hstack([tumor_volumes,
                          total_doses,
                          tumor_count,
                          laterality,
                          data.ajcc8.reshape(-1,1),
                          data.tumor_distances])
    return copy.copy(features)

def get_input_lymph_features(data, num_pca_components = 10):
    return copy.copy(data.lymph_nodes)

def get_dose_clusters(doses):
    kmeans = KMeans(n_clusters = 5, random_state = 0)
    bad_patients = set(ErrorChecker().get_data_outliers(doses))
    good_patients = [i for i in np.arange(doses.shape[0]) if i not in bad_patients]
    kmeans_clusters = kmeans.fit_predict(doses[good_patients])
    return kmeans_clusters, good_patients

def get_nca_features(features, doses, min_components = 5, lmnn = False, k = 4, reg = .25):
    n_patients = doses.shape[0]
    if lmnn is False:
        n_components = min([min_components, features.shape[1]])
    else:
        n_components = features.shape[1]
    output = np.zeros((n_patients,n_components))
    for p in range(n_patients):
        feature_subset = np.delete(features, p, axis = 0)
        dose_subset = np.delete(doses, p , axis = 0)
        nca = get_fitted_nca(feature_subset, dose_subset,
                                          n_components = n_components,
                                          lmnn = lmnn, k = k,
                                          reg = reg)
        tf = nca.transform(features[p,:].reshape(1,-1))
#        output[p,:] = (tf - tf.min())/(tf.max() - tf.min())
        output[p,:] = tf/np.linalg.norm(tf)
        print(output[p,:])
    return output

def get_fitted_nca(features, doses, n_components = 5, lmnn = False, k = 4, reg = .25):
    if lmnn:
        nca = metric_learn.lmnn.python_LMNN(k = k, use_pca = False,
                                            regularization = reg)
    else:
        nca = NeighborhoodComponentsAnalysis(n_components = n_components,
                                         max_iter = 300,
                                         init = 'pca',
                                         random_state = 0)
    kmeans_clusters, good_clusters = get_dose_clusters(doses)
    for col in range(features.shape[1]):
        feature = features[:, col]
        if feature.std() > .0001:
            features[:, col] = (feature - feature.mean())/feature.std()
        else:
            features[:, col] = feature - feature.mean()
    sampler = SMOTE(k_neighbors  = 2)
    resampled_features, resampled_clusters = sampler.fit_resample(
            features[good_clusters, :], kmeans_clusters)
    nca.fit(resampled_features, resampled_clusters)
#    features = nca.transform(features)
    return nca

def get_nca_similarity(db, feature_type = 'tumors', min_nca_components = 4, lmnn = False, k = 4, reg = .25):
    doses = db.doses
    n_patients = doses.shape[0]
    if feature_type == 'tumors':
        input_features = get_input_tumor_features(db)
    elif feature_type in ['distance', 'distances']:
        input_features = get_input_distance_features(db)
    elif feature_type in ['lymph', 'lymph nodes']:
        input_features = get_input_lymph_features(db)
    elif feature_type in ['organ', 'organs']:
        input_features = get_input_organ_features(db)
    nca_features = get_nca_features(input_features, doses,
                                    min_components = min_nca_components,
                                    lmnn = lmnn, k = k, reg = reg)
    similarity = np.zeros((n_patients, n_patients))
    max_similarities = set([])
    mixed_laterality = set(['R','L'])
    for p1 in range(n_patients):
        x1 = nca_features[p1, :]
        for p2 in range(p1+1, n_patients):
            if set(db.lateralities[[p1,p2]]) == mixed_laterality:
                continue
            x2 = nca_features[p2, :]
            if np.linalg.norm(x1 - x2) < .001:
                max_similarities.add((p1,p2))
                continue
            similarity[p1, p2] = 1/np.linalg.norm(x1 - x2)
    similarity += similarity.transpose()
    similarity = .99*(similarity - similarity.min())/(similarity.max() - similarity.min())
    for pair in max_similarities:
        similarity[pair[0], pair[1]] = 1
    for i in range(n_patients):
        similarity[i,i] = 0
    return similarity

def get_test_tumor_similarity(db):
    n_patients = db.get_num_patients()
    new_distance_similarity = np.zeros((n_patients, n_patients))
    for i in range(n_patients):
        gtvs1 = db.gtvs[i]
        for ii in range(n_patients):
            gtvs2 = db.gtvs[ii]
            new_distance_similarity[i,ii] = get_max_tumor_ssim(gtvs1, gtvs2)
    new_distance_similarity += new_distance_similarity.transpose()
    return new_distance_similarity

def centroid_based_tumor_organ_pairs(db):
    for p in range(db.get_num_patients()):
        p_centroids = db.centroids[p,:,:]
        gtv = db.gtvs[p]
        for t_ind in range(len(gtv)):
            t = gtv[t_ind]
            organ = Constants.organ_list[0]
            min_dist = np.linalg.norm(t.position - p_centroids[0])
            for o in range(1, Constants.num_organs):
                dist = np.linalg.norm(t.position - p_centroids[o])
                if dist < min_dist:
                    min_dist = dist
                    organ = Constants.organ_list[o]
            gtv[t_ind] = GTV(t.name, t.volume, t.position, t.doses, t.dists, organ)

def organ_selection(organ_list, db, similarity_function = None,
                    use_classes = False):
    def tsim(x):
        model = TsimModel(organs = [Constants.organ_list.index(o) for o in x],
                                   similarity_function = similarity_function,
                                   use_classes = use_classes)
        return model.get_similarity(db)
    distance_similarity = tsim(organ_list)
    baseline = threshold_grid_search(db, distance_similarity)[0]
    optimal = (organ_list, baseline)
    bad_organs = []
    best_score = 100
    for organ in organ_list:
        organ_subset = copy(organ_list)
        organ_subset.remove(organ)
        distance_subset_sim = tsim(organ_subset)
        best_score, best_threshold, best_min_matches = threshold_grid_search(db, distance_subset_sim, print_out = False)
        if best_score < baseline:
            bad_organs.append((organ, best_score, best_threshold, best_min_matches))
            if best_score < optimal[1]:
                optimal = (organ_subset, best_score)
                print(set(Constants.organ_list) - set(optimal[0]), best_score)
    return optimal

def optimal_organ_search(db, similarity_function = None, use_classes = False):
    optimal_organs = []
    organ_set = Constants.organ_list
    best_score = None
    while True:
        optimal_organs, best = organ_selection(organ_set, db,
                                               similarity_function = similarity_function,
                                               use_classes = use_classes)
        if len(optimal_organs) == len(organ_set):
            break
        organ_set = optimal_organs
        best_score = best
    return optimal_organs, best_score

def get_tumor_organ_vectors(db):
    o_centroids, t_centroids = db.get_transformed_centroids()
    vectors = np.zeros((o_centroids.shape))
    for p in range(db.get_num_patients()):
        o_centers = o_centroids[p]
        t_centers = t_centroids[p]
        distances = np.stack([g.dists for g in db.gtvs[p] if g.volume > 0], axis = 1)
        new_vectors = np.zeros(o_centers.shape)
        for organ in range(new_vectors.shape[0]):
            nearest_tumor_arg = np.argmin(distances[organ])
            nearest_tumor = t_centers[nearest_tumor_arg]
            organ_tumor_vector = nearest_tumor - o_centers[organ]
            new_vectors[organ] = organ_tumor_vector/np.linalg.norm(organ_tumor_vector)
        vectors[p] = new_vectors
    return vectors

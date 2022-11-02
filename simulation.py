import numpy as np
import pandas as pd
import random 
import igraph as ig
from random import seed
from random import randrange
#import matplotlib.pyplot as plt

def matching_simulation(n_feature, n_meta, N, ba_m, ba_power, multi_rate, poisson_rate, max_meta_link, max_feature_link, effective_ratio, link, cov_base):
    """
    The function simulation the uncertainty matching problem between known features and unknown metabolites
    
    INPUT: 
    n_feature/n_meta/N: the number of features/meta/sample
    ba_m/ba_power: the parameters in igraph to simulate the network of metabolites
    multi_rate: the ratio of features having multiple matching to metabolites
    poisson_rate: the parameter of Poisson distribution when determing the number of multiple matching
    max_meta_link/max_feature_link: the cap of number of matchings connected to one metabolite/feature
    efffective_ratio: the ratio of effective metabolites which have true matchings to the features and determine the output Y
    link: the link function from X to Y, can be as "logistic" or "abs"
    cov_base: the parameters for generating covariance matrix for X
    
    OUTPUT: 
    a list containing the following info will be returned:
    graph: the metabolites graph which is using igraph package
    X: the simulated features expression matrix
    dist: the shorest distance matrix for metabolites graph
    feature_meta: the matching between features and metabolites
    feature_meta_true: the true matching between features and metabolites
    pool: the list of effective metabolites
    """
    def simulation_ba(ba_m, num_nodes, ba_power):
        graph = ig.Graph.Barabasi(n=num_nodes, m = ba_m, power=ba_power, implementation='psumtree',outpref=False,
                                  directed = False, zero_appeal= 0)
        dg = np.array(graph.degree())
        dist = np.array(graph.shortest_paths())
        return(graph, dg, dist)

    graph, dg, dist = simulation_ba(ba_m, n_meta, ba_power)
    weight_data = dict(zip(range(n_meta), dg))

    feature_meta = np.zeros([n_feature, n_meta])
    feature_meta_true = np.zeros([n_feature, n_meta])

    connected_feature = []
    connected_meta = []
    pool = np.array(range(n_meta)).tolist()
    
    for i in range(n_feature):
        #if multi[i] == 1:
        #    idx1 = np.random.choice(pool,poisson[i],replace = False).tolist()
        #else:
        idx1 = np.random.choice(pool,1)[0]
        feature_meta_true[i, idx1] = 1
        idx1_neighbor = np.where((dist[idx1,:]<=3) & (dist[idx1,:]>=1))[0]  # the false connection can only connect to the nearby meta nodes
        meta_link = np.sum(feature_meta,0)

        idx1_neighbor = idx1_neighbor[np.where(meta_link[idx1_neighbor]<max_meta_link)[0]]
        connected_meta = connected_meta + [idx1]

        t = ((np.random.poisson(poisson_rate, 1) + 1) * (multi_rate > np.random.random_sample(1)))[0]
        if t!=0:
            idx2 = np.random.choice(idx1_neighbor, min(t,len(idx1_neighbor)))

            feature_meta[i, [idx1] + list(idx2)] = 1
        else:
            feature_meta[i, idx1] = 1
        pool.remove(idx1)

    connected_feature = np.unique(connected_feature).tolist()

    connected_meta_count = feature_meta_true.sum(axis = 0) # the true connection
    meta_conn = np.where(connected_meta_count==1)[0]
    dg_conn = [weight_data[i] for i in meta_conn]


    cores = meta_conn[np.where(dg_conn > np.quantile(dg_conn, 0.99))[0]]  # 0.98
    cores = np.random.choice(cores, 1) # only choose one as the core node
    count = len(cores)
    pool = list(cores)
    cc = 1
    while len(pool)< n_feature * effective_ratio and cc < 90:
        u = []
        for i in cores:
            i_neighbor = (np.where(dist[i,:]==1)[0]).tolist()
            i_neighbor = [i for i in i_neighbor if i in meta_conn] # insure the effective meta is true meta
            random.shuffle(i_neighbor)
            add_neighbor = i_neighbor[0:int(len(i_neighbor) * 0.5)]
            if len(add_neighbor)!=0:
                pool = pool + add_neighbor
                u = u + add_neighbor
        u = np.unique(u).tolist()
        cores = u
        pool = (np.unique(pool)).tolist()
        #print("Fining in length",str(cc),"and have choose",len(pool),"nodes")
        cc = cc + 1
    pool = pool[0:int(n_feature * effective_ratio)]
   
    # generate X matrix
    # use the covariance of the true connection in meta nodes to build up the covariance for features
    feature_meta_true2 = np.zeros([n_feature, n_meta])
    feature_meta_true2[:,pool] = feature_meta_true[:,pool]
    true_features = np.where(feature_meta_true[:,pool] ==1)[0]
    covariance =  np.dot(np.dot(feature_meta_true2, np.power(cov_base, dist)),feature_meta_true2.T) 
    mean = np.zeros([n_feature])-1
    mean[np.random.choice(np.arange(n_feature), int(len(mean)/2))] = 1
    
    X = np.random.multivariate_normal(mean, covariance, N)

    low = 0.7
    high = 0.9
    beta_negative_rate = 0.5
    epsilon = True
    weight_a = 1
    beta = np.array([random.uniform(low, high) for i in range(len(true_features)+1)])  # add one constant beta_0
    negative_index = np.random.randint(len(true_features)+1, size=int(len(true_features)+1*beta_negative_rate))
    beta[negative_index] = beta[negative_index] * (-1)

    Z = np.dot(np.column_stack((np.ones((N,1)),X[:,true_features])), beta)
    
    y = 1 / (weight_a * (1 + np.exp(- Z + np.median(Z))))
    y = (y >= 0.5/weight_a) + 0
    
    print("Finish generating Y with 1 class %d, 0 class %d" %((y==1).sum(), (y==0).sum()))
    X = np.column_stack((X,y))
    res = [graph, X, dist, feature_meta, feature_meta_true, pool]
    return(res)

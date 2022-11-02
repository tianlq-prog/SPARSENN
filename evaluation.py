import numpy as np
from sklearn.utils import shuffle
from sklearn import metrics
def evaluation(net, feature_meta_true, true_meta, degree, count_feature):
    """
    GOAL: evaluation of model performance from the aspect of feature selection
    
    INPUT VARIABLE:
    net: the NN model
    feature_meta_true: the true connection between features and metabolites
    true_meta: the list of predictive metabolites
    degree: the node degree of metabolic network
    count_feature: the number of matching connected to each feature
    
    OUTPUT VARIABLE:  rate, feature_auc, feature_prauc, meta_auc, meta_prauc, degree_auc, degree_prauc
    rate: the ratio of true predicted matching between features and metabolites
    feature_auc: the AUC between true predictive feature and our estimated feature importance
    feature_prauc: the PR-AUC between true predictive feature and our estimated feature importance
    meta_auc: the AUC between true predictive metabolites and our estimated metabolic importance
    meta_prauc: the PR-AUC between true predictive metabolites and our estimated metabolic importance    
    link_auc: the AUC between true predictive metabolites and the degree of metabolic network
    link_prauc: the PR-AUC between true predictive metabolites and the degree of metabolic network
    """
    params = []
    for parameters in net.parameters():
        params.append(parameters.detach().cpu().numpy())

    h0 = params[0].T
    h1 = params[2].T

    f_left = np.sum(abs(h0),0)/count_feature**0.5
    f_right = np.sum(abs(h1),1)/degree**0.5
    f_imp = (f_left/f_left.sum()) + (f_right/f_right.sum())

    n_meta = feature_meta_true.shape[1]
    n_feature = feature_meta_true.shape[0]
    true_meta = np.array(true_meta).astype(int)
    y_meta  = np.zeros([n_meta])
    y_meta[true_meta] = 1

    true_feature = np.where(feature_meta_true[:,true_meta]==1)[0] # index of true features

    true_feature_v = np.zeros([n_feature])
    true_feature_v[true_feature] = 1
    true_feature_v = true_feature_v.astype(int) # vector representation of true features

    # calculate the rate of true link
    idx1,idx2 = np.where(feature_meta_true == 1)
    l = np.zeros(feature_meta_true.shape[0])
    l[idx1] = idx2

    max_idx = np.argmax(abs(h0)/(np.sum(abs(h0),0)+ 0.000001) , axis=1)
    rate = ((l==max_idx).sum()-(l==0).sum())/((l!=0).sum())

    h0_1 = np.zeros(h0.shape)
    h0_1[h0!=0] = 1
    
    h0_new = np.zeros(h0.shape)
    h0_new[range(n_feature), max_idx] = h0[range(n_feature), max_idx]
    h0_new[h0_new!=0] = 1
    feature_imp2 = np.dot(h0_new, f_imp)
    
    # feature auc & prauc
    feature_auc = metrics.roc_auc_score(true_feature_v,feature_imp2)
    feature_prauc = metrics.average_precision_score(true_feature_v,feature_imp2)
    
    # meta auc & prauc
    meta_auc = metrics.roc_auc_score(y_meta,f_imp)
    meta_prauc = metrics.average_precision_score(y_meta,f_imp)
   
    # degree auc & prauc
    degree_auc = metrics.roc_auc_score(y_meta,degree)
    degree_prauc = metrics.average_precision_score(y_meta,degree)

    res = [rate, feature_auc, feature_prauc, meta_auc, meta_prauc, degree_auc, degree_prauc]
    return(res)
import numpy as np
import math
import torch
import torch.nn as nn
from torch.nn.parameter import Parameter
import torch.nn.functional as F
import torch.nn.init as init
import igraph

def getLayerSizeList(partition, threshold_layer_size, sparsify_coefficient):
    """
    Obtain the size of each sparse layer
    
    INPUT:
    partition: the adjacent matrix of metabolic network
    threshold_layer_size: the threshold of sparese layer
    sparsify_coefficient: the coefficient of each sparse level
    
    OUTPUT:
    sparsify_hidden_layer_size_dict: a dictionary indicating the sparse layer
    """
    n_meta = np.shape(partition)[0]
    n_layer = math.floor(np.log10(1.0 * threshold_layer_size / n_meta) / np.log10(sparsify_coefficient)) + 3
    # dict for number of neurons in each layer
    sparsify_hidden_layer_size_dict = {}

    sparsify_hidden_layer_size_dict['n_hidden_0'] = int(n_meta)

    for i in range(1,n_layer):
        sparsify_hidden_layer_size_dict['n_hidden_%d' % (i)] = int(n_meta * (sparsify_coefficient) ** (i-1))
    return sparsify_hidden_layer_size_dict


def getPartitionMatricesList(sparsify_hidden_layer_size_dict, degree_dict, feature_meta, partition):
    """
    Obtain the linkage matrix among two sparse layers
    """
    np.random.seed(1);  # for reproducable result
    g = igraph.Graph.Adjacency((partition).tolist(), mode = "undirected")
    dist = np.array(g.shortest_paths()) # use the shortest distance matrix to assign links
    
    sum_remove_node_list = []  # keep note of which nodes are already removed
    
    partition_mtx_dict = {}

    partition_mtx_dict["p0"] = feature_meta  # first matrix being the connection from features to meta
    partition_mtx_dict["p1"] = partition  # first matrix being the whole adjacency matrix

    for i in range(2, len(sparsify_hidden_layer_size_dict)):
        num_nodes_to_remove = sparsify_hidden_layer_size_dict["n_hidden_%d" % (i-1)] - \
                              sparsify_hidden_layer_size_dict["n_hidden_%d" % (i)]

        sorted_node_degree_list = sorted(degree_dict.items(), key=lambda item: item[
            1])  # sort node degree dict according to number of degrees

        temp_remove_list = []
        max_to_remove_node_degree = sorted_node_degree_list[num_nodes_to_remove - 1][1]
        for j in range(num_nodes_to_remove):  # any node with degree less than `max_to_remove_node_degree` is certain to be removed
            if sorted_node_degree_list[j][1] < max_to_remove_node_degree:
                id_to_remove_node = sorted_node_degree_list[j][0]
                temp_remove_list.append(id_to_remove_node)
            else:
                break  # node with more degrees is not under consideration

        # sample from all nodes that have max_to_remove_node_degree to reach number of nodes required to be removed
        sample_list = []
        for j in range(len(temp_remove_list), len(sorted_node_degree_list)):
            if sorted_node_degree_list[j][1] == max_to_remove_node_degree:
                sample_list.append(sorted_node_degree_list[j]);
            else:
                break  # node with more degrees is not under consideration

        sample_idx_list = sorted(
            np.random.choice(len(sample_list), num_nodes_to_remove - len(temp_remove_list), replace=False))
        for idx in sample_idx_list:
            temp_remove_list.append(sample_list[idx][0])

        # sum up add nodes to be removed
        all_list = np.arange(partition.shape[0])
        previous_layer_list = [x for x in all_list if x not in sum_remove_node_list]
        temp_partition = np.delete(partition, sum_remove_node_list, axis=0)
        sum_remove_node_list += temp_remove_list
        temp_partition = np.delete(temp_partition, sum_remove_node_list, axis=1)
        next_layer_list = [x for x in all_list if x not in sum_remove_node_list]
        # assign each neuron at least one linkage
        for k in range(len(previous_layer_list)):
            if sum(dist[k,next_layer_list]==float("inf"))==len(next_layer_list):
                idx = np.random.choice(len(next_layer_list), 1, replace=False)
            else:
                idx = np.argsort(dist[k,next_layer_list], axis = -1)[0]
            temp_partition[k, idx] = 1
        
        for j in range(len(temp_remove_list)):
            degree_dict.pop(temp_remove_list[j])

        partition_mtx_dict["p%d" % i] = temp_partition

    return partition_mtx_dict

def getNodeDegreeDict(partition):
    """
    Obtain the node degree using the adjacent matrix of metabolic network
    """
    degree_dict = {}
    row, col = partition.shape
    for i in range(row):
        degree_dict[i] = -1  # decrease its own
        for j in range(0, col):
            if partition[i, j] == 1:
                degree_dict[i] += 1

    return degree_dict

############### Sparse-nn function #################
def truncated_normal_(tensor,mean=0,std=0.09):
    with torch.no_grad():
        size = tensor.shape
        tmp = tensor.new_empty(size+(4,)).normal_()
        valid = (tmp < 2) & (tmp > -2)
        ind = valid.max(-1, keepdim=True)[1]
        tensor.data.copy_(tmp.gather(-1, ind).squeeze(-1))
        tensor.data.mul_(std).add_(mean)
        return tensor
    
class myLinear(nn.Module):
    __constants__ = ['bias']
 
    def __init__(self, in_features, out_features, bias=True):
        super(myLinear, self).__init__()
        self.in_features = in_features
        self.out_features = out_features
        self.weight = Parameter(torch.Tensor(out_features, in_features))
        if bias:
            self.bias = Parameter(torch.Tensor(out_features))
        else:
            self.register_parameter('bias', None)
        self.reset_parameters()
 
    def reset_parameters(self):
        self.weight = truncated_normal_(self.weight, mean = 0, std = 0.1)
        if self.bias is not None:
            fan_in, _ = init._calculate_fan_in_and_fan_out(self.weight)
            bound = 1 / math.sqrt(fan_in)
            init.uniform_(self.bias, -bound, bound)
 
    
    def forward(self, input):
        return F.linear(input, self.weight, self.bias)
    
class SparseLinear(nn.Module):
    """
    Define our linear connection layer which enabled sparse connection
    """
    def __init__(self, in_dim, out_dim, m):
        indices_mask = [np.where(m==1)[1].tolist(),np.where(m==1)[0].tolist()]
        
        super(SparseLinear, self).__init__()
 
        def backward_hook(grad):
            # Clone due to not being allowed to modify in-place gradients
            out = grad.clone()
            out[self.mask] = 0
            return out
 
        self.linear = myLinear(in_dim, out_dim, True)
        self.mask = torch.ones([out_dim, in_dim]).byte()
        self.mask[indices_mask] = 0 # create mask
        self.linear.weight.data[self.mask] = 0 # zero out bad weights
        self.linear.weight.register_hook(backward_hook) # hook to zero out bad gradients
 
    def forward(self, input):
        return self.linear(input)

class Net(nn.Module):
    """
    The network structure
    """
    def __init__(self,partition_mtx_dict, num_hidden_layer_neuron_list, keep_prob):
        super(Net,self).__init__()
        layer1 = nn.Sequential()
        for i in range(len(partition_mtx_dict)):
            mtx = partition_mtx_dict["p%d" % i]  # the mask matrix
            layer1.add_module('f'+ str(i), SparseLinear(mtx.shape[0], mtx.shape[1], mtx))
            layer1.add_module("f_relu"+str(i), nn.ReLU(True))
            #layer1.add_module("bn1"+str(i), nn.BatchNorm1d(mtx.shape[1]))
        self.layer1 = layer1
        
        layer2 = nn.Sequential()
        num_hidden_layer_neuron_list = [mtx.shape[1]] + num_hidden_layer_neuron_list + [2]
        for j in range(1, len(num_hidden_layer_neuron_list)-1):
            layer2.add_module('h'+str(j), myLinear(num_hidden_layer_neuron_list[j-1], num_hidden_layer_neuron_list[j]))
            layer2.add_module('h_relu'+str(j), nn.ReLU(True))
            #layer2.add_module("bn2"+str(j), nn.BatchNorm1d(num_hidden_layer_neuron_list[j]))
            layer2.add_module('h_drop'+str(j), nn.Dropout(p=keep_prob))
        layer2.add_module('h'+str(j+1), myLinear(num_hidden_layer_neuron_list[j], num_hidden_layer_neuron_list[j+1]))
        self.layer2 = layer2
        
    def forward(self,input):
        out = self.layer1(input)
        out = self.layer2(out)
        return out
        
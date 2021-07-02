import os
import sys
import numpy as np

import torch

import time


def forloopdists(feats, protos):

    dist = np.empty((feats.shape[0], protos.shape[0]))
    for i in range(feats.shape[0]):
        for k in range(protos.shape[0]):
        dist[i, k] = np.dot(feats[i, :]-protos[k, :],
                            (feats[i, :]-protos[k, :]).transpose())

    return dist


def numpydists(feats, protos):
    dist = -2 * np.matmul(feats, protos.transpose())
    dist += np.sum(feats**2, axis=1)[:, np.newaxis]
    dist += np.sum(protos**2, axis=1)[np.newaxis, :]
    return dist

# def pytorchdists(feats0,protos0,device):

#   ft=torch.from_numpy(feats0).to(device)
#   prot=torch.from_numpy(protos0).to(device)

#   print(ft.size(),prot.size())


#   #dist= torch.sum(torch.pow(ft.unsqueeze(1)-prot.unsqueeze(0),2.0),dim=2) #N,P,D
#   dist= -2*torch.mm(ft,prot.t())
#   dist+= torch.sum(torch.pow(ft,2.0),dim=1).unsqueeze(1) + torch.sum(torch.pow(prot,2.0),dim=1).unsqueeze(0)

#   return dist.cpu().numpy()


def run():

    ########
    ##
    # if you have less than 8 gbyte, then reduce from 250k
    ##
    ###############
    # 5000 instead of 250k for forloopdists
    feats = np.random.normal(size=(250000, 300))
    protos = np.random.normal(size=(500, 300))

    '''
	since = time.time()
	dists0=forloopdists(feats,protos)
	time_elapsed=float(time.time()) - float(since)
	print('Comp complete in {:.3f}s'.format( time_elapsed ))
	'''
    device = torch.device('cpu')

    since = time.time()
    dists1 = pytorchdists(feats, protos, device)
    time_elapsed = float(time.time()) - float(since)

    print('Comp complete in {:.3f}s'.format(time_elapsed))
    print(dists1.shape)

    since = time.time()
    dists2 = numpydists(feats, protos)
    time_elapsed = float(time.time()) - float(since)

    print('Comp complete in {:.3f}s'.format(time_elapsed))
    print(dists2.shape)
    print('df', np.max(np.abs(dists1-dists2)))


if __name__ == '__main__':
    run()

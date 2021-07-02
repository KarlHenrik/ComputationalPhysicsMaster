import os,sys,numpy as np

import torch

import time


class kmcomp():
  def __init__(self):
    pass

  def getinit(self,features,numprotos, init):
    if init == 'randselect':

      randinds = torch.randperm(self.n)[:numprotos]
      protos=torch.index_select(features, 0, randinds.to(features.device))
      return protos

    #elif: implement other ways to select here
    else:
      print( 'init is not randselect')
      raise

  def distance(self,features,centers):
    dist= -2*torch.mm(features,centers.t()) 
    dist+= torch.sum(torch.pow(features,2.0),dim=1).unsqueeze(1) + torch.sum(torch.pow(centers,2.0),dim=1).unsqueeze(0)
    return dist

  def cluster(self,features,numprotos,init='randselect',maxiter=100,convthresh=1e-3):

    self.n= features.shape[0]
    self.dim = features.shape[1]

    #check if all is ok
    if numprotos >= self.n:
      print('numprotos >= self.n')
      raise

    #initialize prototypes
    centers= self.getinit(features,numprotos, init)
  
    print(features.shape,centers.shape)

    for curiter in range(maxiter):

      #compute distances
      dists = self.distance(features,centers)
      #print('dists.shape', dists.shape)

      # for every sample find the nearest center
      inds = torch.argmin(dists,dim=1) 

      #update centers
      meanchange=0
      for c in range(numprotos):

        # find center samples such that c is nearest center:
        samplesbelongto_c = torch.nonzero(inds == c)
        #print(c,samplesbelongto_c.shape)
        newcenter = torch.mean(features[samplesbelongto_c,:])
        meanchange += torch.norm(centers[c,:] - newcenter) / float(numprotos)
        centers[c,:] = newcenter
      #print('centers.shape',centers.shape)
      print('iter {:d} meanchange {:.5f}'.format(curiter, meanchange.item()))

      if meanchange <= convthresh:
        #converged
        print('converged at iter:', curiter)
        break


def run():
  features=torch.tensor(np.random.normal(size=(100,2))) #.to('cuda:0')
  km = kmcomp()

  numprotos  = 10
  res = km.cluster(features,numprotos,init='randselect',maxiter=100,convthresh=1e-3)


if __name__=='__main__':
  run()

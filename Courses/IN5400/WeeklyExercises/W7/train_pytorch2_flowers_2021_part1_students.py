from __future__ import print_function, division

import torch
import torch.nn as nn
import torch.optim as optim
from torch.optim import lr_scheduler
from torch.autograd import Variable
import numpy as np
import torchvision
from torchvision import datasets, models, transforms
#import matplotlib.pyplot as plt
import time
import os

#import skimage.io
import PIL.Image

from torch.utils.data import Dataset, DataLoader
from torchvision import transforms, utils



class dataset_flowers(Dataset):
  def __init__(self, root_dir, trvaltest, transform=None):


    self.root_dir = root_dir

    self.transform = transform
    self.imgfilenames=[]
    self.labels=[]


    if trvaltest==0:
      #load training data
    elif trvaltest==1:
      #load validation data
    elif trvaltest==2:
      #load test data
      #TODO
    else:
      #TODO: print some error + exit() or an exception

  def __len__(self):
      #TODO
      return None

  def __getitem__(self, idx):

    #TODO
    sample = {'image': image, 'label': label, 'filename': filename}

    return sample


def train_epoch(model,  trainloader,  losscriterion, device, optimizer ):

    model.train() # IMPORTANT
 
    losses = list()
    for batch_idx, data in enumerate(trainloader):
      #TODO trains the model

    return np.mean(losses)


def evaluate_acc(model, dataloader, losscriterion, device):

    model.eval() # IMPORTANT

    losses = []
    curcount = 0
    accuracy = 0
    
    with torch.no_grad():
      for ctr, data in enumerate(dataloader):
        #TODO
        #computes predictions on samples from the dataloader
        # computes accuracy, to count how many samples, you can just sum up labels.shape[0]


    return accuracy.item() , np.mean(losses)


def train_model_nocv_sizes(dataloader_train, dataloader_test ,  model ,  losscriterion, optimizer, scheduler, num_epochs, device):

  best_measure = 0
  best_epoch =-1

  for epoch in range(num_epochs):
    print('Epoch {}/{}'.format(epoch, num_epochs - 1))
    print('-' * 10)

    model.train(True)
    losses=train_epoch(model,  dataloader_train,  losscriterion,  device , optimizer )

    if scheduler is not None:
      scheduler.step()

    model.train(False)
    measure, meanlosses = evaluate_acc(model, dataloader_test, losscriterion, device)
    print(' perfmeasure', measure)

    if measure > best_measure: #higher is better or lower is better?
      #TODO
      #bestweights = model.state_dict()
      #update   best_measure, best_epoch

  return best_epoch, best_measure, bestweights



def runstuff_finetunealllayers():

  #someparameters
  batchsize_tr=16
  batchsize_test=16
  maxnumepochs=2 #TODO can change when code runs here

  #TODO depends on what you can use
  device= torch.device('cpu') #torch.device('cuda')

  numcl=102
  #transforms
  data_transforms = {}
  data_transforms['train']=transforms.Compose([
          transforms.Resize(224),
          transforms.RandomCrop(224),
          transforms.RandomHorizontalFlip(),
          transforms.ToTensor(),
          transforms.Normalize([0.485, 0.456, 0.406], [0.229, 0.224, 0.225])
      ])
  data_transforms['val']=transforms.Compose([
          transforms.Resize(224),
          transforms.CenterCrop(224),
          transforms.RandomHorizontalFlip(),
          transforms.ToTensor(),
          transforms.Normalize([0.485, 0.456, 0.406], [0.229, 0.224, 0.225])
      ])
  
  #TODO
  datasets={}
  datasets['train']= None
  datasets['val']= None
  datasets['test']= None

  #TODO
  dataloaders={}
  dataloaders['train']= None
  dataloaders['val']= None
  dataloaders['test']=None
  

  #TODO
  #model
  model = None

  model.to(device)

  #TODO
  criterion = None

  lrates=[0.01, 0.001]

  best_hyperparameter= None
  weights_chosen = None
  bestmeasure = None
  for lr in lrates:

    #TODO
    optimizer = None
    
    best_epoch, best_perfmeasure, bestweights = train_model_nocv_sizes(dataloader_train = dataloaders['train'], dataloader_test = dataloaders['val'] ,  model = model ,  losscriterion = criterion , optimizer = optimizer, scheduler = None, num_epochs = maxnumepochs , device = device)

    if best_hyperparameter is None:
      best_hyperparameter = lr
      weights_chosen = bestweights
      bestmeasure = best_perfmeasure
     
    elif True:  #TODO what criterion here?
      #TODO
      pass

  model.load_state_dict(weights_chosen)

  accuracy, testloss = evaluate_acc(model = model , dataloader  = dataloaders['test'], losscriterion = criterion, device = device)

  print('accuracy val',bestmeasure , 'accuracy test',accuracy  )



def runstuff_finetunelastlayer():

  pass
  #TODO 


def runstuff_fromscratch():

  pass
  #TODO 




if __name__=='__main__':

  #runstuff_fromscratch()
  runstuff_finetunealllayers()
  #runstuff_finetunelastlayer()




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

    """
    Args:

        root_dir (string): Directory with all the images.
        transform (callable, optional): Optional transform to be applied
        on a sample.
    """

    self.root_dir = root_dir

    self.transform = transform
    self.imgfilenames=[]
    self.labels=[]

    if trvaltest==0:
      fn=os.path.join(self.root_dir,'trainfile.txt')
      with open(fn,'r') as f:
        for line in f: 
          v=line.rstrip().split()
          nm=os.path.join(self.root_dir, os.path.join('jpg',v[0]))
          #print(nm)
          self.imgfilenames.append(nm)
          self.labels.append(int(v[1]))
    if trvaltest==1:
      fn=os.path.join(self.root_dir,'valfile.txt')
      with open(fn,'r') as f:
        for line in f: 
          v=line.rstrip().split()
          nm=os.path.join(self.root_dir, os.path.join('jpg',v[0]))
          #print(nm)
          self.imgfilenames.append(nm)
          self.labels.append(int(v[1]))
    if trvaltest==2:
      fn=os.path.join(self.root_dir,'testfile.txt')
      with open(fn,'r') as f:
        for line in f: 
          v=line.rstrip().split()
          nm=os.path.join(self.root_dir, os.path.join('jpg',v[0]))
          #print(nm)
          self.imgfilenames.append(nm)
          self.labels.append(int(v[1]))


  def __len__(self):
      return len(self.imgfilenames)

  def __getitem__(self, idx):

    image = PIL.Image.open(self.imgfilenames[idx]).convert('RGB')

    if self.transform:
      image = self.transform(image)

    label=self.labels[idx]

    sample = {'image': image, 'label': label, 'filename': self.imgfilenames[idx]}

    return sample


def train_epoch(model,  trainloader,  losscriterion, device, optimizer ):

    model.train()
 
    losses = list()
    for batch_idx, data in enumerate(trainloader):

        inputs=data['image'].to(device)
        labels=data['label'].to(device)
        

        optimizer.zero_grad()

        output = model(inputs)
        loss = losscriterion(output, labels)
  
        loss.backward()
        optimizer.step()

        losses.append(loss.item())
        if batch_idx%100==0:
          print('current mean of losses ',np.mean(losses))

    return np.mean(losses)


def evaluate_acc(model, dataloader, losscriterion, device):

    model.eval()

    losses = []
    curcount = 0
    accuracy = 0
    
    with torch.no_grad():
      for ctr, data in enumerate(dataloader):
      
          inputs = data['image'].to(device)        
          outputs = model(inputs)

          labels = data['label']

          cpuout= outputs.to('cpu')
          
          loss = losscriterion(cpuout, labels)
          losses.append(loss.item())
          
          _, preds = torch.max(cpuout, 1)
          labels = labels.float()
          
          corrects = torch.sum(preds == labels.data) / float(labels.shape[0] )
          accuracy = accuracy*( curcount/ float(curcount+ labels.shape[0]) ) + corrects.float()* ( labels.shape[0]/ float(curcount+ labels.shape[0]) )
          curcount+= labels.shape[0]



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
      bestweights= model.state_dict()
      best_measure = measure
      best_epoch = epoch
      print('current best', best_measure, ' at epoch ', best_epoch)

  return best_epoch, best_measure, bestweights




def runstuff_fromscratch():

  #someparameters
  batchsize_tr=32
  batchsize_test=16
  maxnumepochs=5 

  device=torch.device('cuda')

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
  

  datasets={}
  datasets['train']= dataset_flowers(root_dir='./flowers_data',trvaltest=0, transform=data_transforms['train'])  
  datasets['val']= dataset_flowers(root_dir='./flowers_data',trvaltest=1, transform=data_transforms['val'])
  datasets['test']= dataset_flowers(root_dir='./flowers_data',trvaltest=2, transform=data_transforms['val'])


  dataloaders={}
  dataloaders['train']= torch.utils.data.DataLoader(datasets['train'], batch_size=batchsize_tr, shuffle=True) 
  dataloaders['val']= torch.utils.data.DataLoader(datasets['val'], batch_size=batchsize_test, shuffle=False )
  dataloaders['test']= torch.utils.data.DataLoader(datasets['test'], batch_size=batchsize_test, shuffle=False)  
  


  #model
  model = models.resnet18(pretrained=False)
  num_ftrs = model.fc.in_features
  model.fc = nn.Linear(num_ftrs, numcl)
  model.fc.reset_parameters()

  model.to(device)

  criterion = torch.nn.CrossEntropyLoss(weight=None, size_average=None, ignore_index=-100, reduce=None, reduction='mean')

  lrates=[0.01, 0.001]

  best_hyperparameter= None
  weights_chosen = None
  bestmeasure = None
  for lr in lrates:
    print("\nTrying learning rate = {}".format(lr))
    optimizer = optim.SGD(model.parameters(), lr=lr, momentum=0.9)
    best_epoch, best_perfmeasure, bestweights = train_model_nocv_sizes(dataloader_train = dataloaders['train'], dataloader_test = dataloaders['val'] ,  model = model ,  losscriterion = criterion , optimizer = optimizer, scheduler = None, num_epochs = maxnumepochs , device = device)

    if best_hyperparameter is None:
      best_hyperparameter = lr
      weights_chosen = bestweights
      bestmeasure = best_perfmeasure
    elif best_perfmeasure > bestmeasure:
      best_hyperparameter = lr
      weights_chosen = bestweights
      bestmeasure = best_perfmeasure

  model.load_state_dict(weights_chosen)

  accuracy, testloss  = evaluate_acc(model = model , dataloader  = dataloaders['test'], losscriterion = criterion, device = device)

  print('accuracy val',bestmeasure , 'accuracy test',accuracy  )




def runstuff_finetunealllayers():

  #someparameters
  batchsize_tr=24
  batchsize_test=16
  maxnumepochs=5 

  device=torch.device('cpu')

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
  

  datasets={}
  datasets['train']= dataset_flowers(root_dir='./flowers_data',trvaltest=0, transform=data_transforms['train'])  
  datasets['val']= dataset_flowers(root_dir='./flowers_data',trvaltest=1, transform=data_transforms['val'])
  datasets['test']= dataset_flowers(root_dir='./flowers_data',trvaltest=2, transform=data_transforms['val'])


  dataloaders={}
  dataloaders['train']= torch.utils.data.DataLoader(datasets['train'], batch_size=batchsize_tr, shuffle=True) 
  dataloaders['val']= torch.utils.data.DataLoader(datasets['val'], batch_size=batchsize_test, shuffle=False )
  dataloaders['test']= torch.utils.data.DataLoader(datasets['test'], batch_size=batchsize_test, shuffle=False)  
  


  #model
  model = models.resnet18(pretrained=True)
  num_ftrs = model.fc.in_features
  model.fc = nn.Linear(num_ftrs, numcl)
  model.fc.reset_parameters()

  model.to(device)

  criterion = torch.nn.CrossEntropyLoss(weight=None, size_average=None, ignore_index=-100, reduce=None, reduction='mean')

  lrates=[0.01, 0.001]

  best_hyperparameter= None
  weights_chosen = None
  bestmeasure = None
  for lr in lrates:
    print("\nTrying learning rate = {}".format(lr))
    optimizer = optim.SGD(model.parameters(), lr=lr, momentum=0.9)
    best_epoch, best_perfmeasure, bestweights = train_model_nocv_sizes(dataloader_train = dataloaders['train'], dataloader_test = dataloaders['val'] ,  model = model ,  losscriterion = criterion , optimizer = optimizer, scheduler = None, num_epochs = maxnumepochs , device = device)

    if best_hyperparameter is None:
      best_hyperparameter = lr
      weights_chosen = bestweights
      bestmeasure = best_perfmeasure
    elif best_perfmeasure > bestmeasure:
      best_hyperparameter = lr
      weights_chosen = bestweights
      bestmeasure = best_perfmeasure

  model.load_state_dict(weights_chosen)

  accuracy, testloss = evaluate_acc(model = model , dataloader  = dataloaders['test'], losscriterion = criterion, device = device)

  print('accuracy val',bestmeasure , 'accuracy test',accuracy  )

def runstuff_finetunelastlayer():

  #someparameters
  batchsize_tr=32
  batchsize_test=16
  maxnumepochs=5 

  device=torch.device('cuda')

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
  

  datasets={}
  datasets['train']= dataset_flowers(root_dir='./flowers_data',trvaltest=0, transform=data_transforms['train'])  
  datasets['val']= dataset_flowers(root_dir='./flowers_data',trvaltest=1, transform=data_transforms['val'])
  datasets['test']= dataset_flowers(root_dir='./flowers_data',trvaltest=2, transform=data_transforms['val'])


  dataloaders={}
  dataloaders['train']= torch.utils.data.DataLoader(datasets['train'], batch_size=batchsize_tr, shuffle=True) 
  dataloaders['val']= torch.utils.data.DataLoader(datasets['val'], batch_size=batchsize_test, shuffle=False )
  dataloaders['test']= torch.utils.data.DataLoader(datasets['test'], batch_size=batchsize_test, shuffle=False)  
  


  #model
  model = models.resnet18(pretrained=True)
  num_ftrs = model.fc.in_features
  model.fc = nn.Linear(num_ftrs, numcl)
  model.fc.reset_parameters()

  model.to(device)

  criterion = torch.nn.CrossEntropyLoss(weight=None, size_average=None, ignore_index=-100, reduce=None, reduction='mean')

  lrates=[0.01, 0.001]

  best_hyperparameter= None
  weights_chosen = None
  bestmeasure = None
  for lr in lrates:
    print("\nTrying learning rate = {}".format(lr))
    optimizer = optim.SGD(model.fc.parameters(), lr=lr, momentum=0.9)
    best_epoch, best_perfmeasure, bestweights = train_model_nocv_sizes(dataloader_train = dataloaders['train'], dataloader_test = dataloaders['val'] ,  model = model ,  losscriterion = criterion , optimizer = optimizer, scheduler = None, num_epochs = maxnumepochs , device = device)

    if best_hyperparameter is None:
      best_hyperparameter = lr
      weights_chosen = bestweights
      bestmeasure = best_perfmeasure
    elif best_perfmeasure > bestmeasure:
      best_hyperparameter = lr
      weights_chosen = bestweights
      bestmeasure = best_perfmeasure

  model.load_state_dict(weights_chosen)

  accuracy, testloss  = evaluate_acc(model = model , dataloader  = dataloaders['test'], losscriterion = criterion, device = device)

  print('accuracy val',bestmeasure , 'accuracy test',accuracy)

if __name__=='__main__':
  #runstuff_fromscratch()
  #runstuff_finetunealllayers()
  runstuff_finetunelastlayer()




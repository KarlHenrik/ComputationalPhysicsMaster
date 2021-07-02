import torch
import torch.nn as nn
import torch.optim as optim
from torch import Tensor
from torch.optim import lr_scheduler
from torch.autograd import Variable

import torchvision
from torchvision import datasets, models, transforms, utils
from torch.utils.data import Dataset, DataLoader
from torchvision.transforms import ToTensor

import time
import os
import numpy as np

import sklearn.metrics

from vocparseclslabels import PascalVOC

from typing import Callable, Optional

from PIL import Image

import matplotlib.pyplot as plt
import matplotlib.image as mpimg

# Dataset class for getting image tensors and labels of training and validation data
class dataset_voc(Dataset):
    def __init__(self, root_dir, trvaltest, transform=None):
        #self._X = []
        self._Y = []
        self._imgfilenames = []
        self._root_dir = root_dir
        self._transform = transform
        
        trvaltest = ["train", "val", "test"][trvaltest] # going from integer to string identifier
        
        pVOC = PascalVOC(root_dir) # object for getting the correct file names
        im_dict = {}
        for i, im_class in enumerate(pVOC.list_image_sets()):
            files = pVOC.imgs_from_category_as_list(im_class, trvaltest) # all files with class im_class in tr/val/test data
            for file in files:
                if file in im_dict:
                    im_dict[file][i] = 1
                else:
                    im_dict[file] = np.zeros(20)
                    im_dict[file][i] = 1

        for key in im_dict:
            #self._X.append(self._im2tensor(key))
            self._Y.append(im_dict[key])
            self._imgfilenames.append(key)

    def __len__(self):
        return len(self._imgfilenames)
    
    def getNumItems(self):
        return sum(len(item) for item in self._imgfilenames)

    def __getitem__(self, idx):
        image = self._im2tensor(self._imgfilenames[idx])
        #image = self._X[idx]
        label = torch.tensor(self._Y[idx]).float()
        sample = {"image": image, "label": label, "filename": self._imgfilenames[idx]}

        return sample
    
    def _im2tensor(self, im_name):
        img_path = self._root_dir + "JPEGImages/" + im_name +".jpg"
        image = Image.open(img_path)
        if self._transform:
            image = self._transform(image)
        else:
            transforms.ToTensor()(image)
        image = Variable(image)
        return image
    
def evaluate_meanavgprecision(model, dataloader, loss_func, device, numcl):
    model.eval() # Makes the forward pass more efficient for evaluation

    curcount = 0
    accuracy = 0
    idx_start = 0
    
    fnames = []  # filenames as they come out of the dataloader
    num_items = dataloader.dataset.getNumItems()
    pred_arr = torch.from_numpy(np.zeros(shape=(num_items, numcl))) # prediction scores for each class. each row is a list of scores. one score per image
    label_arr = torch.from_numpy(np.zeros(shape=(num_items, numcl))) # labels scores for each class. each row is a list of labels. one label per image

    losses = []
    avgprecs = np.zeros(numcl)  # average precision for each class

    with torch.no_grad(): # We need no gradients for evaluation, so this is faster
        for batch_idx, data in enumerate(dataloader):
            if (batch_idx % 100 == 0) and (batch_idx >= 100):
                print(f"at val batchindex: {batch_idx}")

            # Extracting data
            fname = data["filename"]
            inputs = data["image"].to(device)
            labels = data["label"].to(device)
            # Calculating what we need
            pred = model(inputs.to(device))
            loss = loss_func(pred, labels)
            # Storing values
            batch_size = len(fname)
            idx_end = idx_start + batch_size
            label_arr[idx_start:idx_end] = labels
            pred_arr[idx_start:idx_end] = pred
            idx_start += batch_size

            fnames += fname # concatenating lists
            losses.append(loss.item())

    for c in range(numcl):
        avgprecs[c] = sklearn.metrics.average_precision_score(y_true = label_arr[:, c], y_score = pred_arr[:, c])
        
    return avgprecs, np.mean(losses), label_arr, pred_arr, fnames

def traineval2_model_nocv(dataloader_train, dataloader_test,  model, loss_func, optimizer, scheduler, num_epochs, device, numcl):

    best_measure = 0
    best_epoch = -1

    trainlosses = []
    testlosses = []
    testperfs = []

    for epoch in range(num_epochs):
        print(f"Epoch {epoch}/{num_epochs - 1}")
        print("----------")

        avgloss = train_epoch(model, dataloader_train, loss_func, device, optimizer)
        trainlosses.append(avgloss)

        if scheduler is not None:
            scheduler.step()

        perfmeasure, testloss, concat_labels, concat_pred, fnames = evaluate_meanavgprecision(
            model, dataloader_test, loss_func, device, numcl)
        testlosses.append(testloss)
        testperfs.append(perfmeasure)

        print(f"at epoch: {epoch} classwise perfmeasure {perfmeasure}")

        avgperfmeasure = np.mean(perfmeasure)
        print(f"at epoch: {epoch} avgperfmeasure {avgperfmeasure}")

        if avgperfmeasure > best_measure:
            bestweights = model.state_dict()
            best_measure = avgperfmeasure
            best_epoch = epoch

    return best_epoch, best_measure, bestweights, trainlosses, testlosses, testperfs

def train_epoch(model, trainloader, loss_func, device, optimizer):
    model.train()

    losses = []
    for batch_idx, data in enumerate(trainloader):
        inputs = data["image"].to(device)
        labels = data["label"].to(device)

        # zero the parameter gradients
        optimizer.zero_grad()

        # forward + backward + optimize
        outputs = model(inputs)
        loss = loss_func(outputs, labels)
        loss.backward()
        optimizer.step()
        
        losses.append(loss.item())

    return np.mean(losses)

def runstuff():
    # --------------------- SETUP BLOCK --------------------------

    config = dict()

    # True #TODO change this to True for training on the cluster, eh
    config["use_gpu"] = True
    config["lr"] = 0.005
    config["batchsize_train"] = 16
    config["batchsize_val"] = 64
    config["maxnumepochs"] = 35

    config["scheduler_stepsize"] = 10
    config["scheduler_factor"] = 0.3

    # kind of a dataset property
    config["numcl"] = 20

    # data augmentations
    data_transforms = {
        "train": transforms.Compose([
            transforms.Resize(256),
            transforms.RandomCrop(224),
            transforms.RandomHorizontalFlip(),
            transforms.ToTensor(),
            transforms.Normalize([0.485, 0.456, 0.406], [0.229, 0.224, 0.225])
        ]),
        "val": transforms.Compose([
            transforms.Resize(224),
            transforms.CenterCrop(224),
            transforms.RandomHorizontalFlip(),
            transforms.ToTensor(),
            transforms.Normalize([0.485, 0.456, 0.406], [0.229, 0.224, 0.225])
        ]),
    }

    # datasets
    root_dir = "data/VOCdevkit/VOC2012/"
    image_datasets = {}
    image_datasets["train"] = dataset_voc(root_dir = root_dir, trvaltest=0, transform=data_transforms["train"])
    image_datasets["val"] = dataset_voc(root_dir = root_dir, trvaltest=1, transform=data_transforms["val"])

    # dataloaders
    dataloaders = {}
    dataloaders["train"] = DataLoader(image_datasets["train"], num_workers=1, batch_size = config["batchsize_train"])
    dataloaders["val"] = DataLoader(image_datasets["val"], num_workers=1, batch_size = config["batchsize_val"])

    # device
    if True == config["use_gpu"]:
        device = torch.device("cuda:0")
    else:
        device = torch.device("cpu")

    # model
    model = models.resnet18(pretrained=True) # pretrained resnet18
    model.fc = nn.Sequential( # rewriting last layer to have 20 sigmoid outputs
        nn.Linear(512, 20),
        nn.Sigmoid()
    )
    for param in model.parameters():
        param.requires_grad = False
    for param in model.fc.parameters():
        param.requires_grad = True

    model = model.to(device)

    #loss_func = BCE_custom()
    loss_func = nn.BCELoss()

    optimizer = optim.Adam(params = model.fc.parameters(), lr = config["lr"])

    # Decay LR by a factor of 0.3 every X epochs
    lr_sc = lr_scheduler.StepLR(optimizer, step_size = config["scheduler_stepsize"], gamma = config["scheduler_factor"])

    best_epoch, best_measure, bestweights, trainlosses, testlosses, testperfs = traineval2_model_nocv(
            dataloaders["train"], dataloaders["val"], model, loss_func, optimizer, lr_sc, num_epochs=config["maxnumepochs"], device=device, numcl=config["numcl"])
    
    
if __name__=='__main__':
    runstuff()
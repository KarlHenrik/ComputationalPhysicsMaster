import sys
import torch
import torch.nn as nn
import torch.optim as optim
from torch.optim import lr_scheduler

import torchvision
from torchvision import datasets, models, transforms, utils
from torch.utils.data import Dataset, DataLoader
#import matplotlib.pyplot as plt

from torch import Tensor

import time
import os
import numpy as np
import random

import PIL.Image
import sklearn.metrics
import matplotlib.pyplot as plt

from vocparseclslabels import PascalVOC

from typing import Callable, Optional

seed = 42
torch.manual_seed(seed)
torch.cuda.manual_seed(seed)
np.random.seed(seed)
random.seed(seed)
torch.manual_seed(seed)
torch.backends.cudnn.benchmark = False
torch.backends.cudnn.deterministic = True

############# Stuff ## ###############################################
TRAINING = False # Set false if you only want to evaluate
ROOT = '/itf-fi-ml/shared/IN5400/dataforall/mandatory1/VOCdevkit/VOC2012/'
######################################################################

class dataset_voc(Dataset):
    def __init__(self, root_dir, trvaltest, transform=None):
        """
        Args:
            root_dir (string): path to main file
            trvaltest (int): train=0, val=1, test=2
            transform (torchvision.transforms): transforms to be applied to data
        """
        if trvaltest == 0:
            dataset = 'train'
        elif trvaltest == 1:
            dataset = 'val'
        elif trvaltest == 2:
            dataset = 'test'
        else:
            raise Exception("Strange...")

        self.root_dir = root_dir
        # Dict that maps file names to labels
        self.imgfilenames = {}
        self.transform = transform

        pv = PascalVOC(self.root_dir)
        categories = pv.list_image_sets()

        # Add all images to dict
        for cat in categories:
            for img in pv.imgs_from_category_as_list(cat, dataset):
                self.imgfilenames[img] = np.zeros(len(categories))

        # Label images
        for idx in range(len(categories)):
            imgs_from_cat = pv.imgs_from_category_as_list(
                categories[idx], dataset)
            for img in imgs_from_cat:
                # Update labels
                self.imgfilenames[img][idx] = 1

        # Cast to list, so it is possible to index it
        self.imgfilenames = list(self.imgfilenames.items())

    def __len__(self):
        return len(self.imgfilenames)

    def __getitem__(self, idx):
        image = PIL.Image.open(
            self.root_dir + "JPEGImages/"  +  self.imgfilenames[idx][0] + ".jpg")
        image = self.transform(image)
        label = self.imgfilenames[idx][1]

        sample = {'image': image, 'label': label,
                  'filename': self.imgfilenames[idx][0]}

        return sample


def train_epoch(model,  trainloader,  criterion, device, optimizer):

    model.train()
    losses = []
    running_loss = 0
    for batch_idx, data in enumerate(trainloader):
        optimizer.zero_grad()

        inputs = data['image'].to(device)
        labels = data['label'].to(device)

        outputs = model(inputs)
        loss = criterion(outputs, labels)
        if loss.isnan().any():
            print("NaN found")
            sys.exit()
        loss.backward()
        optimizer.step()

        losses.append(loss.item())

        running_loss += loss.item()

    return np.mean(losses)


def evaluate_meanavgprecision(model, dataloader, criterion, device, numcl):
    print("Evaluating...")
    model.eval()

    prediction_labels = [[] for _ in range(numcl)]
    prediction_scores = [[] for _ in range(numcl)]

    # prediction scores for each class. each numpy array is a list of scores. one score per image
    # Appending to python list should be faster no?
    prediction_scores = [[] for _ in range(numcl)]
    # concat_pred = [np.empty(shape=(0)) for _ in range(numcl)]
    # labels scores for each class. each numpy array is a list of labels. one label per image
    prediction_labels = [[] for _ in range(numcl)]
    # concat_labels = [np.empty(shape=(0)) for _ in range(numcl)]
    avgprecs = np.zeros(numcl)  # average precision for each class
    fnames = []  # filenames as they come out of the dataloader

    with torch.no_grad():
        losses = []
        for batch_idx, data in enumerate(dataloader):

            if (batch_idx % 100 == 0) and (batch_idx >= 100):
                print('at val batchindex: ', batch_idx)

            inputs = data['image'].to(device)
            outputs = model(inputs)

            labels = data['label'].to(device)

            loss = criterion(outputs, labels.to(device))
            losses.append(loss.item())

            cpuout = outputs.to('cpu')
            # Collect scores and label to later calulate average precision
            for i in range(labels.shape[0]):
                # For every output in batch
                single_output = cpuout[i]
                label = labels[i]
                fnames.append(data['filename'][i])
                for j in range(numcl):
                    # Add prediction score for every class
                    prediction_scores[j].append(single_output[j].cpu().numpy())
                    prediction_labels[j].append(label[j].cpu().numpy())

            """this was an accuracy computation
            cpuout= outputs.to('cpu')
            _, preds = torch.max(cpuout, 1)
            labels = labels.float()
            corrects = torch.sum(preds == labels.data)
            accuracy = accuracy*( curcount/ float(curcount+labels.shape[0]) ) + corrects.float()* ( curcount/ float(curcount+labels.shape[0]) )
            curcount += labels.shape[0]"""
    for c in range(numcl):
        avgprecs[c] = sklearn.metrics.average_precision_score(y_true=np.array(
            prediction_labels[c]), y_score=np.array(prediction_scores[c]))

    write_to_file(prediction_scores, fnames)

    return avgprecs, np.mean(losses), prediction_labels, prediction_scores, fnames


def write_to_file(prediction_scores, fnames):
    f = open("results.txt", "w")
    for i in range(20):
        for j in range(len(fnames)):
            if prediction_scores[i][j] > 0:
                f.write(f"{fnames[j]} {prediction_scores[i][j]}\n")
        f.write(",\n")
    f.close()


def traineval2_model_nocv(dataloader_train, dataloader_test,  model,  criterion, optimizer, scheduler, num_epochs, device, numcl):
    start_time = int(time.time())
    best_measure = 0
    best_epoch = -1

    trainlosses = []
    testlosses = []
    testperfs = []

    for epoch in range(num_epochs):
        print('Epoch {}/{}'.format(epoch, num_epochs - 1))
        print('-' * 10)
        avgloss = train_epoch(model,  dataloader_train,
                              criterion,  device, optimizer)
        trainlosses.append(avgloss)

        if scheduler is not None:
            scheduler.step()

        perfmeasure, testloss, concat_labels, concat_pred, fnames = evaluate_meanavgprecision(
            model, dataloader_test, criterion, device, numcl)
        testlosses.append(testloss)
        testperfs.append(perfmeasure)


        print('at epoch: ', epoch, ' classwise perfmeasure ', perfmeasure)

        avgperfmeasure = np.mean(perfmeasure)
        print('at epoch: ', epoch, ' avgperfmeasure ', avgperfmeasure)


        if avgperfmeasure > best_measure:  # higher is better or lower is better?
            bestweights = model.state_dict()
            # Tracking current best performance measure and epoch
            best_measure = avgperfmeasure
            best_epoch = epoch
            # Saving weights and scores
            torch.save({
                'weights': bestweights,
                'AP': avgperfmeasure,
                'epoch': epoch}, f"./models/best{start_time}.pt")

    return best_epoch, best_measure, bestweights, trainlosses, testlosses, testperfs


class yourloss(nn.modules.loss._Loss):

    def __init__(self, reduction: str = 'mean') -> None:
        """ 
        Take a binary cross entropy loss independently 
        for 20 classes and sum/average it, depending on reduction
        """
        super(yourloss, self).__init__()

        def sigmoid(logits):
            return torch.where(logits >= 0, 1/(1+torch.exp(-logits)),
                               torch.exp(logits)/(1+torch.exp(logits)))

        def BCELoss(logits, labels):
            probabilities = sigmoid(logits)
            return -(labels * torch.log(probabilities) + (1-labels) * torch.log(1-probabilities))

        self.loss = BCELoss

    def forward(self, input_: Tensor, target: Tensor) -> Tensor:
        loss = self.loss(input_, target)
        return loss.mean()


def runstuff():
    config = dict()

    # True #TODO change this to True for training on the cluster, eh
    config['use_gpu'] = True
    config['lr'] = 0.0005
    config['batchsize_train'] = 128
    config['batchsize_val'] = 64
    config['maxnumepochs'] = 12

    config['scheduler_stepsize'] = 2
    config['scheduler_factor'] = 0.3

    # kind of a dataset property
    config['numcl'] = 20

    # data augmentations
    data_transforms = {
        'train': transforms.Compose([
            transforms.Resize(256),
            transforms.RandomCrop(224),
            transforms.RandomHorizontalFlip(),
            transforms.ToTensor(),
            transforms.Normalize([0.485, 0.456, 0.406], [0.229, 0.224, 0.225])
        ]),
        'val': transforms.Compose([
            transforms.Resize(224),
            transforms.CenterCrop(224),
            transforms.RandomHorizontalFlip(),
            transforms.ToTensor(),
            transforms.Normalize([0.485, 0.456, 0.406], [0.229, 0.224, 0.225])
        ]),
    }

    # datasets
    image_datasets = {}
    image_datasets['train'] = dataset_voc(
        root_dir=ROOT, trvaltest=0, transform=data_transforms['train'])
    image_datasets['val'] = dataset_voc(
        root_dir=ROOT, trvaltest=1, transform=data_transforms['val'])
    
    # dataloaders
    dataloaders = {}
    dataloaders['train'] = DataLoader(
        image_datasets['train'], config['batchsize_train'], num_workers=1, shuffle=True)
    dataloaders['val'] = DataLoader(
        image_datasets['val'], config['batchsize_val'], num_workers=1)

    # device
    if True == config['use_gpu']:
        device = torch.device('cuda:0')

    else:
        device = torch.device('cpu')

    # model
    model = models.resnet18(pretrained=True)  # pretrained resnet18
    # overwrite last linear layer
    model.fc = nn.Linear(model.fc.in_features, 20)
    model = model.to(device)

    lossfct = yourloss()

    # Observe that all parameters are being optimized
    someoptimizer = optim.Adam(model.parameters(), lr=config['lr'])

    # Decay LR by a factor of 0.3 every X epochs
    somelr_scheduler = optim.lr_scheduler.StepLR(
        someoptimizer, config['scheduler_stepsize'], config['scheduler_factor'])

    training = True
    if TRAINING == False:  # Produce stuff that is needed for report
        # This loads best model, evaluates, and writes results to file to be used for GUI

        model.load_state_dict(torch.load(
            "models/best1614952859.pt")['weights'])  # Load best weights
        perfmeasure, testloss, concat_labels, concat_pred, fnames = evaluate_meanavgprecision(
            model, dataloaders['val'], lossfct, device, config['numcl'])
        write_to_file(concat_pred, fnames)

        max_t = np.max(np.array(concat_pred), axis=1).min()
        t_dist = np.arange(20) * max_t/20
        tailacc = np.zeros(20)
        for c in range(20):
            scores = np.array(concat_pred[c])
            labels = np.array(concat_labels[c])
            for t_idx in range(len(t_dist)):
                tailacc[t_idx] += 1/(scores[scores > t_dist[t_idx]]
                                     ).shape[0] * labels[scores > t_dist[t_idx]].sum()

        print(f"AP Scores: {perfmeasure}, mAP: {perfmeasure.mean()}")
        plt.plot(t_dist, tailacc/20)
        plt.xlabel("t")
        plt.ylabel("Tailacc(t)")
        plt.show()

    elif TRAINING == True:
        best_epoch, best_measure, bestweights, trainlosses, testlosses, testperfs = traineval2_model_nocv(
            dataloaders['train'], dataloaders['val'],  model,  lossfct, someoptimizer, somelr_scheduler, num_epochs=config['maxnumepochs'], device=device, numcl=config['numcl'])


if __name__ == '__main__':
    runstuff()

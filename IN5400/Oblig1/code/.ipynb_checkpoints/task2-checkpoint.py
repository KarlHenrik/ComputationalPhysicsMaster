# for training av validation
import torch
import torch.nn as nn
import torch.optim as optim

from torchvision import datasets, models, transforms
from torch.utils.data import Dataset, DataLoader
from torchvision.transforms import ToTensor

import copy
import time
import numpy as np

import sklearn.metrics

from vocparseclslabels import PascalVOC

from PIL import Image
# custom functions
from vocTrainEval import dataset_voc, evaluate_meanavgprecision, traineval2_model_nocv, train_epoch

# --------------- Mostly unchanged given functions --------------------
# change targetmodule to value
def setbyname2(targetmodel, name, value):
    def iteratset(obj, components, value, nametail=[]):

        if not hasattr(obj, components[0]):
            return False
        elif len(components) == 1:
            if not hasattr(obj, components[0]):
                print('object has not the component:', components[0])
                print('nametail:', nametail)
                exit()
            setattr(obj, components[0], value)
            #print('found!!', components[0])
            # exit()
            return True
        else:
            nextobj = getattr(obj, components[0])

            newtail = nametail
            newtail.append(components[0])
            #print('components ',components, nametail, newtail)
            # print(type(obj),type(nextobj))

            return iteratset(nextobj, components[1:], value, nametail=newtail)

    components = name.split('.')
    success = iteratset(targetmodel, components, value, nametail=[])
    return success

# compare outputs of two models
def comparetwomodeloutputs(model1, model2, dataloader, device):

    model1.eval()
    model2.eval()

    curcount = 0
    avgdiff = 0

    with torch.no_grad():
        for batch_idx, data in enumerate(dataloader):
            inputs = data['image'].to(device)
            outputs1 = model1(inputs)
            outputs2 = model2(inputs)

            diff = torch.mean(torch.abs((outputs1-outputs2).flatten()))

            labels = data['label']
            avgdiff = avgdiff * (curcount / float(curcount+labels.shape[0])) + diff.item() * (
                labels.shape[0] / float(curcount+labels.shape[0]))

            curcount += labels.shape[0]

    return avgdiff

# ----------------------- Network change ---------------
class wsconv2(nn.Conv2d):
    def __init__(self, in_channels, out_channels, kernel_size, stride,
                 padding, weight, dilation, groups, bias, eps=1e-12):
        super(wsconv2, self).__init__(in_channels, out_channels,
                                      kernel_size, stride, padding, dilation, groups, bias)

        self.weight = weight
        self.eps = eps

    def forward(self, x):
        # torch.nn.functional.conv2d documentation tells about weight shapes
        n_c = torch.sqrt(torch.var(self.weight, axis = (1, 2, 3)) + self.eps)
        scaled = torch.zeros_like(self.weight)
        for i in range(len(n_c)):
            scaled[i] = self.weight[i] / n_c[i]
            
        return torch.nn.functional.conv2d(x, weight = scaled, stride = self.stride, padding = self.padding, dilation = self.dilation, groups = self.groups)

def bntoWSconverter(model, case):
    lastwasconv2 = False
    for nm, module in model.named_modules():
        if isinstance(module, nn.Conv2d):
            lastwasconv2 = True

            usedeps = 1e-12  # use 1e-12 if you add it to a variance term, and 1e-6 if you add it to a standard deviation term
            
            new_module = wsconv2(module.in_channels, module.out_channels, module.kernel_size, module.stride,
                 module.padding, module.weight, module.dilation, module.groups, module.bias, usedeps)
            
            setbyname2(model, nm, new_module)
            n_c = torch.sqrt(torch.var(module.weight.detach(), axis = (1, 2, 3)) + usedeps) # for next batchnorm

        elif isinstance(module, nn.BatchNorm2d):

            if lastwasconv2 == False:
                print('got disconnected batchnorm??')
                exit()
                
            if case == "A":
                module.running_mean /= n_c
                module.running_var /= n_c**2
            elif case == "B":
                module.bias = torch.nn.Parameter(module.bias + (module.weight.detach() * n_c - module.weight.detach()) * module.running_mean / (module.running_var + 1e-12)**0.5)
                module.weight = torch.nn.Parameter(module.weight * n_c)

            lastwasconv2 = False

        else:
            lastwasconv2 = False
            
# ------------------ Setup ----------------
# routine to test that your copied model at evaluation time works as intended

config = dict()

config['use_gpu'] = True
config['lr'] = 0.0005
config['batchsize_train'] = 16
config['batchsize_val'] = 64
config["maxnumepochs"] = 2

config["scheduler_stepsize"] = 2
config["scheduler_factor"] = 0.3
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
root_dir = "VOCdevkit/VOC2012/"
image_datasets = {}
image_datasets["train"] = dataset_voc(root_dir = root_dir, trvaltest=0, transform=data_transforms["train"])
image_datasets["val"] = dataset_voc(root_dir = root_dir, trvaltest=1, transform=data_transforms["val"])

# dataloaders
dataloaders = {}
dataloaders["train"] = DataLoader(image_datasets["train"], num_workers=0, batch_size = config["batchsize_train"], shuffle=True)
dataloaders["val"] = DataLoader(image_datasets["val"], num_workers=0, batch_size = config["batchsize_val"])

# device
if True == config["use_gpu"]:
    device = torch.device("cuda:0")
else:
    device = torch.device("cpu")
    
# model
model = models.resnet18(pretrained=True)
model.fc = nn.Sequential( # rewriting last layer to have 20 sigmoid outputs
    nn.Linear(model.fc.in_features, 20),
    nn.Sigmoid()
)

# ------------------ Testing -------------------
for case in ["A", "B"]:
    model2 = copy.deepcopy(model.to(device))

    bntoWSconverter(model2, case = case) # changes model2 inplace
    model = model.to(device)
    model2 = model2.to(device)

    avgdiff = comparetwomodeloutputs(model, model2, dataloaders["val"], device)

    # order 1e-3 is okay, 1e-2 is still okay.
    print(f"model checking averaged difference {avgdiff} for case {case}")
    
# -------------------- Gave up on training ----------------
"""
Training is super slow and bad for both cases

loss_func = nn.BCELoss()

optimizer = optim.Adam(params = model.parameters(), lr = config["lr"])

lr_sc = optim.lr_scheduler.StepLR(optimizer, step_size = config["scheduler_stepsize"], gamma = config["scheduler_factor"]) # Decay LR by a factor of 0.3 every X epochs

best_all, trainlosses, testlosses, testperfs = traineval2_model_nocv(
        dataloaders["train"], dataloaders["val"], model2, loss_func, optimizer, lr_sc, num_epochs=config["maxnumepochs"], device=device, numcl=config["numcl"])
"""
# for training av validation
import torch
import torch.nn as nn
import torch.optim as optim

from torchvision import datasets, models, transforms
from torch.utils.data import Dataset, DataLoader
from torchvision.transforms import ToTensor

import time
import numpy as np

import sklearn.metrics

from vocparseclslabels import PascalVOC

from PIL import Image
# custom functions
from vocTrainEval import dataset_voc, evaluate_meanavgprecision, traineval2_model_nocv, train_epoch
# for plotting
import matplotlib.pyplot as plt

# --------------------- SETUP BLOCK --------------------------
config = dict()

config["use_gpu"] = True
config["lr"] = 0.0005
config["batchsize_train"] = 16
config["batchsize_val"] = 64
config["maxnumepochs"] = 10

config["scheduler_stepsize"] = 2
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
    
model = models.resnet18(pretrained=True) # pretrained resnet18
model.fc = nn.Sequential( # rewriting last layer to have 20 sigmoid outputs
    nn.Linear(model.fc.in_features, 20),
    nn.Sigmoid()
)
model = model.to(device)

loss_func = nn.BCELoss()

optimizer = optim.Adam(params = model.parameters(), lr = config["lr"])

lr_sc = optim.lr_scheduler.StepLR(optimizer, step_size = config["scheduler_stepsize"], gamma = config["scheduler_factor"]) # Decay LR by a factor of 0.3 every X epochs
# training
best_all, trainlosses, testlosses, testperfs = traineval2_model_nocv(
        dataloaders["train"], dataloaders["val"], model, loss_func, optimizer, lr_sc, num_epochs=config["maxnumepochs"], device=device, numcl=config["numcl"])
# saving best model and its results
modelfile = "models/bestMdl.pt"
torch.save(best_all, modelfile)
print(f"model saved to {modelfile}")

# Plotting results from training

# train and test losses
plt.plot(trainlosses, label="Train")
plt.plot(testlosses, label="Test")
plt.legend()
plt.title("Train and test losses")
plt.xlabel("Epoch")
plt.show()
# test mAP
plt.plot(np.mean(testperfs, axis=1))
plt.title("Test mAP")
plt.xlabel("Epoch")
plt.show()
# tail accuracy
def tailAcc(preds, labels, t):
    preds_true = np.greater(preds, t) # predicted to be true
    correct_pred = np.logical_and(preds_true, labels) # correct predictions!
    
    return np.sum(correct_pred) / np.sum(preds_true)

preds = best_all["preds"].cpu().numpy()
labels = best_all["labels"].cpu().numpy()

max_t = np.min(np.max(preds, axis = 0)) - 1e-6
ts = np.linspace(0, max_t, 20)
mean_tail_accs = np.zeros(20) # mean tail acc for each t

for i, t in enumerate(ts):
    tail_accs = np.zeros(config["numcl"]) # tail accs for each class for this t
    for cl in range(config["numcl"]):
        tail_accs[cl] = tailAcc(preds[:, cl], labels[:, cl], t)
    
    mean_tail_accs[i] = np.mean(tail_accs) # mean tail acc for this t

plt.plot(ts, mean_tail_accs)
plt.title("Tailacc averaged over all classes")
plt.xlabel("t")
plt.show()

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

# Dataset class for getting image tensors and labels of training and validation data
class dataset_voc(Dataset):
    def __init__(self, root_dir, trvaltest, transform=None):
        self._labels = []
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
            self._labels.append(im_dict[key])
            self._imgfilenames.append(key)

    def __len__(self):
        return len(self._imgfilenames)

    def __getitem__(self, idx):
        img_path = self._root_dir + "JPEGImages/" + self._imgfilenames[idx] +".jpg"
        image = Image.open(img_path)
        if self._transform:
            image = self._transform(image)
        else:
            transforms.ToTensor()(image)
        
        label = torch.from_numpy(self._labels[idx]).float()
        sample = {"image": image, "label": label, "filename": self._imgfilenames[idx]}

        return sample
    
def evaluate_meanavgprecision(model, dataloader, loss_func, device, numcl):
    model.eval() # Makes the forward pass more efficient for evaluation

    curcount = 0
    accuracy = 0
    idx_start = 0
    
    num_items = len(dataloader.dataset)
    pred_arr = torch.from_numpy(np.zeros(shape=(num_items, numcl))).to(device) # prediction scores for each class. each row is a list of scores. one score per image
    label_arr = torch.from_numpy(np.zeros(shape=(num_items, numcl))).to(device) # labels scores for each class. each row is a list of labels. one label per image
    fnames = []  # filenames as they come out of the dataloader
    losses = []
    
    with torch.no_grad(): # We need no gradients for evaluation, so this is faster
        for batch_idx, data in enumerate(dataloader):
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
            
    avgprecs = np.zeros(numcl)  # average precision for each class
    for c in range(numcl):
        avgprecs[c] = sklearn.metrics.average_precision_score(y_true = label_arr[:, c].to("cpu"), y_score = pred_arr[:, c].to("cpu"))
        
    return avgprecs, np.mean(losses), label_arr, pred_arr, fnames

def traineval2_model_nocv(dataloader_train, dataloader_test,  model, loss_func, optimizer, scheduler, num_epochs, device, numcl):
    best_measure = 0
    best_epoch = -1

    trainlosses = []
    testlosses = []
    testperfs = []

    for epoch in range(num_epochs):
        print(f"Epoch {epoch + 1}/{num_epochs}")
        avgloss = train_epoch(model, dataloader_train, loss_func, device, optimizer) # training
        if scheduler is not None:
            scheduler.step()
        perfmeasure, testloss, label_arr, pred_arr, fnames = evaluate_meanavgprecision(model, dataloader_test, loss_func, device, numcl) # evaluation
        
        trainlosses.append(avgloss)
        testlosses.append(testloss)
        testperfs.append(perfmeasure)
        
        avgperfmeasure = np.mean(perfmeasure)
        if avgperfmeasure > best_measure: # if this epoch is best, store all results from this epoch
            best_measure = avgperfmeasure
            best_all = {"preds": pred_arr,
                        "labels": label_arr,
                        "fnames": fnames,
                        "epoch": epoch,
                        "weights": model.state_dict(),
                        "cwAP": perfmeasure}
            
        print(f"Classwise perfmeasure: {[np.round(pm, 3) for pm in perfmeasure]}")
        print(f"Avgperfmeasure: {avgperfmeasure:.3f}, train {avgloss:.3f}, test {testloss:.3f}\n----------")

    return best_all, trainlosses, testlosses, testperfs

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
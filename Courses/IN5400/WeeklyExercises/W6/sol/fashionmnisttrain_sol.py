import torch
import torch.nn as nn
import torch.nn.functional as F
from utils.utility_functions import datasetFashionMNIST
import numpy as np
import matplotlib.pyplot as plt
from sklearn import metrics


import os
#%matplotlib inline


def datasetup():

  config = {
            'batch_size': 32,
            'use_cuda': True,       #True=use Nvidia GPU | False use CPU
            'log_interval': 100,    #How often to display (batch) loss during training
            'epochs': 10,          #Number of epochs
            'learningRate': 0.001
           }
           


  dataPath = './data/fashionmnist/' #'/projects/in5400/MNIST_fashion/'

  # Create dataset objects
  train_dataset = datasetFashionMNIST(dataPath=dataPath, train=True)
  val_dataset   = datasetFashionMNIST(dataPath=dataPath, train=False)

  train_loader = torch.utils.data.DataLoader(train_dataset, batch_size=config['batch_size'], shuffle=True, num_workers=1)
  val_loader   = torch.utils.data.DataLoader(val_dataset, batch_size=config['batch_size'], shuffle=False, num_workers=1)


  classes = ['T-shirt / top', 'Trouser', 'Pullover', 'Dress', 'Coat', 'Sandal', 'Shirt', 'Sneaker', 'Bag', 'Ankle boot']

  return config, train_dataset, val_dataset, train_loader, val_loader, classes


class Model(nn.Module):
    def __init__(self):
        super(Model, self).__init__()
        self.cnn1 = nn.Conv2d(in_channels=1,  out_channels=16, kernel_size=3, stride=1, padding=1)
        self.cnn2 = nn.Conv2d(in_channels=16, out_channels=16, kernel_size=3, stride=1, padding=1)
        self.cnn3 = nn.Conv2d(in_channels=16, out_channels=16, kernel_size=3, stride=1, padding=1)
        self.maxpool1 = nn.MaxPool2d(kernel_size=2)
        self.cnn4 = nn.Conv2d(in_channels=16, out_channels=32, kernel_size=3, stride=1, padding=1)
        self.cnn5 = nn.Conv2d(in_channels=32, out_channels=32, kernel_size=3, stride=1, padding=1)
        self.cnn6 = nn.Conv2d(in_channels=32, out_channels=32, kernel_size=3, stride=1, padding=1)
        self.maxpool2 = nn.MaxPool2d(kernel_size=2)

        self.act= nn.ReLU()#nn.LeakyReLU(negative_slope=0.1)#nn.ReLU() #nn.Tanh()#nn.ReLU()


        self.fc1 = nn.Linear(32*7*7, 10)

    def forward(self, x):
        out = x.reshape(x.size(0), 1, 28, 28)
        out = self.act(self.cnn1(out))
        out = self.act(self.cnn2(out))
        out = self.act(self.cnn3(out))
        out = self.maxpool1(out)
        out = self.act(self.cnn4(out))
        out = self.act(self.cnn5(out))
        out = self.act(self.cnn6(out))
        out = self.maxpool2(out)
        out = out.view(out.size(0), -1)
        out = self.fc1(out)
        return out
        
        

    
    
def run_epoch(model, epoch, data_loader, optimizer, loss_fct, is_training, config, attach_hook):
    """
    Args:
        model        (obj): The neural network model
        epoch        (int): The current epoch
        data_loader  (obj): A pytorch data loader "torch.utils.data.DataLoader"
        optimizer    (obj): A pytorch optimizer "torch.optim"
        is_training (bool): Whether to use train (update) the model/weights or not. 
        config      (dict): Configuration parameters

    Intermediate:
        totalLoss: (float): The accumulated loss from all batches. 
                            Hint: Should be a numpy scalar and not a pytorch scalar

    Returns:
        loss_avg         (float): The average loss of the dataset
        accuracy         (float): The average accuracy of the dataset
        confusion_matrix (float): A 10x10 matrix
    """

    if is_training==True: 
        model.train()
    else:
        model.eval()

    total_loss       = 0 
    correct          = 0 
    confusion_matrix = np.zeros(shape=(10,10))
    labels_list      = [0,1,2,3,4,5,6,7,8,9]

    for batch_idx, data_batch in enumerate(data_loader):
        if config['use_cuda'] == True:
            images = data_batch[0].to('cuda') # send data to GPU
            labels = data_batch[1].to('cuda') # send data to GPU
        else:
            images = data_batch[0]
            labels = data_batch[1]

        if not is_training:
            with torch.no_grad():
                prediction = model(images)
                loss        = loss_fct(prediction, labels)
                total_loss += loss.item()    
            
        elif is_training:
        
            if attach_hook and (batch_idx<50):# we dont need to attach gross amounts of gradients
                handles=[]
                outputpath= './tmp/'
                filestub= 'batchindex{:d}'.format(batch_idx)
                for ind,(name,module) in enumerate(model.named_modules()):
                    print('name: {}'.format(name) )
                    if isinstance(module,torch.nn.Conv2d):
                        h=module.register_backward_hook( hook_factory(outputpath = outputpath , filestub = filestub+'_'+name ))
                        handles.append(h)   
        
            prediction = model(images)
            loss        = loss_fct(prediction, labels)
            total_loss += loss.item()

            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
            
            
            if attach_hook:
                # remove handles, as they are tied to a specific filename
                for h in handles:
                    h.remove()
                handles=[]
            

        # Update the number of correct classifications and the confusion matrix
        predicted_label  = prediction.max(1, keepdim=True)[1][:,0]
        correct          += predicted_label.eq(labels).cpu().sum().numpy()
        confusion_matrix += metrics.confusion_matrix(labels.cpu().numpy(), predicted_label.cpu().numpy(), labels=labels_list)

        # Print statistics
        #batchSize = len(labels)
        if batch_idx % config['log_interval'] == 0:
            print(f'Epoch={epoch} | {(batch_idx+1)/len(data_loader)*100:.2f}% | loss = {loss:.5f}')

    loss_avg         = total_loss / len(data_loader)
    accuracy         = correct / len(data_loader.dataset)
    confusion_matrix = confusion_matrix / len(data_loader.dataset)

    return loss_avg, accuracy, confusion_matrix


def loss_fn(prediction, labels):
    """Returns softmax cross entropy loss."""
    loss = F.cross_entropy(input=prediction, target=labels)
    return loss

def run():

    config, train_dataset, val_dataset, train_loader, val_loader, classes = datasetup()

    model = Model()
    if config['use_cuda'] == True:
        model.to('cuda')
        
    optimizer = torch.optim.SGD(model.parameters(), lr = config['learningRate'])        


    train_loss = np.zeros(shape=config['epochs'])
    train_acc  = np.zeros(shape=config['epochs'])
    val_loss   = np.zeros(shape=config['epochs'])
    val_acc    = np.zeros(shape=config['epochs'])
    train_confusion_matrix = np.zeros(shape=(10,10,config['epochs']))
    val_confusion_matrix   = np.zeros(shape=(10,10,config['epochs']))

    for epoch in range(config['epochs']):
    
        if epoch+1 == config['epochs']:
            attach_hook = True
        else:
            attach_hook = False      
        
        train_loss[epoch], train_acc[epoch], train_confusion_matrix[:,:,epoch] = \
                                   run_epoch(model, epoch, train_loader, optimizer, loss_fct=loss_fn, is_training=True, config=config, attach_hook = attach_hook)

        val_loss[epoch], val_acc[epoch], val_confusion_matrix[:,:,epoch]     = \
                                   run_epoch(model, epoch, val_loader, optimizer, loss_fct=loss_fn, is_training=False, config=config, attach_hook = attach_hook)



    doplot=False
    
    if doplot:

      #plot losses
      plt.figure(figsize=(18, 16), dpi= 80, facecolor='w', edgecolor='k')
      ax = plt.subplot(2, 1, 1)
      # plt.subplots_adjust(hspace=2)
      ax.plot(train_loss, 'b', label='train loss')
      ax.plot(val_loss, 'r', label='validation loss')
      ax.grid()
      plt.ylabel('Loss', fontsize=18)
      plt.xlabel('Epochs', fontsize=18)
      ax.legend(loc='upper right', fontsize=16)

      ax = plt.subplot(2, 1, 2)
      plt.subplots_adjust(hspace=0.4)
      ax.plot(train_acc, 'b', label='train accuracy')
      ax.plot(val_acc, 'r', label='validation accuracy')
      ax.grid()
      plt.ylabel('Accuracy', fontsize=18)
      plt.xlabel('Epochs', fontsize=18)
      val_acc_max = np.max(val_acc)
      val_acc_max_ind = np.argmax(val_acc)
      plt.axvline(x=val_acc_max_ind, color='g', linestyle='--', label='Highest validation accuracy')
      plt.title('Highest validation accuracy = %0.1f %%' % (val_acc_max*100), fontsize=16)
      ax.legend(loc='lower right', fontsize=16)
      #plt.ion()
      plt.show()

    #accuracy
    ind = np.argmax(val_acc)
    class_accuracy = val_confusion_matrix[:,:,ind]
    for ii in range(len(classes)):
        acc = val_confusion_matrix[ii,ii,ind] / np.sum(val_confusion_matrix[ii,:,ind])
        print(f'Accuracy of {str(classes[ii]).ljust(15)}: {acc*100:.01f}%')
    print('best average validation accuracy',val_acc[ind])



def save_feature_hook(module, input_ , output, outputpath, filestub):

  outname= os.path.join(outputpath, filestub )

  print(outname)
  for el in input_:
    print('shape',el.shape)
  #np.save(outname,  input_[0].clone().to('cpu').data.numpy())
  np.save(outname,  output[0].clone().to('cpu').data.numpy())


def hook_factory(outputpath, filestub):

    if not os.path.isdir(outputpath):
      os.makedirs(outputpath)

    # define the function with the right signature to be created
    def ahook(module, input_, output):
        # instantiate it by taking a parametrized function,
        # and fill the parameters
        # return the filled function
        return save_feature_hook(module, input_, output, outputpath = outputpath, filestub = filestub)

    # return the hook function as if it were a string
    return ahook

if __name__=='__main__':
    run()

   

import os,sys,numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

def run():
 
  outputpath= './tmp/'
  
  batchinds= np.arange(0,50)
  layernames=['cnn{:d}'.format(i) for i in range(1,7)]

  quantiles=[50,70,90]

  hists= [[] for i in range(len(layernames))] 
  valsperquantile = np.zeros(( len(quantiles), len(layernames) ))
  
  fig, axs = plt.subplots(1, len(quantiles))
  maxval=-1
  for li,layer in enumerate(layernames):
    
    for batch_idx in batchinds:
      filestub='batchindex{:d}'.format(batch_idx)+'_'+layer+'.npy'
      outname= os.path.join(outputpath, filestub )
      gs = np.load( outname)
      print(gs.shape)
      
      hists[li].append( np.abs(gs).flatten().tolist() )
      #hists[li].append( np.mean(np.abs(gs),axis=(2,3)).flatten().tolist() )
      
    vals = np.percentile(hists[li],quantiles)
    maxval = max (maxval, vals[-1])
    for i,q in enumerate(quantiles):
      valsperquantile[i,li]=vals[i]

  for i,q in enumerate(quantiles):
    axs[i].plot( [int(i+1) for i in range(len(layernames))], valsperquantile[i,:],'g+')
    axs[i].xaxis.set_major_locator(mticker.MultipleLocator(1))

    
  plt.show()
if __name__=='__main__':
  run()

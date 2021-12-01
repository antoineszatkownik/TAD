import tensorflow as tf
from tensorflow.keras import datasets, layers, models
from keras.models import Sequential
from keras.utils import np_utils
from sklearn.manifold import MDS 
import keras

import matplotlib.pyplot as plt

import pandas as pd
import numpy as np

from scipy.spatial import ConvexHull
from scipy import sparse, ndimage, misc, signal

import HiCtoolbox

import sys
import csv

print("Num GPUs Available: ", len(tf.config.list_physical_devices('GPU')))

def sld_windows(mat,arr_head,fsize,pos):
    dicx, dicy = [],[]
    if fsize%2 != 0:
        fsize+=1
    for i in pos:
        if i + 1 + fsize//2 in arr_head:
            label = 1
        else: label = 0
        dicx.append(mat[i:i+fsize,i:i+fsize])
        dicy.append(label)
    
    return (np.array(dicx), np.array(dicy))

def HiC(chrom = 1, R = 25000, alpha = 0.227):
    HiCfilename= r'dataforstudent\HiC\GM12878\25kb_resolution_intrachromosomal\chr{}_25kb.RAWobserved'.format(chrom)
    
    #Build matrix
    A=np.loadtxt(HiCfilename)
    A=np.int_(A)
    A=np.concatenate((A,np.transpose(np.array([A[:,1],A[:,0],A[:,2]]))), axis=0)#build array at pb resolution
    A = sparse.coo_matrix( (A[:,2], (A[:,0],A[:,1])))
    binned_map=HiCtoolbox.bin2d(A,R,R) #!become csr sparse array
    del A #keep space
    contact_map=HiCtoolbox.SCN(binned_map.copy()) 

    return np.asarray(contact_map)**alpha #now we are not sparse at all

def arr(csv, chrom):
    return csv[csv.chr1==chrom]

print("Read HIC")
Hic = HiC()
print("DONE!")
print()

class DataGenerator(tf.keras.utils.Sequence):
    'Generates data for Keras'
    def __init__(self,Hic,chrom=1, batch_size=32, dim=(32,32,32), n_channels=1,
                 n_classes=2, shuffle=True):
        'Initialization'
        self.dim = dim
        self.batch_size = batch_size
        self.n_channels = n_channels
        self.n_classes = n_classes
        self.shuffle = shuffle
        self.chrom = chrom
        self.hic = Hic
        
        tmp_rec = pd.read_csv("dataforstudent\\Arrowhead\\GSE63525_GM12878_primary+replicate_Arrowhead_domainlist.txt", sep="\t")
        rec1 = tmp_rec[tmp_rec.chr1 ==str(chrom)]                   #Chr choice
        arrec = np.array(rec1[['x1','x2']])          #get coords,only x1&x2 as x1=y1 and x2=y2
        arrec = (arrec//25000)*25000                 #As 25kb res
        self.rec = np.ceil(arrec/25000)                  #Treated as HiC with rescaled factor R
        
        self.list_IDs = [i for i in range(len(self.hic)-32)]
        self.on_epoch_end()

    def __len__(self):
        'Denotes the number of batches per epoch'
        return int(np.floor(len(self.list_IDs) / self.batch_size))

    def __getitem__(self, index):
        'Generate one batch of data'
        # Generate indexes of the batch
        indexes = self.indexes[index*self.batch_size:(index+1)*self.batch_size]

        # Find list of IDs
        list_IDs_temp = [self.list_IDs[k] for k in indexes]

        # Generate data
        X, y = self.__data_generation(list_IDs_temp)

        return X, y

    def on_epoch_end(self):
        'Updates indexes after each epoch'
        self.indexes = np.arange(len(self.list_IDs))
        if self.shuffle == True:
            np.random.shuffle(self.indexes)

    def __data_generation(self, list_IDs_temp):
        'Generates data containing batch_size samples' # X : (n_samples, *dim, n_channels)
        # Initialization
        #arrec = arr(self.rec,self.chrom)
        X = np.empty((self.batch_size, self.n_channels,*self.dim))
        y = np.empty((self.batch_size), dtype=int)
        for i, ID in enumerate(list_IDs_temp):
            X[i,], y[i] = sld_windows(self.hic,self.rec,32,[ID])
            
            
            #print(X.shape)
        return X, y


# Parameters
params = {'dim': (32,32),
          'batch_size': 32,
          'n_classes': 2,
          'n_channels': 1,
          'shuffle': True}


# Generators
training_generator = DataGenerator(Hic, **params)
#validation_generator = DataGenerator(Hic,**params)
"Generator created"

# Design model
#
model = Sequential()
model.add(layers.Conv2D(32, (3, 3), activation='relu',data_format='channels_first', input_shape=(1,32,32)))
model.add(layers.MaxPool2D((2, 2), data_format='channels_first' ))
model.add(layers.Dropout(0.25))

model.add(layers.Conv2D(64, (3, 3), activation='relu'))
model.add(layers.MaxPool2D((2, 2)))
model.add(layers.Dropout(0.25))

model.add(layers.Conv2D(128, (3, 3), activation='relu'))
model. add(layers.MaxPool2D(2,2))
model.add(layers.Dropout(0.25))

model.add(layers.Flatten())
model.add(layers.Dense(64,activation = 'relu'))
model.add(layers.Dropout(0.5))
model.add(layers.Dense(1, activation = 'sigmoid'))

model.compile(loss = keras.losses.binary_crossentropy, optimizer = 'adam', metrics=['accuracy'])

print()
print("Fit model:")
# Train model on dataset
model.fit(training_generator, validation_data=training_generator,epochs = 10, verbose =2)
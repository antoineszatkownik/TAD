import tensorflow as tf
from tensorflow.keras import  layers
from keras.models import Sequential
import keras

import numpy as np
import os

import pandas as pd

class DataGenerator(tf.keras.utils.Sequence):
    'Generates data for Keras'
    def __init__(self, list_IDs, batch_size=32, dim=(32,32,32), n_channels=1,
                 n_classes=2, shuffle=True):
        'Initialization'
        self.dim = dim
        self.batch_size = batch_size
        self.list_IDs = list_IDs
        self.n_channels = n_channels
        self.n_classes = n_classes
        self.shuffle = shuffle
        self.on_epoch_end()

    def __len__(self):
        'Denotes the number of batches per epoch'
        return int(np.floor(len(self.list_IDs) / self.batch_size))

    def __getitem__(self, index):
        'Generate one batch of data'
        # Generate indexes of the batch
        indexes = self.indexes[index*self.batch_size:(index+1)*self.batch_size]

        # Find list of IDs
        list_IDs_temp = [list(self.list_IDs.keys())[k] for k in indexes]

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
        X = np.empty((self.batch_size, self.n_channels, *self.dim))
        y = np.empty((self.batch_size), dtype=int)

        # Generate data
        for i, ID in enumerate(list_IDs_temp):
            # Store sample
            #print(ID)
            X[i,] = np.array(pd.read_csv('data/'+str(self.list_IDs[ID])+'/' + ID,
                                         header=None,
                                         sep=" "))

            # Store class
            y[i] = self.list_IDs[ID]

        return X, y

params = {'dim': (33,33),
          'batch_size': 32,
          'n_classes': 2,
          'n_channels': 1,
          'shuffle': True}

train_part = {}
val_part = {}

for root, dirs, files in os.walk("data", topdown=False):
    for name in files:
        l = name.split("_")
        if int(l[1]) < 16:
            train_part[name] = int(root[-1])
        
        elif int(l[1]) >= 16 & int(l[1]) <18:
            val_part[name] = int(root[-1])

training_generator = DataGenerator( train_part, **params)
validation_generator = DataGenerator(val_part, **params)

model = Sequential()
model.add(layers.Conv2D(32, (3, 3), activation='relu',data_format='channels_first', input_shape=(1,33,33)))
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

callback = tf.keras.callbacks.EarlyStopping(monitor='loss', patience=5)
print()
print("Fit model:")
# Train model on dataset
history = model.fit(training_generator, validation_data=training_generator,epochs = 100, verbose =2)

hist_df = pd.DataFrame(history.history) 
hist_csv_file = 'history.csv'
with open(hist_csv_file, mode='w') as f:
    hist_df.to_csv(f)

model.save("model.h5", save_format = "h5")
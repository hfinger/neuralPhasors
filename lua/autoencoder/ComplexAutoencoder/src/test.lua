require 'unsup';
require 'nn';
require 'gnuplot';
require 'Encoder';
require 'Decoder';
require 'image';
require 'optim';
require 'functions'
require 'randomkit'
require 'pl'


train = torch.load('/net/store/nbp/projects/phasesim/src_kstandvoss/lua/autoencoder/mnist.t7/train_32x32.t7', 'ascii')
  
  train = train.data


  traindata = torch.Tensor(train:size()[1],1,32,32)
  for i = 1,train:size()[1] do   
    traindata[i] = train[i]
  end
  traindata:div(255)
  print(traindata[1])
--print(torch.sqrt( torch.pow(autoencoder.output[1],2) + torch.pow(autoencoder.output[2],2)))











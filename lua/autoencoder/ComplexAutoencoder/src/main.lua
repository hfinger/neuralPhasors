package.path = package.path .. ";/net/store/nbp/projects/phasesim/src_kstandvoss/lua/autoencoder/ComplexAutoencoder/src/?.lua";
require 'nn';
require 'gnuplot';
require 'Encoder';
require 'Decoder';
require 'image';
require 'optim';
require 'functions';
require 'eval';
require 'randomkit';
require 'pl';
require 'mattorch';
local matFiles = require 'matFiles';
local start = require 'start'

local lfs = require"lfs";
loaded = false
function cuda()
    require 'cutorch';
    require 'cunn';
end
if pcall(cuda) then loaded = true else print('Cuda cannot be required') end
    
opt = lapp[[
   -f,--full          (default 5000)        use the full dataset or only n samples
   -o,--optimization  (default "SGD")       optimization: SGD | adadelta | adam
   -r,--learningRate  (default 0.05)        learning rate, for SGD only
   -b,--batchSize     (default 10)          batch size
   -m,--momentum      (default 0)           momentum, for SGD only
   -d,--coefL2        (default 0)           L2 penalty on the weights
   -h,--hidden        (default 50)          number of hidden units
   -k,--kernel        (default 7)           kernelsize
   -p,--phase                               training with/without phase 
   -n,--noise         (default 0)           noise level for denoising         
   -c,--criterion     (default MSE)         Loss Function: MSE | BCE | Complex
   -l,--coefL1        (default 1e-4)        L1 penalty on the weights
   -u,--cuda                                Use Cuda
   -w,--workdir       (default no)          working Directory
   -i,--fileId        (default no)          file identifier
   -s,--dataSet       (default LabelMe)     Dataset
   -e,--epochs        (default 1)           Number of Dataset Iterations
   -t,--stacked                             use second layer 
]]

steps = opt.full
workdir = opt.workdir
fileid = opt.fileId
epochs = opt.epochs or 1
algo = opt.optimization
--load Data
print('Load Data')

local traindata, testdata
local inpD
if opt.dataSet == 'LabelMe' then
  traindata = 'LabelMe'
  testdata = loadData('LabelMe',10,1)
  inpD = 6
  w = 400
  h = 300
  steps = steps - 20
  opt.full = steps
elseif opt.dataSet == 'MNIST' then
  traindata, testdata = loadData(false, 'mnist')
  traindata = createData(steps+100, traindata, 5)
  testdata = createData(100, testdata, 5)
  inpD = 1
  w = 100
  h = 100
elseif opt.dataSet == 'white' then
  traindata = createData(steps+100, 'white', 6)
  testdata = createData(100, 'white', 6)
  inpD = 1
  w = 32
  h = 32  
end

conf = {inputDim = inpD, width=w, height=h}

--[[
autoencoder, error = start.run(traindata, opt, opt.stacked , conf, epochs)
autoencoder:evaluate()
evaluate(autoencoder, opt.cuda, error, steps, workdir, i, testdata, epochs, opt.batchSize, opt.dataSet)
--matFiles.run(autoencoder, opt.workdir)

autoencoder:clearState()
if opt.stacked then
  name = 'stackedModel.net'
else 
  name = 'model.net'
end   
--torch.save(name,autoencoder)

--]]

set = opt.dataSet
--coefL1 = {0 , 1e-2, 1e-3, 1e-4, 1e-5, 1e-6}
--noise = {0, 0.1, 0.3, 0.5, 0.8}

kernel = {5, 7, 9 , 11}
coefL2 = {20, 50, 100}
rates = {0.0,0.1, 0.3, 0.5, 0.8}


i = 1
for k1, l1 in pairs(kernel) do
  for k2, l2 in pairs(coefL2) do
    for k3, lr in pairs(rates) do
      if i > 13 and i < 30 then
        collectgarbage()       
        opt = {full = steps, coefL2 = 0.001, coefL1 = 0, hidden = l2, cuda = true, noise = lr, batchSize = opt.batchSize, learningRate = 0.001, kernel = l1, criterion = 'MSE', workdir = i, optimization = opt.optimization}
        autoencoder, error = start.run(traindata, opt, opt.stacked, conf, epochs)
        autoencoder:evaluate()
        autoencoder:clearState()
        --evaluate(autoencoder, opt.cuda, error, steps, i .. '/', 1, testdata, params, epochs, opt.batchSize, set)
        matFiles.run(autoencoder, i)
        params = {l2 = l2, rate = lr, kernel = l1}
        torch.save('params'..i..'.dat',params)
      end  
      i = i + 1
    end
  end
end



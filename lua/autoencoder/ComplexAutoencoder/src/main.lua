package.path = package.path .. ";/net/store/nbp/projects/phasesim/src_kstandvoss/lua/autoencoder/ComplexAutoencoder/src/?.lua";
require 'unsup';
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
   -o,--optimization  (default "SGD")       optimization: SGD | ADAGRAD
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
epochs = opt.epochs
--load Data
print('Load Data')

local traindata, testdata
local inpD
if opt.dataSet == 'LabelMe' then
  traindataPos = torch.Tensor(3,300,400)
  traindataNeg = torch.Tensor(3,300,400)
  splitdata = torch.Tensor(steps,6,300,400):zero()
  print('load')
  --f = 1254
  f = 1010735
  for i=1,steps do
      f = f + 1
      function get()
        --test = mattorch.load('/net/store/nbp/projects/phasesim/workdir/kstandvoss/labelMeWhite/'..i..'/static_newyork_city_urban/IMG_'..f..'.jpg/act1.mat')   
        test = mattorch.load('/net/store/nbp/projects/phasesim/workdir/20130726_Paper/Autoencoder/labelMeWhite/05june05_static_street_boston/p' .. f..'.jpg/act1.mat')
        traindataPos = test['act']:transpose(2,3)
        traindataNeg:copy(traindataPos)
        traindataPos[traindataPos:le(0)] = 0 
        traindataNeg[traindataNeg:ge(0)] = 0
        splitdata[i]:sub(1,3):add(100,traindataPos)
        splitdata[i]:sub(4,6):add(100,traindataNeg)
      end  
      if pcall(get) then else f = f + 1; i = i -1; end
  end   

  traindata = splitdata
  testdata = traindata
  inpD = 6
  w = 300
  h = 400
elseif opt.dataSet == 'MNIST' then
  traindata, testdata = loadData(false, 'mnist')
  traindata = createData(steps+100, traindata, 2)
  testdata = createData(steps+100, testdata, 5)
  inpD = 1
  w = 100
  h = 100
elseif opt.dataSet == 'white' then
  traindata = createData(steps+100, 'white', 2)
  testdata = createData(steps+100, 'white', 4)
  inpD = 1
  w = 32
  h = 32  
end

conf = {inputDim = inpD, width=w, height=h}
autoencoder, error = start.run(traindata, testdata, opt, opt.stacked , conf, epochs)
autoencoder:evaluate()
if opt.stacked then epochs = epochs * 2 end
evaluate(autoencoder, opt.cuda, error, steps, workdir, fileid, testdata, params, epochs, opt.batchSize, opt.dataSet)

--[[
coefL1 = {0.1, 0.3, 0.5, 0.7, 0.9}
coefL2 = {10, 20, 25}
rates = {5, 7, 9}
i = 1
for k1, l1 in pairs(coefL1) do
  for k2, l2 in pairs(coefL2) do
    for k3, lr in pairs(rates) do     
      opt = {full = steps, coefL2 = 1e-6, coefL1 = 1e-3, hidden = 30, cuda = true, noise = l1, batchSize = l2, learningRate = 0.1, kernel = lr, criterion = 'MSE', workdir = i}
      start.run(traindata, testdata, opt)
      i = i + 1
    end
  end
end
]]--



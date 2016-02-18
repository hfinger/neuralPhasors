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

local lfs = require"lfs";
loaded = false
function cuda()
    require 'cutorch';
    require 'cunn';
end
if pcall(cuda) then loaded = true else print('Cuda cannot be required') end
    

opt = lapp[[
   -f,--full          (default 5000)        use the full dataset or only n samples
   -o,--optimization  (default "SGD")       optimization: SGD | LBFGS 
   -r,--learningRate  (default 0.05)        learning rate, for SGD only
   -b,--batchSize     (default 10)          batch size
   -m,--momentum      (default 0)           momentum, for SGD only
   -d,--weightDecay   (default 5e-5)        weight decay for SGD
   -h,--hidden        (default 50)          number of hidden units
   -k,--kernel        (default 7)           kernelsize
   -p,--phase                               training with/without phase 
   -n,--noise         (default 0)           noise level for denoising         
   -c,--criterion     (default MSE)         Loss Function: MSE | BCE | Complex
   -l,--coefL1        (default 1e-4)        L1 penalty on the weights
   -u,--cuda                                Use Cuda
   -w,--workdir       (default no)          working Directory
   -i,--fileId        (default no)          file identifier
   --coefL2           (default 0)           L2 penalty on the weights
]]

local steps = opt.full
local batchSize = opt.batchSize
local l1 = opt.coefL1
local noise = opt.noise


local fileid
if opt.fileid == 'no' then
    fileid = false
else    
    fileid = opt.fileId
end
local workdir
if opt.workdir == 'no' then
    workdir = false
else    
    workdir = opt.workdir
end
--load Data
--traindata, testdata = loadData()
local traindata = createData(steps+100):div(255)
local testdata = createData(steps+100):div(255)

local n_input = 32*32
local config = {learningRate = opt.learningRate,
            weightDecay = opt.weightDecay,
            momentum = opt.momentum,
            learningRateDecay = 5e-5,
            nesterov = true,
            dampening = 0}
         
--initialize Encode
local enc = nn.Encoder(1,opt.hidden,opt.kernel)
local dec = nn.Decoder(opt.hidden,1,opt.kernel)



local autoencoder = nn.Sequential()
autoencoder:add(enc)
autoencoder:add(dec)

    
--set Criterion 
local criterion = nn.ParallelCriterion()

if opt.criterion == 'MSE' then
    xcrit = nn.MSECriterion() --compare x_in and x_out
    ycrit = nn.MSECriterion() --compare y_in and y_out
elseif opt.criterion == 'BCE' then
    xcrit = nn.BCECriterion() --compare x_in and x_out
    ycrit = nn.BCECriterion() --compare y_in and y_out
else
    xcrit = nn.ComplexCrit() --compare x_in and x_out
    ycrit = nn.ComplexCrit() --compare y_in and y_out
end
criterion:add(xcrit,0.5)
criterion:add(ycrit,0.5)


if loaded and opt.cuda then
    print('Use cuda')
    enc:cuda()
    dec:cuda()
    autoencoder:cuda()
    criterion:cuda()
    traindata:cuda()
end
local parameters, gradParameters = autoencoder:getParameters()
timer = torch.Timer()

print("start training")
error = trainModel(steps, batchSize, traindata, autoencoder, criterion, parameters, gradParameters, config, opt.phase, noise, opt.cuda, l1)
name = 'model' .. opt.full .. '.net'
print('finished training after ' .. timer:time().real)
evaluate(autoencoder, opt.cuda, false, steps, workdir, fileid)

--[[

-- stacked
enc1 = enc:clone()
enc1:evaluate()
enc2 = nn.Encoder(50,10,7,1e-5) --eventuell keine Convultion sondern Linear und eventuell pooling
dec2 = nn.Decoder(10,1,7)
autoencoder2 = nn.Sequential()
autoencoder2:add(enc1)
autoencoder2:add(enc2)
autoencoder2:add(dec2)

if loaded and opt.cuda then
    enc2:cuda()
    enc2.encoder:cuda()
    dec2:cuda()
    dec2.decoder:cuda()
    autoencoder2:cuda()
    testdata:cuda()
end

parameters, gradParameters = autoencoder2:getParameters()
print("start training")
steps = opt.full
batchSize = opt.batchSize
trainModel(steps, batchSize, testdata, autoencoder2, criterion, parameters, gradParameters, config, opt.phase, opt.noise, opt.cuda)
name = 'stackedModel' .. opt.full .. '.net'
print('finished training after ' .. timer:time().real)
evaluate(autoencoder2, opt.cuda, false, steps)

if opt.cuda then 
    enc = enc:float()
    enc.encoder = cudnnNetToCpu(enc.encoder)
    dec = dec:float()
    dec.decoder = cudnnNetToCpu(dec.decoder)
    ae = autoencoder:float()
end   


--]]

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

--load Data
traindata, testdata = loadData{normalizeMean = false} 

opt = lapp[[
   -f,--full          (default 5000)        use the full dataset or only n samples
   -o,--optimization  (default "SGD")       optimization: SGD | LBFGS 
   -r,--learningRate  (default 0.05)        learning rate, for SGD only
   -b,--batchSize     (default 10)          batch size
   -m,--momentum      (default 0)           momentum, for SGD only
   -d,--weightDecay   (default 5e-5)        weight decay for SGD
   -h,--hidden        (default 50)          number of hidden units
   -k,--kernel        (default 7)           kernelsize
   -p,--phase         (default false)       training with/without phase
   -n,--noise         (default 0.1)         noise level for denoising         
   -c, --criterion    (default MSE)         Loss Function: MSE | BCE | Complex
   -w
   --coefL1           (default 0)           L1 penalty on the weights
   --coefL2           (default 0)           L2 penalty on the weights
]]


n_input = 32*32
config = {learningRate = opt.learningRate,
         weightDecay = opt.weightDecay,
         momentum = opt.momentum,
         learningRateDecay = 5e-7}
         
--initialize Encode
enc = nn.Encoder(1,opt.hidden,opt.kernel)
dec = nn.Decoder(opt.hidden,1,opt.kernel)


autoencoder = nn.Sequential()
autoencoder:add(enc)
autoencoder:add(dec)

--set Criterion 
criterion = nn.ParallelCriterion()

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

parameters, gradParameters = autoencoder:getParameters()

timer = torch.Timer()
print("start training")
steps = opt.full
batchSize = opt.batchSize
trainModel(steps, batchSize, traindata, autoencoder, criterion, parameters, gradParameters, config, opt.phase, opt.noise)
name = 'model' .. opt.full .. '.net'
print('finished training after ' .. timer:time().real .. ' - save model ' .. name)
torch.save(name,autoencoder)




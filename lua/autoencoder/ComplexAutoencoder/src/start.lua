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

local R = {}

local function run(traindata, testdata, opt, stacked, conf, epochs)
    steps = opt.full
    batchSize = opt.batchSize
    l1 = opt.coefL1
    l2 = opt.coefL2
    
   
    noise = opt.noise
    
    inputDim = conf.inputDim
    width = conf.width
    height = conf.height
     
    fileid = null
    if opt.fileId == 'no' then
        fileid = false
    else    
        fileid = opt.fileId
    end
    workdir = null
    if opt.workdir == 'no' then
        workdir = false
    else    
        workdir = '/net/store/nbp/projects/phasesim/src_kstandvoss/lua/autoencoder/ComplexAutoencoder/Data/'
        workdir = workdir .. opt.workdir .. '/'
        os.execute("mkdir " .. workdir)
    end
    local config = {learningRate = opt.learningRate, momentumDecay=0.99, updateDecay=0.99}
    --local config = {learningRate = opt.learningRate, learnungRateDecay=1e-4, momentum=0.9, nesterov=true, dampening=0}
    --initialize Encode
    local enc = nn.Encoder(inputDim,opt.hidden,opt.kernel,false, width, height)
    local dec = nn.Encoder(opt.hidden,inputDim,opt.kernel,false, width, height)
    
    
    
    local autoencoder = nn.Sequential()
    autoencoder:add(enc)
    autoencoder:add(dec)
    autoencoder:training()
    
    --set Criterion 
    local criterion = nn.ParallelCriterion()
    
    if opt.criterion == 'MSE' then
        xcrit = nn.MSECriterion() --compare x_in and x_out
        ycrit = nn.MSECriterion() --compare y_in and y_out
    elseif opt.criterion == 'BCE' then
        xcrit = nn.SmoothL1Criterion() --compare x_in and x_out
        ycrit = nn.SmoothL1Criterion() --compare y_in and y_out
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
    autoencoder:get(1).convolution2:share(autoencoder:get(1).convolution1,'weight','bias','gradWeight','gradBias')
    autoencoder:get(1).convolution3:share(autoencoder:get(1).convolution1,'weight','bias','gradWeight','gradBias')
    autoencoder:get(2).convolution2:share(autoencoder:get(2).convolution1,'weight','bias','gradWeight','gradBias')
    autoencoder:get(2).convolution3:share(autoencoder:get(2).convolution1,'weight','bias','gradWeight','gradBias')
    
    lightModel = autoencoder:clone('weight','bias','running_mean','running_std')
    
    timer = torch.Timer()
    
    errors = torch.Tensor(epochs, steps/batchSize)       
    print("start training")
    for i = 1,epochs do
      n = noise
      error = trainModel(steps, batchSize, traindata, autoencoder, criterion, parameters, gradParameters, config, opt.phase, noise, opt.cuda, l1, l2)
      noise = n
      errors[i] = error
    end
    errors = errors:view(-1)
    name = 'model' .. opt.full .. '.net'
    print('finished training after ' .. timer:time().real)
    torch.save(name, lightModel)
    params = {
        d = opt.coefL2,
        l = opt.coefL1,
        r = opt.learningRate
    }
    
    if not stacked then
        return autoencoder, errors
    else
        dec = autoencoder:get(2):clone()
        enc = autoencoder:get(1):clone()

        enc2 = nn.Encoder(4 * opt.hidden,10,7,false, 400, 300) --eventuell keine Convultion sondern Linear und eventuell pooling
        dec2 = nn.Encoder(10, 4 * opt.hidden,7,false, 400, 300)
        autoencoder2 = nn.Sequential()
        autoencoder2:add(enc2)
        autoencoder2:add(dec2)
        autoencoder2:training()
        
        
        stackedEnc = nn.Sequential()
        stackedEncPar = nn.ParallelTable()
        stackedEncSeq = nn.Sequential()
        stackedEncSeq:add(nn.Reshape(opt.hidden, 150, 2, 200, 2))
        stackedEncSeq:add(nn.Transpose({2,3},{3,5},{4,5}))
        stackedEncSeq:add(nn.Reshape(opt.hidden*4,150,200))
        stackedEncPar:add(stackedEncSeq)
        stackedEncPar:add(stackedEncSeq)
        stackedEnc:add(enc)
        stackedEnc:add(stackedEncPar)
        stackedEnc = stackedEnc:cuda()
        
        stackedDec = nn.Sequential()
        stackedDecPar = nn.ParallelTable()
        stackedDecSeq = nn.Sequential()
        stackedDecSeq:add(nn.Reshape(opt.hidden,2,2,150,200))
        stackedDecSeq:add(nn.Transpose({4,5},{3,5},{2,3}))
        stackedDecSeq:add(nn.Reshape(opt.hidden, 300, 400))
        stackedDecPar:add(stackedDecSeq)
        stackedDecPar:add(stackedDecSeq)
        stackedDec:add(stackedDecPar)
        stackedDec:add(dec)
        stackedDec = stackedDec:cuda()
        
        if loaded and opt.cuda then
            enc:cuda()
            enc2:cuda()
            enc2.encoder:cuda()
            dec2:cuda()
            dec2.encoder:cuda()
            autoencoder2:cuda()
            criterion:cuda()
        end

        local config2 = {learningRate = opt.learningRate, momentumDecay= 0.99, updateDecay=0.99}

        parameters, gradParameters = autoencoder2:getParameters()
        autoencoder2:get(1).convolution2:share(autoencoder2:get(1).convolution1,'weight','bias','gradWeight','gradBias')
        autoencoder2:get(1).convolution3:share(autoencoder2:get(1).convolution1,'weight','bias','gradWeight','gradBias')
        autoencoder2:get(2).convolution2:share(autoencoder2:get(2).convolution1,'weight','bias','gradWeight','gradBias')
        autoencoder2:get(2).convolution3:share(autoencoder2:get(2).convolution1,'weight','bias','gradWeight','gradBias')
      
        print("start training")
        steps = opt.full
        batchSize = opt.batchSize
        timer2 = torch.Timer()
        
        errors = torch.Tensor(epochs * 2, steps/batchSize)
        for i = 1,2*epochs do
          n = noise
          error = trainModel(steps, batchSize, traindata, autoencoder2, criterion, parameters, gradParameters, config2, opt.phase, 0.3, opt.cuda, 0, 1e-5, stackedEnc)
          noise = n
          errors[i] = error
        end  
        errors = errors:view(-1)
        name = 'stackedModel' .. opt.full .. '.net'
        print('finished training after ' .. timer2:time().real)

        net = nn:Sequential()
        net:add(stackedEnc)
        net:add(enc2)
        net:add(dec2)
        net:add(stackedDec)

        --torch.save(name, net)
        return net, errors
    end    
end

R.run = run

return R

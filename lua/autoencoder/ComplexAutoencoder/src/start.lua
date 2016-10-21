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
    print(noise)
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
    local enc = nn.Encoder(inputDim,opt.hidden,opt.kernel,false, width, height, 'forward')
    local dec = nn.Encoder(opt.hidden,inputDim,opt.kernel,false, width, height, 'backward')
    
    
    
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
        traindata = traindata:cuda()
        testdata = testdata:cuda()
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
        m = opt.hidden
        enc2 = nn.Encoder( m , 10 ,7,false, width, height) --eventuell keine Convultion sondern Linear und eventuell pooling
        dec2 = nn.Encoder(10, m, 7,false, width, height)
        
        con1 = nn.SpatialConvolution(opt.hidden, m, 1, 1):cuda() --nn.Identity() 
        con2 = nn.SpatialConvolution(opt.hidden, m, 1, 1):cuda()
        con2:share(con1,'weight','bias','gradWeight','gradBias')
        

        autoencoder2 = nn.Sequential()
        --autoencoder2:add(nn.ParallelTable():add(con1):add(con2)) 
        autoencoder2:add(enc2)
        autoencoder2:add(dec2)
                
        con3 = nn.SpatialConvolution(m, opt.hidden, 1, 1):cuda()
        con4 = nn.SpatialConvolution(m, opt.hidden, 1, 1):cuda()
        con4:share(con3,'weight','bias','gradWeight','gradBias')
        --autoencoder2:add(nn.ParallelTable():add(con3):add(con4)) 
        
        autoencoder2:training()

        --enc.trans:activate()
        if loaded and opt.cuda then
            enc:cuda()
            dec:cuda()
            enc2:cuda()
            enc2.encoder:cuda()
            dec2:cuda()
            dec2:cuda()
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
          error = trainModel(steps, batchSize, traindata, autoencoder2, criterion, parameters, gradParameters, config2, opt.phase, 0.6, opt.cuda, 0, 1e-5, enc)
          noise = n
          errors[i] = error
        end  
        errors = errors:view(-1)
        name = 'stackedModel' .. opt.full .. '.net'
        print('finished training after ' .. timer2:time().real)
        
        --dec2.trans:activate() 
        
        net = nn:Sequential()
        net:add(enc)
        net:add(autoencoder2)
        net:add(dec)
       -- print(dec:listModules())
        if loaded and opt.cuda then
            enc:cuda()
            enc.encoder:cuda()
            enc2:cuda()
            enc2.encoder:cuda()
            dec2.encoder:cuda()
            dec2:cuda()
            dec.encoder:cuda()
            dec:cuda()
            dec.encoder:cuda()
            net:cuda()
        end
        --torch.save(name, net)
        return net, errors
    end    
end
--]]
R.run = run

return R

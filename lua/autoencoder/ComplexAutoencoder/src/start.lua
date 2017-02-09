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

local function run(traindata, opt, stacked, conf, epochs)
    steps = opt.full
    batchSize = opt.batchSize
    l1 = opt.coefL1
    l2 = opt.coefL2
    algo = opt.optimization
   
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
    --local config = {learningRate = opt.learningRate, learnungRateDecay=1e-4, momentum=0.9, nesterov=true, dampening=0}
    --initialize Encode
    pool1 = nn.SpatialAveragePooling(2,2,2,2)--,padW=3,padH=3}
    pool2 = nn.SpatialAveragePooling(2,2,2,2)--,padW=3,padH=3}
          
    down = nn.ParallelTable():add(pool1):add(pool2):cuda()
    up = nn.ParallelTable():add(nn.SpatialUpSamplingNearest(2)):add(nn.SpatialUpSamplingNearest(2)):cuda()
    if stacked then 
      autoencoder = torch.load('/net/store/nbp/projects/phasesim/src_kstandvoss/lua/autoencoder/ComplexAutoencoder/src/model.net')   
      getEnc = 1
      getDec = 4
    else

      local enc = nn.Encoder(inputDim,opt.hidden,opt.kernel,false, width, height)
      local dec = nn.Encoder(opt.hidden,inputDim,opt.kernel,false, width, height)
      
      getEnc = 1
      getDec = 4
      
      autoencoder = nn.Sequential()
      autoencoder:add(enc)
      
      autoencoder:add(down)
      autoencoder:add(up)
      --end  
      autoencoder:add(dec)
      autoencoder:training()
      
      
      if loaded and opt.cuda then
          print('Use cuda')
          enc:cuda()
          dec:cuda()
          autoencoder:cuda()
          criterion:cuda()
          --traindata = traindata:cuda()
      end
      
      
      local parameters, gradParameters = autoencoder:getParameters()
      autoencoder:get(getEnc).convolution2:share(autoencoder:get(getEnc).convolution1,'weight','bias','gradWeight','gradBias')
      autoencoder:get(getEnc).convolution3:share(autoencoder:get(getEnc).convolution1,'weight','bias','gradWeight','gradBias')
      autoencoder:get(getDec).convolution2:share(autoencoder:get(getDec).convolution1,'weight','bias','gradWeight','gradBias')
      autoencoder:get(getDec).convolution3:share(autoencoder:get(getDec).convolution1,'weight','bias','gradWeight','gradBias')
      
      --lightModel = autoencoder:clone('weight','bias','running_mean','running_std')
      
      timer = torch.Timer()
      
      errors = torch.Tensor(epochs, steps/batchSize)       
      print("start training")
      for i = 1,epochs do
        print('epoch: ' .. i)
        n = noise
        error = trainModel(steps, batchSize, traindata, autoencoder, criterion, parameters, gradParameters, config, opt.phase, noise, opt.cuda, l1, l2, getEnc, getDec, algo)
        noise = n
        errors[i] = error
      end
      errors = errors:view(-1)
      name = 'model' .. opt.full .. '.net'
      print('finished training after ' .. timer:time().real)
      --torch.save(name, lightModel)
      params = {
          d = opt.coefL2,
          l = opt.coefL1,
          r = opt.learningRate
      }
    end  
    if not stacked then
        return autoencoder, errors
    else       
        dec = autoencoder:get(getDec):clone()
        enc = autoencoder:get(getEnc):clone()
        m = opt.hidden
        enc2 = nn.Encoder( m , 20 ,opt.kernel,false, width, height) --eventuell keine Convultion sondern Linear und eventuell pooling
        dec2 = nn.Encoder(20, m, opt.kernel,false, width, height)
        
                
        autoencoder2 = nn.Sequential()
        autoencoder2:add(enc2)
        autoencoder2:add(dec2)

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
        getEnc = 1
        getDec = 2
        local config2 = {learningRate = opt.learningRate, momentumDecay= 0.99, updateDecay=0.99}
        
        parameters, gradParameters = autoencoder2:getParameters()
        autoencoder2:get(getEnc).convolution2:share(autoencoder2:get(getEnc).convolution1,'weight','bias','gradWeight','gradBias')
        autoencoder2:get(getEnc).convolution3:share(autoencoder2:get(getEnc).convolution1,'weight','bias','gradWeight','gradBias')
        autoencoder2:get(getDec).convolution2:share(autoencoder2:get(getDec).convolution1,'weight','bias','gradWeight','gradBias')
        autoencoder2:get(getDec).convolution3:share(autoencoder2:get(getDec).convolution1,'weight','bias','gradWeight','gradBias')
      
        print("start training")
        steps = opt.full
        batchSize = opt.batchSize
        timer2 = torch.Timer()
       

        errors = torch.Tensor(epochs, steps/batchSize)
        for i = 1,epochs do
          print('epoch: ' .. i)
          n = noise
          error = trainModel(steps, batchSize, traindata, autoencoder2, criterion, parameters, gradParameters, config2, opt.phase, 0, opt.cuda, l1, l2, getEnc, getDec, algo, nn.Sequential():add(enc):add(down))
          noise = n
          errors[i] = error
        end  
        errors = errors:view(-1)
        name = 'stackedModel' .. opt.full .. '.net'
        print('finished training after ' .. timer2:time().real)
        
        --dec2.trans:activate() 
        
        net = nn:Sequential()
        net:add(enc)
        net:add(down)
        net:add(autoencoder2)
        net:add(up)
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
            net:cuda()
        end
        --torch.save(name, net)
        return net, errors
    end    
end
--]]
R.run = run

return R

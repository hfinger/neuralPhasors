
require 'rmsprop'
function loadData(normalizeMean, set)
  
  workdir = '/net/store/nbp/projects/phasesim/src_kstandvoss/lua/autoencoder'
  if set == 'mnist' then
      train = torch.load(workdir .. '/mnist.t7/train_32x32.t7', 'ascii')
      test = torch.load(workdir .. '/mnist.t7/test_32x32.t7', 'ascii')  
  elseif set == 'LabelMe' then
      train = image.load('labelMe.dat')
  end
  train = train.data
  test = test.data

  traindata = torch.Tensor(train:size()[1],1,32,32)
  for i = 1,train:size()[1] do   
    traindata[i] = train[i]
  end
  traindata:div(255)

  
  testdata = torch.Tensor(test:size()[1],1,32,32)
  for i = 1,test:size()[1] do   
    testdata[i] = test[i]
  end
  testdata:div(255)

  
  if normalizeMean then
    traindata:add(-traindata:mean())
    testdata:add(-testdata:mean())
  end

  return traindata, testdata
end


--TODO store results during training
function trainModel(steps, batchsize, data, model, criterion, parameters, gradParameters, config, usePhase, noiseLevel, cuda, l1, l2, stacked)
    local epoch = epoch or 1
    local l2error = torch.Tensor(steps/batchsize)
    local k = 1
    local shuffle = torch.randperm((#data)[1]):type('torch.LongTensor')
    data = data:index(1,shuffle)
    for i = 1,steps,batchsize do
        --data = createData(batchsize+1, traindata, 2)
        print(i / steps * 100 ..'%')
        local feval = function(x)
            if x ~= parameters then
                parameters:copy(x)
            end

            -- reset gradients
            gradParameters:zero()
            local f = 0 
            for j = 0,batchsize-1 do       
               
                local activity = data[i+j]
                local phase = torch.Tensor(#activity):zero()
                if usePhase then 
                    phase = (torch.rand(#activity)*2*math.pi)-math.pi
                end
                local input = {}
                local target = {}
                if stacked then
                    input = {torch.cmul(activity,torch.cos(phase)),torch.cmul(activity,torch.sin(phase))}    
                    if cuda then
                        activity = activity:cuda()
                        input[1] = input[1]:cuda()
                        input[2] = input[2]:cuda()
                    end   
                    --print(stacked)
                    stacked:forward(input)
                    activity = torch.sqrt(torch.pow(stacked.output[1],2) + torch.pow(stacked.output[2],2))
                    --phase = torch.atan2(stacked.output[2],stacked.output[1])  
                    phase = torch.Tensor(#activity):zero()
                end    
                
                if noiseLevel then
                    --noise = torch.Tensor(#activity):bernoulli(1-noiseLevel)
                    noise = randomkit.normal(torch.Tensor(#activity),0,noiseLevel)
                    if cuda then
                        noise = noise:cuda()
                        activity = activity:cuda()
                        phase = phase:cuda()
                    end   
                    --corrupted = torch.cmul(activity,noise)
                    corrupted = activity + noise
                    input = {torch.cmul(corrupted,torch.cos(phase)),torch.cmul(corrupted,torch.sin(phase))}
                    target = {torch.cmul(activity,torch.cos(phase)),torch.cmul(activity,torch.sin(phase))}
                else
                    input = {torch.cmul(activity,torch.cos(phase)),torch.cmul(activity,torch.sin(phase))}
                    target = input
                end
                
                if cuda then
                    activity = activity:cuda()
                    input[1] = input[1]:cuda()
                    input[2] = input[2]:cuda()
                    target[1] = target[1]:cuda()
                    target[2] = target[2]:cuda()
                end             
                local output = model:forward(input)

                if cuda then
                    output[1] = output[1]:cuda()
                    output[2] = output[2]:cuda()
                end    
                local err = criterion:forward(output,target) 
                local df_dw = criterion:backward(output,target)               
                local a_out = torch.sqrt(torch.pow(output[1],2) + torch.pow(output[2],2))
                f = f + err 
                
                model:backward(input,df_dw)
                 
                if l1 ~= 0 or l2 ~= 0 then
                  -- locals:
                  local norm,sign= torch.norm,torch.sign
                  -- Loss:
                  weightsEnc = model:get(1).convolution1.weight
                  weightsDec = model:get(2).convolution1.weight
                  f = f + l1 * norm(weightsEnc,1) + l1 * norm(weightsDec,1) 
                  f = f + l2 * norm(weightsEnc,2)^2/2 + l2 * norm(weightsDec,2)^2/2
                  -- Gradients:
                  wSize = weightsEnc:view(-1):size(1)
                  bSize = model:get(1).convolution1.bias:view(-1):size(1)
                  gradParameters:sub(1,wSize):add(sign(weightsEnc:view(-1)):mul(l1) + weightsEnc:view(-1):clone():mul(l2))
                  gradParameters:sub(wSize+bSize+1,2*wSize+bSize):add(sign(weightsDec:view(-1)):mul(l1) + weightsDec:clone():view(-1):mul(l2))    
                end
                
            end
            
            if not stacked and (epoch == 1 or epoch==5) then
              display(model) 
            end 
            gradParameters:div(batchsize)
            f = f/batchsize
            l2error[k] = f
            k = k+1
            return f, gradParameters
        end 
        optim.adadelta(feval, parameters, config)

    end 
    epoch = epoch + 1 
    return l2error
end

iter = 0
function display(model)
   iter = iter or 0
   require 'image'
   if iter % 10 == 0 then
     local conv = model:get(1).convolution1.weight:clone()
     local weights
     if conv:size(2) > 1 then
      weights = conv:sub(1,10,1,3) - conv:sub(1,10,4,6)
      else 
        weights = conv
      end  
      local test = image.toDisplayTensor{input = weights, zoom=4, padding = 2, scaleeach=true}
      image.save('/net/store/nbp/projects/phasesim/src_kstandvoss/lua/autoencoder/ComplexAutoencoder/Data/weights' .. iter  .. '.png', test)
    --mat = weights:view(weights:size(1), -1):double()
    --mattorch.save(workdir .. '/ComplexAutoencoder/Data/MatlabData/weights' .. iter .. '.mat', mat)     
   end
   iter = iter + 1
end

function trainStripes(number)
    white = torch.Tensor(1,32,32):fill(0)
    which = torch.random(0,100)
    n = torch.random(1,number)
    if which < 50 then 
        for i=1,n do
            row = torch.random(1,32)
            white[{1,row,{}}]:fill(1)
        end
    else
        for i=1,n do
            col = torch.random(1,32)
            white[{1,{},col}]:fill(1)
        end
    end    
    return white 
end

function whiteStripes(number)
    white = torch.Tensor(1,32,32):fill(0)
    n = torch.random(0,number)
    for i=1,n do
        row = torch.random(1,32)
        white[{1,row,{}}]:fill(1)
    end
 
    for i=1,number-n do
        col = torch.random(1,32)
        white[{1,{},col}]:fill(1)
    end
    return white 
end    

function createData(amount, set, number)
    if set == 'white' then
        data = torch.Tensor(amount,1,32,32)
        for i=1,amount do
            data[i] = whiteStripes(number)
        end
    else
        data = torch.Tensor(amount,1,100,100)
        for i=1,amount do
            data[i] = mnist(number,set)
        end
    end        
    return data
end   

function mnist(number,data)
    big = torch.Tensor(1,100,100):zero()
    n = torch.random(1,number)
    for i = 1,n do
        exp = torch.random(1,data:size()[1])
        stim = data[exp]
        stim = stim[{{1},{2,29},{2,29}}] --crop mnist numbers to get rid of zero padding
        x = torch.random(1,72)
        y = torch.random(1,72)
        big[{{1},{x,x+27}, {y,y+27}}] = torch.cmax(big[{{1},{x,x+27}, {y,y+27}}], stim)
    end
    return big
end    



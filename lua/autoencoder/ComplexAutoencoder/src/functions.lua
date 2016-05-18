
require 'rmsprop'
workdir = '/net/store/nbp/projects/phasesim/src_kstandvoss/lua/autoencoder'
function loadData(normalizeMean)
  
  train = torch.load(workdir .. '/mnist.t7/train_32x32.t7', 'ascii')
  test = torch.load(workdir .. '/mnist.t7/test_32x32.t7', 'ascii')
  
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
function trainModel(steps, batchsize, data, model, criterion, parameters, gradParameters, config, usePhase, noiseLevel, cuda, l1, stacked)
    l2error = torch.Tensor(steps/10 + 3)
    k = 1
    for i = 1,steps,batchsize do
    
        local feval = function(x)
            if x ~= parameters then
                parameters:copy(x)
            end

            -- reset gradients
            gradParameters:zero()
            local f = 0 
            for j = 0,batchsize-1 do        
                local activity = data[i+j]
                local phase = torch.Tensor(1,activity:size()[2],activity:size()[3]):zero()
                if usePhase then 
                    phase = (torch.rand(1,activity:size()[2],activity:size()[3])*2*math.pi)-math.pi
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
                    stacked:forward(input)
                    activity = torch.sqrt(torch.pow(stacked.output[1],2) + torch.pow(stacked.output[2],2))
                    phase = torch.atan2(stacked.output[2],stacked.output[1])  
                end    
                if noiseLevel then
                    noise = torch.Tensor(#activity):bernoulli(1-noiseLevel)--randomkit.normal(torch.Tensor(#activity),0,noiseLevel)
                    if cuda then
                        noise = noise:cuda()
                        activity = activity:cuda()
                        phase = phase:cuda()
                    end    
                    corrupted = torch.cmul(activity,noise)
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
                --local a_out = torch.sqrt(torch.pow(output[1],2) + torch.pow(output[2],2))
                if not stacked then
                    display(input[1],model)
                end 
--                if (i+j) % 10 == 0 and not stacked  then
--                    l2error[k] = torch.norm(activity-a_out)
--                   k = k + 1
--               end
                f = f + err +  l1 * torch.norm(x,1)
                model:backward(input,df_dw)
            end
            gradParameters:div(batchsize)
            f = f/batchsize
            return f, gradParameters
        end 
        rmsprop(feval, parameters, config)

    end 
    return l2error
end

iter = 0
function display(input,model)
   iter = iter or 0
   require 'image'
   if iter % 1000 == 0 then
    gnuplot.pngfigure(workdir .. '/ComplexAutoencoder/Data/weights' .. iter .. '.png')
    local weights = model:get(1):getWeights()
    local test = image.toDisplayTensor{input = weights, padding = 2}
    gnuplot.imagesc(test:view(test:size(2), test:size(3)))   
    gnuplot.plotflush()
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


function replaceModules(net, orig_class_name, replacer)
  local nodes, container_nodes = net:findModules(orig_class_name)
  for i = 1, #nodes do
    for j = 1, #(container_nodes[i].modules) do
      if container_nodes[i].modules[j] == nodes[i] then
        local orig_mod = container_nodes[i].modules[j]
        container_nodes[i].modules[j] = replacer(orig_mod)
      end
    end
  end
end

function cudnnNetToCpu(net)
  local net_cpu = net:clone():float()

  replaceModules(net_cpu, 'nn.FlattenTable', 
    function(orig_mod)
      local cpu_mod = nn.FlattenTable()
      return cpu_mod
    end)
  
  replaceModules(net_cpu, 'nn.ParallelTable', 
    function(orig_mod)
      local cpu_mod = nn.ParallelTable()
      return cpu_mod
    end)
    
    replaceModules(net_cpu, 'nn.ConcatTable', 
    function(orig_mod)
      local cpu_mod = nn.ConcatTable()
      return cpu_mod
    end)
    
    replaceModules(net_cpu, 'nn.Sequential', 
    function(orig_mod)
      local cpu_mod = nn.Sequential()
      return cpu_mod
    end)
    
    replaceModules(net_cpu, 'nn.SpatialZeroPadding', 
    function(orig_mod)
      local cpu_mod = nn.SpatialZeroPadding()
      return cpu_mod
    end)
    
 replaceModules(net_cpu, 'nn.SpatialConvolution', 
    function(orig_mod)
      local cpu_mod = nn.SpatialConvolutionMM(orig_mod.nInputPlane, orig_mod.nOutputPlane,
          orig_mod.kW, orig_mod.kH, orig_mod.dW, orig_mod.dH, orig_mod.padW, orig_mod.padH)
      cpu_mod.weight:copy(orig_mod.weight)
      cpu_mod.bias:copy(orig_mod.bias)
      return cpu_mod
    end)

  replaceModules(net_cpu, 'nn.ReLU', function() return nn.ReLU() end)

  return net_cpu
end


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
function trainModel(steps, batchsize, data, model, criterion, parameters, gradParameters, config, usePhase, noiseLevel, cuda, l1)
    --l2error = torch.Tensor(100)
    --k = 1
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
                local phase = torch.Tensor(1,32,32):zero()
                if usePhase then
                    phase = (torch.rand(1,32,32)*2*math.pi)-math.pi
                end
                local input = {}
                local target = {}
                if noiseLevel then
                    corrupted = activity + randomkit.normal(torch.Tensor(#activity),0,noiseLevel)
                    input = {torch.cmul(corrupted,torch.cos(phase)),torch.cmul(corrupted,torch.sin(phase))}
                    target = {torch.cmul(activity,torch.cos(phase)),torch.cmul(activity,torch.sin(phase))}
                else
                    input = {torch.cmul(activity,torch.cos(phase)),torch.cmul(activity,torch.sin(phase))}
                    target = input
                end
                if cuda then
                   -- l2error = l2error:cuda()
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
                local df_dw = criterion:backward(model.output,target)               
                local a_out = torch.sqrt(torch.pow(output[1],2) + torch.pow(output[2],2))
                display(input[1],model)
                --if k < 100 then
                --    l2error[k] = torch.norm(activity-a_out)
                --    k = k + 1
                --end
                f = f + err +  l1 * torch.norm(x,1)
                model:backward(input,df_dw)
            end
            
            gradParameters:div(batchsize)
            f = f 
            f = f/batchsize
            return f, gradParameters
        end 
        optim.sgd(feval, parameters, config)
        
    end 
    --return l2error
end

iter = 0
function display(input,model)
   iter = iter or 0
   require 'image'
   if iter % 1000 == 0 then
    gnuplot.pngfigure(workdir .. '/ComplexAutoencoder/Data/weights' .. iter .. '.png')
    local weight = model:get(1):parameters()[1]
    local test = image.toDisplayTensor{input = weight, padding = 2, min=weight:min(), max = weight:max()}
    gnuplot.imagesc(test:view(test:size(2), test:size(3)))   
    gnuplot.plotflush()     
   end
   iter = iter + 1
end

function whiteStripes(number)
    white = torch.Tensor(1,32,32):fill(0)
    for i=1,torch.random(1,number) do
        row = torch.random(1,32)
        col = torch.random(1,32)
        white[{1,row,{}}]:fill(255)
        white[{1,{},col}]:fill(255)
    end
    return white    
end    

function createData(amount)
    data = torch.Tensor(amount,1,32,32)
    for i=1,amount do
        data[i] = whiteStripes(5)
    end
    return data
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
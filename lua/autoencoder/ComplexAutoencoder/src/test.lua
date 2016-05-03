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


--local traindata, testdata = loadData(true, 'cifar')
--model = torch.load('/net/store/nbp/projects/phasesim/src_kstandvoss/lua/autoencoder/ComplexAutoencoder/src/model180.net')

local errors = torch.Tensor(10,54)
for i = 1,10 do
  errors[i] = torch.load('../Data/Error' .. i .. '.dat')
end  
x = torch.range(1,54)--torch.linspace(0,54)


ymax = errors:max(1)
ymin = errors:min(1)

test = torch.cat(x,ymin[1],2)
test = torch.cat(test,ymax[1],2)

ymean = errors:mean(1)
gnuplot.pdffigure('../ErrorInterval.pdf')
gnuplot.xlabel('Training Steps')
gnuplot.ylabel('Error')
gnuplot.plot({test, 'filledcurves'}, {x, ymax[1],'~'}, {x, ymin[1],'~'}, {x, ymean[1],'~'})
gnuplot.plotflush()

--[[
autoencoder = nn.Atan2()
criterion = nn.MSECriterion()


function f(x)
  parameters:copy(x)
  -- Do the forward prop
  autoencoder:zeroGradParameters()
  local err = 0
  local output = autoencoder:forward(test)
  err = err + criterion:forward(output, torch.atan2(test[2],test[1]))
  local gradOutput = criterion:backward(output, torch.atan2(test[2],test[1]))
  autoencoder:backward(test, gradOutput)
  return err, grads
end

parameters, grads = autoencoder:getParameters()

test = torch.Tensor(2,10,10):fill(1)
--test1 = torch.Tensor(2,10,10):fill(1+1e-4)
--test2 = torch.Tensor(2,10,10):fill(1-1e-4)
--grad = (autoencoder:forward(test1) + autoencoder:forward(test2)) / (2*1e-4)
 
--out = autoencoder:forward(test)
--crit = criterion:forward(out, torch.atan2(test[2],test[1]))
--gradOut = autoencoder:backward(test, criterion:backward(out,torch.atan2(test[2],test[1])))
--print(nn.JoinTable(1):forward(gradOut))
--print(#grad)
--local err = grad - gradOut[1]
--print(err)
err = nn.Jacobian.testJacobian(autoencoder, test, 0,1, 1e-4)
--local err = optim.checkgrad(f, parameters:clone())
print(err)

--evaluate(model,true)


model = torch.load('/net/store/nbp/projects/phasesim/src_kstandvoss/lua/autoencoder/ComplexAutoencoder/src/model30000.net')



gnuplot.pngfigure(workdir .. '/ComplexAutoencoder/Data/MnistSingle/weights.png')
local weights = model:get(1):getWeights()
local test = image.toDisplayTensor{input = weights, scaleeach=true, padding = 2}
gnuplot.imagesc(test:view(test:size(2), test:size(3)))   
gnuplot.plotflush()

--]]
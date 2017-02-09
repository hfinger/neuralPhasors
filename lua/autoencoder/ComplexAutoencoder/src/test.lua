package.path = package.path .. ";/net/store/nbp/projects/phasesim/src_kstandvoss/lua/autoencoder/ComplexAutoencoder/src/?.lua";
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
require 'cunn'
require 'cutorch'
require 'Atan2'

--[[
--local traindata, testdata = loadData(true, 'cifar')
steps = 10
dataSet = 'MNIST'
traindataPos = torch.Tensor(3,300,400)
traindataNeg = torch.Tensor(3,300,400)
splitdata = torch.Tensor(steps,6,300,400):zero()
print('load')
if dataSet == 'LabelMe' then
  traindataPos = torch.Tensor(3,300,400)
  traindataNeg = torch.Tensor(3,300,400)
  testdata = torch.Tensor(steps,6,300,400):zero()
  print('load')
  --f = 1254
  f = 1010735
  for i=1,steps do
      f = f + 1
      function get()
        print(i)
        --test = mattorch.load('/net/store/nbp/projects/phasesim/workdir/kstandvoss/labelMeWhite/'..i..'/static_newyork_city_urban/IMG_'..f..'.jpg/act1.mat')   
        test = mattorch.load('/net/store/nbp/projects/phasesim/workdir/20130726_Paper/Autoencoder/labelMeWhite/05june05_static_street_boston/p' .. f..'.jpg/act1.mat')
        traindataPos = test['act']:transpose(2,3)
        traindataNeg:copy(traindataPos)
        traindataPos[traindataPos:le(0)] = 0 
        traindataNeg[traindataNeg:ge(0)] = 0
        testdata[i]:sub(1,3):add(100,traindataPos)
        testdata[i]:sub(4,6):add(100,traindataNeg)
      end  
      if pcall(get) then else f = f + 1; i = i -1; end
  end   
  inpD = 6
  w = 400
  h = 300
elseif dataSet == 'MNIST' then
  traindata, testdata = loadData(false, 'mnist')
  testdata = createData(100, testdata, 4, false)
  inpD = 1
  w = 100
  h = 100
elseif dataSet == 'white' then
  traindata = createData(steps+100, 'white', 6)
  testdata = createData(steps+100, 'white', 6)
  inpD = 1
  w = 32
  h = 32  
end
conf = {inputDim = inpD, width=w, height=h}
workdir = '/net/store/nbp/projects/phasesim/src_kstandvoss/lua/autoencoder/ComplexAutoencoder/'
print('Load Model')
model = torch.load('/net/store/nbp/projects/phasesim/src_kstandvoss/lua/autoencoder/ComplexAutoencoder/src/model.net')
print('Model loaded')

evaluate(model,true, false, false, workdir .. 'Data/', false, testdata, 5, 20, dataSet)
--]]
local errors = torch.Tensor(100,125)
for i = 1,100 do
  errors[i] = torch.load('../Data/Errors/Error' .. i .. '.dat')
end  
x = torch.range(1,125)--torch.linspace(0,154)

ymean = errors:mean(1)

std = errors:std()

ymax = ymean + std--errors:max(1)
ymin = ymean - std--errors:min(1)

test = torch.cat(x,ymin[1],2)
test = torch.cat(test,ymax[1],2)

one = errors[34]
two = errors[12]
three = errors[78]
four = errors[1]
five = errors[99]

gnuplot.pdffigure('../ErrorInterval.pdf')
gnuplot.xlabel('Training batch')
gnuplot.ylabel('Error')
gnuplot.plot({'Std', test, 'with filledcurves fill transparent solid 0.5'}, {x, ymax[1],'~'}, {x, ymin[1],'~'}, {'Mean',x, ymean[1],'lt 0 lw 12'}, {'Example1', one, '-'}, {'Example2', two, '-'}, {'Example3', three, '-'}, {'Example4', four, '-'}, {'Example5', five, '-'})
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
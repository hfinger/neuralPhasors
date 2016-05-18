require 'Cos'
require 'Atan2'
require 'Sin'
require 'AddBias'

local Encoder, Parent = torch.class('nn.Encoder', 'nn.Module')

--[[
Torch Neural Network Module for creating the encoding part of an Autoencoder
@n_input input dimension
@n_hidden number of hidden units / dimension of convolution 
@kernel_size dimension of kernel for spatial convolution
]]--
function Encoder:__init(n_input, n_hidden, kernel_size, useLinear)
   Parent.__init(self)
   self.input = n_input 
   self.hidden = n_hidden
   self.kernel_size = kernel_size
   if useLinear then
    self.useLinear = useLinear
   end
   self.encoder = self:build_enc()
   self.train = true
end

--@inp Tensor of size n_input with input to be encoded
function Encoder:updateOutput(inp)
    self.encoder:forward(inp)
    self.output = self.encoder.output 
    return self.output
end

--@input Tensor with input of Backward-call of previous module
--@gradOutput Gradient output of previous module
function Encoder:updateGradInput(inp, gradOutput)
    self.encoder:backward(inp, gradOutput)
    self.gradInput = self.encoder.gradInput
    return self.gradInput
end

function Encoder:cuda()
    self.encoder:cuda()
    self.atan:cuda()
end    

function Encoder:float()
    self.encoder:float()
    self.atan:float()
    self.Bias:float()
end    

function Encoder:updateParameters(learningRate)
    if self.train then       
        self.encoder:updateParameters(learningRate)
        self.convolution1.bias:zero()
        self.convolution1.gradBias:zero()
        self.Bias:updateParameters(learningRate)
    end    
end    

function Encoder:parameters()
    encParams, encGradParams = self.encoder:parameters()
    biasParams, biasGradParams = self.Bias:parameters()
    return {encParams, biasParams}, {encGradParams, biasGradParams}
end 

function Encoder:getWeights()
    return self.convolution1.weight
end

function Encoder:evaluate()
    self.train = false
    self.encoder:evaluate()
    self.Bias:evaluate()
end

function Encoder:training()
    self.train = true
    self.encoder:training()
    self.Bias:training()
end    
    
--Function to build encoder Network
function Encoder:build_enc()
    
    local enc = nn.Sequential()
    local encpar = nn.ConcatTable()

    --imaginary activation part
    local img = nn.Sequential()

    local par1 = nn.ParallelTable()
    local zeroPad = (self.kernel_size-1)/2
    
    if self.useLinear then
        self.convolution1 = nn.Linear(self.input, self.hidden)
        self.convolution2 = nn.Linear(self.input, self.hidden)
    else
        self.convolution1 = nn.SpatialConvolution(self.input, self.hidden, self.kernel_size, self.kernel_size,1,1, zeroPad, zeroPad)
        self.convolution2 = nn.SpatialConvolution(self.input, self.hidden, self.kernel_size, self.kernel_size,1,1, zeroPad, zeroPad)
    end        
    --set bias terms to zero
    self.convolution2:share(self.convolution1,'weight', 'bias', 'gradWeight', 'gradBias')
    self.convolution1.bias:zero() 
    self.convolution1.gradBias:zero() 

    
    local conv1 = nn.Sequential()
   
    conv1:add(self.convolution1)
--    conv1:add(nn.SpatialDropout(0.1))   
    local conv2 = nn.Sequential()
   
    conv2:add(self.convolution2)
--   conv2:add(nn.SpatialDropout(0.1))   
    
    par1:add(conv1) -- w*x
    par1:add(conv2) -- w*y

    img:add(par1)

    local partable = nn.ConcatTable()
    
    local par2 = nn.ParallelTable()
    local seq1 = nn.Sequential()
    

    par2:add(nn.Square()) -- (wx)^2
    par2:add(nn.Square()) -- (wy)^2

    seq1:add(par2)
    seq1:add(nn.CAddTable()) -- (wx)^2 + (wy)^2
    seq1:add(nn.Sqrt()) 
    seq1:add(nn.MulConstant(0.5,true))


    self.atan = nn.Atan2()
    partable:add(self.atan) -- atan(wy,wx)
    partable:add(seq1)

    img:add(partable)


    --real activation part
    local real = nn.Sequential()

    local par2 = nn.ParallelTable()

    par2:add(nn.Square())
    par2:add(nn.Square())


    real:add(par2)
    real:add(nn.CAddTable())  -- x^2 + y^2
    real:add(nn.Sqrt())
    if self.useLinear then
        self.convolution3 = nn.Linear(self.input, self.hidden)
    else
        self.convolution3 = nn.SpatialConvolution(self.input, self.hidden, self.kernel_size, self.kernel_size,1 ,1, zeroPad, zeroPad)
    end    
    self.convolution3:share(self.convolution1,'weight', 'bias', 'gradWeight', 'gradBias')

    local conv3 = nn.Sequential()
    conv3:add(self.convolution3)
--    conv3:add(nn.SpatialDropout(0.1))    
    
    real:add(conv3)
    real:add(nn.MulConstant(0.5,true))

    --
    encpar:add(img)
    encpar:add(real)

    enc:add(encpar)
    enc:add(nn.FlattenTable())
    -------------------------------------------------
    local seperate = nn.ConcatTable()

    --select first column, i.e. phase

    seperate:add(nn.SelectTable(1))

    --select second and third column and add them up
    local activation = nn.Sequential()
    local concat = nn.ConcatTable()
    concat:add(nn.SelectTable(2))                       --Magic in order to seperate and rejoin table 
    concat:add(nn.SelectTable(3))

    activation:add(concat)
    activation:add(nn.CAddTable()) -- add imaginary and real part

    seperate:add(activation)
    --seperate:add(nn.SelectTable(2))
   
    enc:add(seperate)
    -------------------------------------------------

    local relu = nn.ParallelTable()
    relu:add(nn.Identity())
    
    local bias = nn.Sequential()
    self.Bias = nn.AddBias(self.hidden)
    bias:add(self.Bias)
    bias:add(nn.ReLU())
    bias:add(nn.Reshape(self.hidden*100*100))
    bias:add(nn.Normalize(math.huge))
    bias:add(nn.Reshape(self.hidden,100,100))
    --bias:add(nn.Sigmoid())
    --bias:add(nn.SpatialDropout(0.5))
    relu:add(bias) -- use rectified Linear activation

    enc:add(relu)
    --  Transform back to coordinate form
    local transform = nn.ConcatTable()
    x = nn.Sequential()
    y = nn.Sequential()

    local xpar = nn.ParallelTable()
    local ypar = nn.ParallelTable()

    xpar:add(nn.Cos())
    xpar:add(nn.Identity())

    ypar:add(nn.Sin())
    ypar:add(nn.Identity())

    x:add(xpar)
    x:add(nn.CMulTable())

    y:add(ypar)
    y:add(nn.CMulTable())

    transform:add(x)
    transform:add(y)
    enc:add(transform)
    --outputs a table
    --first dimension = x
    --second = y
    return enc
end
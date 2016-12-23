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
function Encoder:__init(n_input, n_hidden, kernel_size, useLinear, input_height, input_width)
   Parent.__init(self)
   self.input = n_input 
   self.hidden = n_hidden
   self.kernel_size = kernel_size
   self.input_height = input_height
   self.input_width = input_width
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
    phase = self.atan.output
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
    encParams, encGradParams = self.convolution1.weight, self.convolution1.gradWeight
    biasParams, biasGradParams = self.Bias:parameters()
    return {encParams, biasParams}, {encGradParams, biasGradParams}
    --return encParams, encGradParams
end 

function Encoder:getWeights()
    return self.convolution1.weight
end

function Encoder:evaluate()
    self.train = false
    self.encoder:evaluate()
    self.Bias:evaluate()
    self.mul1 = nn.MulConstant(0.5,true)
    self.mul2 = nn.MulConstant(0.5,true)
end

function Encoder:training()
    self.train = true
--    self.encoder:training()
--    self.Bias:training()
    self.mul1 = nn.MulConstant(0,true)
    self.mul2 = nn.MulConstant(1,true)
end    

function Encoder:zeroGradParameters()
    self.encoder:zeroGradParameters()
end

    
--Function to build encoder Network
function Encoder:build_enc()
    
    local enc = nn.Sequential()
    local encpar = nn.ConcatTable()

    --imaginary activation part
    local img = nn.Sequential()

    local par1 = nn.ParallelTable()
    local zeroPad = (self.kernel_size-1)/2
    

    self.convolution1 = nn.SpatialConvolution(self.input, self.hidden, self.kernel_size, self.kernel_size,1,1, zeroPad, zeroPad)
    self.convolution1.bias:zero() 
    self.convolution1.gradBias:zero() 
    
    self.convolution2 = self.convolution1:clone('weight', 'bias', 'gradWeight', 'gradBias')


    par1:add(self.convolution1) -- w*x
    par1:add(self.convolution2) -- w*y
 

    img:add(par1)

    local partable = nn.ConcatTable()
    
    local par2 = nn.ParallelTable()
    local seq1 = nn.Sequential()
    

    par2:add(nn.Square()) -- (wx)^2
    par2:add(nn.Square()) -- (wy)^2

    seq1:add(par2)
    seq1:add(nn.CAddTable()) -- (wx)^2 + (wy)^2
    seq1:add(nn.Sqrt()) 
    self.mul1 = nn.MulConstant(0.5,true) 
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
    self.convolution3 = self.convolution1:clone('weight', 'bias', 'gradWeight', 'gradBias')


    real:add(self.convolution3)
    self.mul2 = nn.MulConstant(0.5,true) 
    real:add(self.mul2)

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
    local addImgReal = nn.Sequential()
    local concat = nn.ConcatTable()
    concat:add(nn.SelectTable(2))                       --Magic in order to seperate and rejoin table 
    concat:add(nn.SelectTable(3))

    addImgReal:add(concat)
    addImgReal:add(nn.CAddTable()) -- add imaginary and real part

    seperate:add(addImgReal)
   
    enc:add(seperate)
    -------------------------------------------------

    local relu = nn.ParallelTable()
    relu:add(nn.Identity())
    
    local activation = nn.Sequential()
    self.Bias = nn.AddBias(self.hidden) --nn.Add(self.hidden, self.input_width,self.input_height) --
    --activation:add(nn.L1Penalty(0.0001))
    activation:add(self.Bias)
    activation:add(nn.ReLU())
    --activation:add(nn.SpatialDropout(0.1))
    relu:add(activation) -- use rectified Linear activation

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
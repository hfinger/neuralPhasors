local Decoder, Parent = torch.class('nn.Decoder', 'nn.Module')

function Decoder:__init(n_input, n_hidden, kernel_size,input_width,input_height)
   Parent.__init(self)
   self.input = n_input 
   self.hidden = n_hidden
   self.kernel_size = kernel_size
   self.input_width = input_width
   self.input_height = input_height
   self.decoder = self:build_dec()
end

--@inp Tensor of size n_input with input to be encoded
function Decoder:updateOutput(inp)
    self.decoder:forward(inp)
    self.output = self.decoder.output 
    return self.output
end

--@input Tensor with input of Backward-call of previous module
--@gradOutput Gradient output of previous module
function Decoder:updateGradInput(inp, gradOutput)
    self.decoder:backward(inp, gradOutput)
    self.gradInput = self.decoder.gradInput
    return self.gradInput
end

function Decoder:updateParameters(learningRate)
    self.decoder:updateParameters(learningRate)
end  

function Decoder:parameters()
    encParams, encGradParams =  self.convolution1.weight, self.convolution1.gradWeight
    biasParams, biasGradParams = self.Bias:parameters()
    return {encParams, biasParams}, {encGradParams, biasGradParams}
end 

function Decoder:cuda()
    self.decoder:cuda()
end  

function Decoder:float()
    self.decoder:float()
end  

--Function to build decoder Network
function Decoder:build_dec()

    --initialize Decoder\
    local enc = nn.Sequential()
    local encpar = nn.ConcatTable()

    --imaginary activation part
    local img = nn.Sequential()

    local par1 = nn.ParallelTable()
    local zeroPad = (self.kernel_size-1)/2
    

    self.convolution1 = nn.SpatialConvolution(self.input, self.hidden, self.kernel_size, self.kernel_size,1,1, zeroPad, zeroPad)
    self.convolution2 = nn.SpatialConvolution(self.input, self.hidden, self.kernel_size, self.kernel_size,1,1, zeroPad, zeroPad)
       
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
    self.convolution3 = nn.SpatialConvolution(self.input, self.hidden, self.kernel_size, self.kernel_size,1 ,1, zeroPad, zeroPad)
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
    self.Bias = nn.AddBias(self.hidden) --vnn.Add(self.hidden,self.input_width,self.input_height)--
    activation:add(self.Bias)
    activation:add(nn.ReLU())
    relu:add(activation)
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
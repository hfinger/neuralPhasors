local Decoder, Parent = torch.class('nn.Decoder', 'nn.Module')

function Decoder:__init(n_input, n_hidden, kernel_size)
   Parent.__init(self)
   self.input = n_input 
   self.hidden = n_hidden
   self.kernel_size = kernel_size
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
    convolution1.bias:zero()
    convolution1.gradBias:zero() 
end  

function Decoder:parameters()
    return self.decoder:parameters()
end 

function Decoder:cuda()
    self.decoder:cuda()
end  

--Function to build decoder Network
function Decoder:build_dec()

    --initialize Decoder
    local dec = nn.ParallelTable()
    convolution1 = nn.SpatialConvolution(self.input, self.hidden, self.kernel_size, self.kernel_size)
    convolution2 = nn.SpatialConvolution(self.input,self.hidden, self.kernel_size, self.kernel_size)
    
    local zeroPad = (self.kernel_size-1)/2
    local decConv1 = nn.Sequential()
    decConv1:add(nn.SpatialZeroPadding(zeroPad,zeroPad,zeroPad,zeroPad))
    decConv1:add(convolution1)
    
    local decConv2 = nn.Sequential()
    decConv2:add(nn.SpatialZeroPadding(zeroPad,zeroPad,zeroPad,zeroPad))
    decConv2:add(convolution2)
    
    convolution2:share(convolution1,'weight', 'gradWeight', 'bias', 'gradBias')
    convolution1.bias:zero()
    convolution1.gradBias:zero() 

    
    dec:add(decConv1)
    dec:add(decConv2)
    
    return dec
end
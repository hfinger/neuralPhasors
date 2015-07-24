local AddBias, Parent = torch.class('nn.AddBias', 'nn.Module')

--[[
Torch Neural Network Module to add learnable Bias-term
@size dimension of input Tensor
]]--
function AddBias:__init(size)
   Parent.__init(self)
   self.bias = torch.Tensor(size,1,1):zero()
   self.gradBias = torch.Tensor(size,1,1):zero() 
end

--@input Tensor with first dimension of size @size 
function AddBias:updateOutput(input)
    exB = torch.expand(self.bias,(#self.bias)[1],(#input)[2],(#input)[3])
    self.output = torch.add(exB,input)
    return self.output
end
 
--@input Tensor with input of Backward-call of previous module
--@gradOutput Gradient output of previous module
function AddBias:updateGradInput(input, gradOutput)
    self.gradInput = gradOutput:clone() --gradient stays the same
    return self.gradInput
end

--@input Tensor with input of Backward-call of previous module
--@gradOutput Gradient output of previous module 
--@scale Factor for parameter accumulation
function AddBias:accGradParameters(input, gradOutput, scale)
   scale = scale or 1
   local gradOutput = gradOutput:view((#self.bias)[1], -1) --take gradient output with dimensions like bias
   self.gradBias:add(scale,gradOutput:sum(2)) --accumulate gradient output into gradient-bias
end
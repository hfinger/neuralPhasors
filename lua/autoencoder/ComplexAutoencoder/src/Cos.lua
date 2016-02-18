local Cos, Parent = torch.class('nn.Cos', 'nn.Module')

--[[
Torch Neural Network Module to apply cosine on input Tensor
]]--
function Cos:__init()
   Parent.__init(self)
end

--Apply cosine
--@input Tensor 
function Cos:updateOutput(input)
    self.output = torch.cos(input) 
    return self.output
end

--Gradient
--@input Tensor with input of previous module
--@gradOutput Tensor with gradient of previous module
function Cos:updateGradInput(input, gradOutput)
    self.gradInput = torch.cmul(gradOutput,-torch.sin(input))  -- grad(cos) = -sin
    return self.gradInput
end
 
 
function Cos:reset()
end
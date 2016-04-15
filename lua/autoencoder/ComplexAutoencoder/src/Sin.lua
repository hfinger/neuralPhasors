local Sin, Parent = torch.class('nn.Sin', 'nn.Module')

--[[
Torch Neural Network Module to apply sine on input Tensor
]]--
function Sin:__init()
   Parent.__init(self)
end

--Apply sine
--@input Tensor 
function Sin:updateOutput(input)
    self.output = torch.sin(input) 
    return self.output
end

--Gradient
--@input Tensor with input of previous module
--@gradOutput Tensor with gradient of previous module
function Sin:updateGradInput(input, gradOutput)
    self.gradInput = torch.cmul(gradOutput,torch.cos(input)) -- grad(sin) = cos
    return self.gradInput
end
 
 
function Sin:reset()
end
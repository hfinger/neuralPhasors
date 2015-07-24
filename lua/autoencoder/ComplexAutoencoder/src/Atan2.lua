local Atan2, Parent = torch.class('nn.Atan2', 'nn.Module')

--[[
Torch Neural Network Module to apply atan2 on input Tensor
]]--
function Atan2:__init()
   Parent.__init(self)
end

--Apply Atan2
--@input with second dimension = 2 for x and y
function Atan2:updateOutput(input)
    self.output = torch.atan2(input[2],input[1]) 
    return self.output
end
 
--Gradient
--@input Tensor with input of previous module
--@gradOutput Tensor with gradient of previous module
function Atan2:updateGradInput(input, gradOutput)
    
    -- xgrad(Atan2) = -y/(x^2+y^2)
    addx = torch.add(torch.pow(input[1],2),torch.pow(input[2],2)) 
    divx = -torch.cdiv(input[2],addx)  
    x = torch.cmul(gradOutput,divx)
    
    -- ygrad(Atan2) = x/(x^2+y^2)
    addy = torch.add(torch.pow(input[1],2),torch.pow(input[2],2))
    divy = torch.cdiv(input[1],addy)
    y = torch.cmul(gradOutput,divy)
    self.gradInput = {x,y}
    return self.gradInput
end
 
 
function Atan2:reset()
end
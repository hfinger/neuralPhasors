require 'nn'

local Atan2, Parent = torch.class('nn.Atan2', 'nn.Module')
 
function Atan2:__init()
   Parent.__init(self)
end

function Atan2:updateOutput(input)
    self.output = torch.atan2(input[2],input[1]) 
    return self.output
end
 
function Atan2:updateGradInput(input, gradOutput)
    addx = torch.add(torch.pow(input[1],2),torch.pow(input[2],2))
    divx = -torch.cdiv(input[2],addx)  
    x = torch.cmul(gradOutput,divx)

    addy = torch.add(torch.pow(input[1],2),torch.pow(input[2],2))
    divy = torch.cdiv(input[1],addy)
    y = torch.cmul(gradOutput,divy)
    self.gradInput = {x,y}
    return self.gradInput
end
 
 
function Atan2:reset()
end
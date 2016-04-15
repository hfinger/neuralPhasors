local Atan2, Parent = torch.class('nn.Atan2', 'nn.Module')

--[[
Torch Neural Network Module to apply atan2 on input Tensor
]]--
function Atan2:__init()
   Parent.__init(self)
   self.useCuda = false
end

--Apply Atan2
--@input with second dimension = 2 for x and y
function Atan2:updateOutput(input)

    self.output = torch.atan2(input[2],input[1]) 
    if self.output:ne(self.output):sum() > 1 then
        self.output = (torch.rand(#self.output)*2*math.pi)-math.pi
    end
    if self.useCuda then
        self.output = self.output:cuda()
    end    
    return self.output
end

function Atan2:cuda()
    self.useCuda = true
end 

function Atan2:float()
    self.useCuda = false
    self.output = self.output:float()
end 

--Gradient
--@input Tensor with input of previous module
--@gradOutput Tensor with gradient of previous module
function Atan2:updateGradInput(input, gradOutput)
    if torch.all(input[1]:eq(0)) and torch.all(input[2]:eq(0)) then
        x = torch.Tensor(#input[1]):zero()
        y = torch.Tensor(#input[2]):zero()
        if self.useCuda then
            x = x:cuda()
            y = y:cuda()
        end
    else
        -- xgrad(Atan2) = -y/(x^2+y^2)
        addx = torch.add(torch.pow(input[1],2),torch.pow(input[2],2)) 
        divx = -torch.cdiv(input[2],addx)  
        x = torch.cmul(gradOutput,divx)
        if x:ne(x):sum() > 1 then
            x = x:zero()
        end
        -- ygrad(Atan2) = x/(x^2+y^2)
        addy = torch.add(torch.pow(input[1],2),torch.pow(input[2],2))
        divy = torch.cdiv(input[1],addy)
        y = torch.cmul(gradOutput,divy) 
        if y:ne(y):sum() > 1 then
            y =y:zero()
        end      
    end
    self.gradInput = {x,y}
    return self.gradInput
end
 
 
function Atan2:reset()
end
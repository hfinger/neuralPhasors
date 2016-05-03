local Atan2, Parent = torch.class('nn.Atan2', 'nn.Module')

--[[
Torch Neural Network Module to apply atan2 on input Tensor
]]--
function Atan2:__init()
   Parent.__init(self)
   self.useCuda = true
end

--Apply Atan2
--@input with second dimension = 2 for x and y
function Atan2:updateOutput(input)
    wx = input[1]:clone()
    wy = input[2]:clone()
  
    zeros = torch.cmul(wx:eq(0),wy:eq(0)):byte()
    self.output = torch.atan2(wy,wx)
    random = ((torch.rand(zeros:sum())*2*math.pi)-math.pi):cuda()
    self.output[zeros] = random     
    if self.output:ne(self.output):sum() > 1 then
            print('Forward Phase undefined')
            undef = self.output[self.output:ne(self.output)]
            random = ((torch.rand(#undef)*2*math.pi)-math.pi):cuda()
            self.output[self.output:ne(self.output)] = random
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
    
    wx = input[1]:clone()
    wy = input[2]:clone()

    if torch.all(wx:eq(0)) and torch.all(wy:eq(0)) then
        x = torch.Tensor(#input[1]):zero()
        y = torch.Tensor(#input[2]):zero()
        if self.useCuda then
            x = x:cuda()
            y = y:cuda()
        end
        print('here')
    else
        -- xgrad(Atan2) = -y/(x^2+y^2)
        addx = torch.add(torch.pow(wx,2),torch.pow(wy,2)) 
        divx = -torch.cdiv(wy,addx)
        x = torch.cmul(gradOutput,divx)
        if x:ne(x):sum() > 1 then
            print('Phase undefined')            
            x = x:zero()
        end
        -- ygrad(Atan2) = x/(x^2+y^2)
        addy = torch.add(torch.pow(wx,2),torch.pow(wy,2))
        divy = torch.cdiv(wx,addy)
        y = torch.cmul(gradOutput,divy)
        if y:ne(y):sum() > 1 then
            print('Phase undefined')
            y =y:zero()
        end      
    end
    self.gradInput = {x,y}
    return self.gradInput
end
 
 
function Atan2:reset()
end
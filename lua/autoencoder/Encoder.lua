require 'nn'
require 'image'
require 'unsup'
require 'optim'

local Encoder, Parent = torch.class('nn.Encoder', 'nn.Module')
 
function Encoder:__init()
   Parent.__init(self)
    
    
end
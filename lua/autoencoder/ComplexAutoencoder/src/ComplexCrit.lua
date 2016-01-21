local ComlexCrit, Parent = torch.class('nn.ComplexCrit', 'nn.Criterion')
 
function ComplexCrit:__init()
   Parent.__init(self)
end

function ComplexCrit:updateOutput(input,target)
    
    a_in = torch.sqrt(torch.add(torch.pow(input[1],2),torch.pow(input[2],2)))
    a_tar = torch.sqrt(torch.add(torch.pow(target[1],2),torch.pow(target[2],2)))
    
    phase_in = torch.acos(torch.div(input[1],a_in))
    phase_tar = torch.acos(torch.div(target[1],a_tar))
    
    a_err = torch.pow((a_tar - a_in),2) 
    a_err = xerr:sum()
    phase_err = (-torch.cos(phase_tar-phase_in)+1)*a_in*a_tar
    phase_err = phase_err:sum()
    
    self.output = (1/(#input[1])) * (a_rr+phase_err)
    return self.output
end
 
function ComplexCrit:updateGradInput(input, gradOutput)
    self.gradInput = torch.cmul(gradOutput,-torch.sin(input)) 
    return self.gradInput
end
 
 
function Cos:reset()
end
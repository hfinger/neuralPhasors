require 'unsup';
require 'nn';
require 'gnuplot';
require 'Encoder';
require 'Decoder';
require 'nngraph';
require 'image';
require 'optim';
require 'functions'
require 'image'

traindata, testdata = loadData{normalizeMean = false} 
autoencoder = torch.load('model500.net')


activity = testdata[101]
activities = torch.Tensor(12,32,32)
activities[1] = activity
phases = torch.Tensor(10,3,32,32)
reconstructionErrors = torch.Tensor(10)
circularVariances = torch.Tensor(10)
print('start iterating')
for i = 0,9 do
    
    local phase = (torch.rand(1,32,32)*2*math.pi)-math.pi
    local inp = {torch.cmul(activity,torch.cos(phase)),torch.cmul(activity,torch.sin(phase))}

    local out = autoencoder:forward(inp)
    local phase_out = torch.atan2(out[2],out[1])
    
    for j = 1,(10*i) do
        out = autoencoder:forward({torch.cmul(activity,torch.cos(phase_out)),torch.cmul(activity,torch.sin(phase_out))})
        phase_out = torch.atan2(out[2],out[1])
    end
    
    local hsv = torch.Tensor(3,32,32)
    hsv[1]:copy((phase_out+math.pi)/(2*math.pi))
    hsv[2] = torch.Tensor(32,32):fill(0.5)
    hsv[3] = torch.Tensor(32,32):fill(0.5)
    local rgb = image.hsv2rgb(hsv)

    phases[i+1] = rgb
    
    local a_out = torch.sqrt( torch.pow(out[1],2) + torch.pow(out[2],2))
    --print(a_out)
    activities[i+2] = a_out
    reconstructionErrors[i+1] = torch.norm(activity-a_out)
    
    local phi = phase_out[torch.gt(activity,0)] -- find pixels of stimulus: activity>0
    local var = 1 - (1/phi:size(1))* torch.sqrt((torch.sum(torch.cos(phi))^2 + torch.sum(torch.sin(phi))^2))
    circularVariances[i+1] = var
end

local phase = torch.Tensor(#activity):zero()
local inp = {torch.cmul(activity,torch.cos(phase)),torch.cmul(activity,torch.sin(phase))}
local out = autoencoder:forward(inp)
local a_out = torch.sqrt( torch.pow(out[1],2) + torch.pow(out[2],2 ))
print('Error Zero Phase: ' .. torch.norm(activity-a_out))
activities[12] = a_out

print('create Plot')
--
gnuplot.pdffigure('PhaseIterations.pdf')
gnuplot.xlabel('Phase iterations * 10e^1')
gnuplot.ylabel('L2 Norm')
x = torch.linspace(0,9,10)
gnuplot.plot(x,reconstructionErrors)
gnuplot.plotflush()
--
gnuplot.pdffigure('circularVariances.pdf')
gnuplot.xlabel('Phase iterations * 20e^1')
gnuplot.ylabel('circular Variance')
gnuplot.plot(x,circularVariances)
gnuplot.plotflush()
--
image.save('Reconstructions.png', image.toDisplayTensor(activities))
image.save('Phases.png', image.toDisplayTensor(phases))

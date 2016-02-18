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
require 'mattorch';
loaded = false
function cuda()
    require 'cutorch';
    require 'cunn';
end
if pcall(cuda) then loaded = true else print('Cuda cannot be required') end

   
function evaluate(model, useCuda, trainerror, steps, workdir, fileId)

    traindata, testdata = loadData{normalizeMean = false} 
    --traindata = createData(1000):div(255)
    --testdata = createData(1000):div(255)
    if model then
        autoencoder = model
    else
        autoencoder = torch.load('model.net')
    end
    
    if not workdir then
        workdir = '/net/store/nbp/projects/phasesim/src_kstandvoss/lua/autoencoder/ComplexAutoencoder/Data/'
    end
    
    if not fileId then
        fileId = ''
    end
   
    activity = whiteStripes(5)
    --activity = torch.cmax(testdata[143],image.translate(testdata[72],4,2))
    activities = torch.Tensor(12,32,32)
    activities[1] = activity
    phases = torch.Tensor(10,3,32,32)
    reconstructionErrors = torch.Tensor(10)
    circularVariances = torch.Tensor(10)
    weights = torch.Tensor(10)
    if loaded and useCuda then
        activity = activity:cuda()
        activities = activities:cuda()
        phases = phases:cuda()
        reconstructionErrors = reconstructionErrors:cuda()
        circularVariances = circularVariances:cuda()
        weights = weights:cuda()
    end    
    print('start iterating')
    for i = 0,9 do
        
        local phase = (torch.rand(1,32,32)*2*math.pi)-math.pi
        if loaded and useCuda then
            phase = phase:cuda()
        end
        local inp = {torch.cmul(activity,torch.cos(phase)),torch.cmul(activity,torch.sin(phase))}

        local out = autoencoder:forward(inp)
        local phase_out = torch.atan2(out[2],out[1])
        

        
        for j = 1,(10*i) do
            out = autoencoder:forward({torch.cmul(activity,torch.cos(phase_out)),torch.cmul(activity,torch.sin(phase_out))})
            phase_out = torch.atan2(out[2],out[1])
        end
        
        local a_out = torch.sqrt( torch.pow(out[1],2) + torch.pow(out[2],2))   
          
        local hsv = torch.Tensor(3,32,32)
        hsv[1]:copy((phase_out[1]+math.pi)/(2*math.pi))
        hsv[2] = torch.Tensor(32,32):fill(0.5)
        hsv[3] = torch.Tensor(32,32):fill(0.5)
        local rgb = image.hsv2rgb(hsv)
    
        phases[i+1] = torch.cmul(rgb, a_out:double():expand(3,32,32))
        gnuplot.pngfigure(workdir .. 'weights'..fileId..'.png')
        local weights = autoencoder:get(1):parameters()[1]
        gnuplot.imagesc(weights[1]:reshape(7,7),'gray')
        gnuplot.imagesc(weights[2]:reshape(7,7),'gray')
        gnuplot.plotflush()
        --print(a_out)
        activities[i+2] = a_out
        reconstructionErrors[i+1] = torch.norm(activity-a_out)
        
        local phi = phase_out[torch.gt(activity,0)] -- find pixels of stimulus: activity>0
        local var = 1 - (1/phi:size(1))* torch.sqrt((torch.sum(torch.cos(phi))^2 + torch.sum(torch.sin(phi))^2))
        circularVariances[i+1] = var
    end
    local phase = torch.Tensor(#activity):zero()
    if loaded and useCuda then
            phase = phase:cuda()
    end
    local inp = {torch.cmul(activity,torch.cos(phase)),torch.cmul(activity,torch.sin(phase))}
    local out = autoencoder:forward(inp)
    local a_out = torch.sqrt(torch.pow(out[1],2) + torch.pow(out[2],2))
    print('Error Zero Phase: ' .. torch.norm(activity-a_out))
    activities[12] = a_out
    
    print('create Plot')
    --
    reconstructionErrors = reconstructionErrors:float()
    --mattorch.save('reconstructionErrors.mat',{reconstructionErrors})
    gnuplot.pdffigure(workdir .. 'PhaseIterations'..fileId..'.pdf')
    gnuplot.xlabel('Phase iterations * 10e^1')
    gnuplot.ylabel('L2 Norm')
    x = torch.linspace(0,9,10)
    gnuplot.plot(x,reconstructionErrors)
    gnuplot.plotflush()

    
    --
    gnuplot.pdffigure(workdir .. 'circularVariances'..fileId..'.pdf')
    gnuplot.xlabel('Phase iterations * 10e^1')
    gnuplot.ylabel('circular Variance')
    gnuplot.plot(x,circularVariances:float())
    gnuplot.plotflush()
    --
    if trainerror then
        gnuplot.pdffigure(workdir .. 'Training_Error'..fileId..'.pdf')
        gnuplot.xlabel('Training Steps')
        gnuplot.ylabel('L2 Norm')
        x = torch.linspace(0,steps,steps/100)
        gnuplot.plot(x,trainerror:float())
        gnuplot.plotflush()
    end
    --    
    image.save(workdir .. 'Reconstructions'..fileId..'.png', image.toDisplayTensor(activities:float()))
    image.save(workdir .. 'Phases'..fileId..'.png', image.toDisplayTensor(phases:float()))
end
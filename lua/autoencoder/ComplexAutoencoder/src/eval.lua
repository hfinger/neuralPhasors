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

   
function evaluate(model, useCuda, trainerror, steps, workdir, fileId, testdata)


    if model then
        autoencoder = model
    else
        autoencoder = torch.load('/net/store/nbp/projects/phasesim/src_kstandvoss/lua/autoencoder/ComplexAutoencoder/src/model50000.net')
    end
    
    if not workdir then
        workdir = '/net/store/nbp/projects/phasesim/src_kstandvoss/lua/autoencoder/ComplexAutoencoder/Data/'
    end
    
    if not fileId then
        fileId = ''
    end
   
    activity = testdata[23]
    activities = torch.Tensor(12,activity:size()[2],activity:size()[3])
    activitiesNoPhase = torch.Tensor(12,activity:size()[2],activity:size()[3])
    activities[1] = activity 
    phases = torch.Tensor(10,3,activity:size()[2],activity:size()[3])
    reconstructionErrors = torch.Tensor(10)
    reconstructionErrorsNoPhase = torch.Tensor(10)
    circularVariances = torch.Tensor(10)

    
    if loaded and useCuda then
        activity = activity:cuda()
        activities = activities:cuda()
        activitiesNoPhase = activitiesNoPhase:cuda()
        phases = phases:cuda()
        reconstructionErrors = reconstructionErrors:cuda()
        reconstructionErrorsNoPhase = reconstructionErrorsNoPhase:cuda()
        circularVariances = circularVariances:cuda()
    end    
    
    print('start iterating')
    for i = 0,9 do
        
        -- Phase iterations 
        local phase = (torch.rand(1,activity:size()[2],activity:size()[3])*2*math.pi)-math.pi

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
        
        gnuplot.pdffigure('hist' .. i .. '.pdf')
        gnuplot.hist(phase_out[torch.gt(activity,0)])
        gnuplot.plotflush()       
        
        local a_out = torch.sqrt( torch.pow(out[1],2) + torch.pow(out[2],2))     
          
        local hsv = torch.Tensor(3,activity:size()[2],activity:size()[3])
        hsv[1]:copy((phase_out[1]+math.pi)/(2*math.pi))
        hsv[2] = torch.Tensor(activity:size()[2],activity:size()[3]):fill(0.5)
        hsv[3] = torch.Tensor(activity:size()[2],activity:size()[3]):fill(0.5)
        local rgb = image.hsv2rgb(hsv)

        rgb = torch.cmul(rgb, a_out:double():expand(3,activity:size()[2],activity:size()[3])) 
        phases[i+1] = rgb 
        activities[i+2] = a_out
        reconstructionErrors[i+1] = torch.norm(activity-a_out)
        
        local phi = phase_out[torch.gt(activity,0)] -- find pixels of stimulus: activity>0
        
        
        local var = 1 - (1/phi:size(1))* torch.sqrt((torch.sum(torch.cos(phi))^2 + torch.sum(torch.sin(phi))^2))
        circularVariances[i+1] = var
        
        --Iterations zero phase
        zeroPhase = torch.Tensor(#activity):zero()
        if loaded and useCuda then
                zeroPhase = zeroPhase:cuda()
        end 
        
        inp = {torch.cmul(activity,torch.cos(zeroPhase)),torch.cmul(activity,torch.sin(zeroPhase))}
        out = autoencoder:forward(inp)

        
        for j = 1,(10*i) do
            out = autoencoder:forward({torch.cmul(activity,torch.cos(zeroPhase)),torch.cmul(activity,torch.sin(zeroPhase))})   
        end
        a_out = torch.sqrt( torch.pow(out[1],2) + torch.pow(out[2],2))   
        reconstructionErrorsNoPhase[i+1] = torch.norm(activity-a_out)
        activitiesNoPhase[i+1] = a_out
    end
    local inp = {torch.cmul(activity,torch.cos(zeroPhase)),torch.cmul(activity,torch.sin(zeroPhase))}
    local out = autoencoder:forward(inp)
    local a_out = torch.sqrt(torch.pow(out[1],2) + torch.pow(out[2],2))
    print('Error Zero Phase: ' .. torch.norm(activity-a_out))
    activities[12] = a_out
    
    print('create Plot')
    --
    reconstructionErrors = reconstructionErrors:double()
    reconstructionErrorsNoPhase = reconstructionErrorsNoPhase:double()
    --mattorch.save('reconstructionErrors.mat',{reconstructionErrors})
    gnuplot.pdffigure(workdir .. 'PhaseIterations'..fileId..'.pdf')
    gnuplot.xlabel('Phase iterations * 10e^1')
    gnuplot.ylabel('L2 Norm')
    x = torch.linspace(0,9,10)
    gnuplot.plot({'Phase', x,reconstructionErrors},{'No phase', x, reconstructionErrorsNoPhase})
    gnuplot.plotflush()
    
    

    
    --
    gnuplot.pdffigure(workdir .. 'circularVariances'..fileId..'.pdf')
    gnuplot.xlabel('Phase iterations * 10e^1')
    gnuplot.ylabel('circular Variance')
    gnuplot.plot(x,circularVariances:double())
    gnuplot.plotflush()
    --
    if trainerror then
        trainerror = trainerror:double()
        gnuplot.pdffigure(workdir .. 'ErrorTrain'..fileId..'.pdf')
        gnuplot.xlabel('Training Steps')
        gnuplot.ylabel('L2 Norm')
        x = torch.linspace(0,steps,steps/10 + 3)
        gnuplot.plot(x,trainerror)
        gnuplot.plotflush()
    end
    --    
    image.save(workdir .. 'Reconstructions'..fileId..'.png', image.toDisplayTensor(activities:double()))
    image.save(workdir .. 'NoPhase'..fileId..'.png', image.toDisplayTensor(activitiesNoPhase:double()))
    image.save(workdir .. 'Phases'..fileId..'.png', image.toDisplayTensor(phases:double()))
end
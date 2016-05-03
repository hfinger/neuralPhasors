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

   
function evaluate(model, useCuda, trainerror, steps, workdir, fileId, testdata, params, epochs, batchSize, set)    

    if not epochs then epochs = 1 else epochs = epochs end
    
    if model then
        autoencoder = model
    end
    
    if not workdir then
        workdir = '/net/store/nbp/projects/phasesim/src_kstandvoss/lua/autoencoder/ComplexAutoencoder/Data/'
    end
    
    if not fileId then
        fileId = ''
    end
    
    if not testdata then
        testdata = torch.load('/net/store/nbp/projects/phasesim/src_kstandvoss/lua/autoencoder/ComplexAutoencoder/src/labelMe4.dat')
    end
    
    
    --testdata = createData(100,testdata,5)
    --activity = testdata[23]
    --test = mattorch.load('/net/store/nbp/projects/phasesim/workdir/20130726_Paper/Autoencoder/labelMeWhite/05june05_static_street_boston/p1010736.jpg/act1.mat')
    
    activity = testdata[1] 
    local indim
    if set == 'LabelMe' then
      indim = 3
      act = function(activity, mode) if mode == 'plus' then return activity:sub(1,3)+ activity:sub(4,6) else return activity:sub(1,3) - activity:sub(4,6) end end
    elseif set == 'MNIST' then
      indim = 1
      act = function(activity) return activity end
    elseif set == 'white' then
      indim = 1
      act = function(activity) return activity end 
    end     
    activities = torch.Tensor(12, indim, activity:size()[2], activity:size()[3])
    activities[1] = act(activity, 'plus') 
    phases = torch.Tensor(10,3, activity:size()[2],activity:size()[3])
    reconstructionErrors = torch.Tensor(10)
    reconstructionErrorsNoPhase = torch.Tensor(10)
    circularVariances = torch.Tensor(10)
    synchronization = torch.Tensor(10)
    
    if loaded and useCuda then
        activity = activity:cuda()
        activities = activities:cuda()
        phases = phases:cuda()
        reconstructionErrors = reconstructionErrors:cuda()
        reconstructionErrorsNoPhase = reconstructionErrorsNoPhase:cuda()
        circularVariances = circularVariances:cuda()
    end    
    
    print('start iterating')
    local phase = (torch.rand(#activity)*2*math.pi)-math.pi
    if loaded and useCuda then
      phase = phase:cuda()
    end
    local atan = nn.Atan2()
    local inp = {torch.cmul(activity,torch.cos(phase)), torch.cmul(activity,torch.sin(phase))}
    local out = autoencoder:forward(inp)
    local a_out = torch.sqrt(torch.pow(out[1],2) + torch.pow(out[2],2))    
    local phase_out = torch.atan2(out[2],out[1])
    
    j = 1
    for i = 0,90 do
        -- Phase iterations 
        
        if i%10 == 0 then 
          --[[
          gnuplot.pdffigure('../Data/histPhase' .. j .. '.pdf')
          gnuplot.hist(phase_out)
          gnuplot.plotflush() 
          gnuplot.pdffigure('../Data/histHidPhase' .. j .. '.pdf')
          gnuplot.hist(autoencoder:get(1).atan.output)
          gnuplot.plotflush()]]--
          local sin
          local cos
          local circMean
          if set == 'LabelMe' then
            sin = torch.cmul(a_out[1],torch.sin(phase_out[1]))
            cos = torch.cmul(a_out[1],torch.cos(phase_out[1]))
            for i=1,a_out:size(1) do 
              sin = sin + torch.cmul(a_out[i],torch.sin(phase_out[i])) 
              cos = cos + torch.cmul(a_out[i],torch.cos(phase_out[i])) 
            end  
            circMean = torch.atan2(sin,cos)
          else
            circMean = phase_out
          end    
          local hsv = torch.Tensor(3,activity:size(2),activity:size(3))
          hsv[1]:copy((circMean+math.pi)/(2*math.pi))
          hsv[2] = torch.Tensor(activity:size(2),activity:size(3)):fill(0.5)
          hsv[3] = torch.Tensor(activity:size(2),activity:size(3)):fill(0.5)
          local rgb = image.hsv2rgb(hsv)
          if set ~= 'LabelMe' then
            rgb = torch.cmul(rgb, a_out:expand(3, a_out:size(2), a_out:size(3)):double())
          end  
          a = act(a_out, 'minus')
          

          phases[j] = rgb 
          activities[j+1] = a
          reconstructionErrors[j] = torch.norm(activity-a_out)
          --synchronization[j] = torch.abs(torch.sum(autoencoder:get(1).output[1]+autoencoder:get(1).output[2])/ torch.sum(torch.sqrt(torch.pow(autoencoder:get(1).output[1],2) + torch.pow(autoencoder:get(1).output[2],2))))       
          j = j+1
        end
        --local phi = phase_out[torch.gt(activity,0)] -- find pixels of stimulus: activity>0  
        --local var = 1 - (1/phi:size(1))* torch.sqrt((torch.sum(torch.cos(phi))^2 + torch.sum(torch.sin(phi))^2))
        --circularVariances[i+1] = var
        
       out = autoencoder:forward({torch.cmul(activity,torch.cos(phase_out)),torch.cmul(activity,torch.sin(phase_out))})
       a_out = torch.sqrt(torch.pow(out[1],2) + torch.pow(out[2],2))  
       phase_out = torch.atan2(out[2],out[1])
    end
    --torch.save(workdir..'synchronization.dat', synchronization)
    --print(torch.sum(synchronization))
    --[[
    zeroPhase = torch.Tensor(#activity):zero()
    if loaded and useCuda then
            zeroPhase = zeroPhase:cuda()
    end 
    local inp = {torch.cmul(activity,torch.cos(zeroPhase)),torch.cmul(activity,torch.sin(zeroPhase))}
    local out = autoencoder:forward(inp)
    local a_out = torch.sqrt(torch.pow(out[1],2) + torch.pow(out[2],2))
    print('Error Zero Phase: ' .. torch.norm(activity-a_out))
    activities[12] = a_out:add(0.5)
    --]]
    print('create Plot')
    --[[
    reconstructionErrors = reconstructionErrors:double()
    reconstructionErrorsNoPhase = reconstructionErrorsNoPhase:double()
    --mattorch.save('reconstructionErrors.mat',{reconstructionErrors})
    gnuplot.pdffigure(workdir .. 'PhaseIterations'..fileId..'.pdf')
    gnuplot.xlabel('Phase iterations * 10e^1')
    gnuplot.ylabel('L2 Norm')
    x = torch.linspace(0,9,10)
    gnuplot.plot({'Phase', x,reconstructionErrors},{'No phase', x, reconstructionErrorsNoPhase})
    gnuplot.plotflush()
    --]]
    

    
    --[[
    gnuplot.pdffigure(workdir .. 'circularVariances'..fileId..'.pdf')
    gnuplot.xlabel('Phase iterations * 10e^1')
    gnuplot.ylabel('circular Variance')
    gnuplot.plot(x,circularVariances:double())
    gnuplot.plotflush()]]--
    torch.save(workdir .. 'Error' ..fileId .. '.dat', trainerror)
    image.save(workdir .. 'Reconstructions'..fileId..'.png', image.toDisplayTensor{input=activities:double(), scaleeach = true, padding=2})
    image.save(workdir .. 'Phases'..fileId..'.png', image.toDisplayTensor{input=phases:double(),scaleeach=true, padding=2})
    
    if trainerror then
        trainerror = trainerror:double()
        gnuplot.pdffigure(workdir .. 'ErrorTrain'..fileId..'.pdf')
        gnuplot.xlabel('Training Steps')
        gnuplot.ylabel('L2 Norm')
        x = torch.linspace(0,steps * epochs, steps*epochs/batchSize)
        gnuplot.plot(x, trainerror)
        gnuplot.plotflush()
    end
    
end
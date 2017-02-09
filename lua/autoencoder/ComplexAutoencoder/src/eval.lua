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

   
function evaluate(model, useCuda, trainerror, steps, workdir, fileId, testdata, epochs, batchSize, set)    
    if not epochs then epochs = 1 else epochs = epochs end
    
    if model then
        autoencoder = model
    end
    
    workdir = '/net/store/nbp/projects/phasesim/src_kstandvoss/lua/autoencoder/ComplexAutoencoder/Data/'
    --if not workdir then
     --   workdir = w
    --end
    
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
    --local activity = torch.Tensor(1,32,32):fill(1)
    activity = activity:cuda()
    local indim
    if set == 'LabelMe' then
      indim = 3
      act = function(activity, mode)
                    local result
                    if mode == 'plus' then 
                        result = activity:sub(1,3) + activity:sub(4,6)
                    elseif mode == 'minus' then
                        result = activity:sub(1,3) - activity:sub(4,6)
                    end 
                    return result
             end
    elseif set == 'MNIST' then
      indim = 1
      act = function(activity) return activity end
    elseif set == 'white' then
      indim = 1
      act = function(activity) return activity end 
    end     

    activities, phases = iterate(activity, autoencoder, indim, 20, set, useCuda and loaded)
    
    print('create Plot')
    if trainerror then
      torch.save(workdir .. 'Error' ..fileId .. '.dat', trainerror)
    end
    print(workdir .. 'Reconstructions'..fileId..'.png')
    image.save(workdir .. 'Reconstructions'..fileId..'.png', image.toDisplayTensor{input=activities:double(), scaleeach = true, padding=2})
    image.save(workdir .. 'Phases'..fileId..'.png', image.toDisplayTensor{input=phases:double(),scaleeach=true, padding=2,nrow=5})
  

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


function iterate(activity, autoencoder, indim, iterations, set, cuda)
    local activities = torch.Tensor(iterations/1 + 1, indim, activity:size()[2], activity:size()[3]):cuda()
    local phases = torch.Tensor(iterations/1, 3, activity:size()[2],activity:size()[3])
    local result = act(activity, 'minus') 
    activities[1] = result
    
    local phase = (torch.rand(#activity)*2*math.pi)-math.pi
    
    if cuda then
        activity = activity:cuda()
        activities = activities:cuda()
        phases = phases:cuda()
        phase = phase:cuda()
    end    
    
    print('start iterating')
    

    --[[ gnuplot.pdffigure('../Data/histInit.pdf')
     gnuplot.hist(phase)
     gnuplot.plotflush() ]]-- 
    local inp = {torch.cmul(activity,torch.cos(phase)), torch.cmul(activity,torch.sin(phase))}
    local out = autoencoder:forward(inp)
    local a_out = torch.sqrt(torch.pow(out[1],2) + torch.pow(out[2],2))    
    local phase_out = nn.Atan2():forward(out)--torch.atan2(out[2],out[1])

    
    j = 1
    for i = 1,iterations do
        -- Phase iterations 
        phase_out:add(torch.Tensor(#phase_out):normal(0,math.pi/8):cuda())
        out = autoencoder:forward({torch.cmul(activity,torch.cos(phase_out)),torch.cmul(activity,torch.sin(phase_out))})
        a_out = torch.sqrt(torch.pow(out[1],2) + torch.pow(out[2],2))  
        phase_out = nn.Atan2():forward(out)
        if i%1 == 0 then 
          local sin
          local cos
          local circMean
          sin = 0
          cos = 0
          for i=1,a_out:size(1) do 
              sin = sin + torch.cmul(a_out[i],torch.sin(phase_out[i])) 
              cos = cos + torch.cmul(a_out[i],torch.cos(phase_out[i])) 
          end  
          circMean = nn.Atan2():forward({sin,cos})  
          local hsv = torch.Tensor(3,activity:size(2),activity:size(3))
          hsv[1]:copy((circMean+math.pi)/(2*math.pi))
          hsv[2] = torch.Tensor(activity:size(2),activity:size(3)):fill(0.5)
          hsv[3] = torch.Tensor(activity:size(2),activity:size(3)):fill(0.5)
          rgb = image.hsv2rgb(hsv)
          if set ~= 'LabelMe' then
            rgb = torch.cmul(rgb, a_out:expand(3, a_out:size(2), a_out:size(3)):double())
          end  
          a = act(a_out, 'minus')

          phases[j] = rgb 
          activities[j+1] = a
          local phi = phase_out[torch.gt(activity,0)] -- find pixels of stimulus: activity>0  
          local var = 1 - (1/phi:numel())*torch.sqrt((torch.sum(torch.cos(phi))^2 + torch.sum(torch.sin(phi))^2))
          local allvar = 1 - (1/phase_out:numel())* torch.sqrt((torch.sum(torch.cos(phase_out))^2 + torch.sum(torch.sin(phase_out))^2))
          print(var/allvar)
          
          j = j+1
     
                 --torch.save('phases.dat',phase_out:double())
        --gnuplot.pdffigure('../Data/histPhase' .. j .. '.pdf')
        --gnuplot.hist(circMean)
        --gnuplot.hist(autoencoder:get(1).output)       
        --gnuplot.plotflush() 
        end
        
    end
   
    return activities, phases
    
end
function segIndex(testdata)
  
end


--reconstructionErrors = torch.Tensor(10)
--reconstructionErrorsNoPhase = torch.Tensor(10)
--circularVariances = torch.Tensor(10)
--synchronization = torch.Tensor(10)

--synchronization[j] = torch.abs(torch.sum(autoencoder:get(1).output[1]+autoencoder:get(1).output[2])/ torch.sum(torch.sqrt(torch.pow(autoencoder:get(1).output[1],2) + torch.pow(autoencoder:get(1).output[2],2))))       
--circularVariances[i+1] = var

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
--[[
gnuplot.pdffigure('../Data/histPhase' .. j .. '.pdf')
gnuplot.hist(phase_out)
gnuplot.plotflush() 
gnuplot.pdffigure('../Data/histHidPhase' .. j .. '.pdf')
gnuplot.hist(autoencoder:get(1).atan.output)
gnuplot.plotflush()]]--
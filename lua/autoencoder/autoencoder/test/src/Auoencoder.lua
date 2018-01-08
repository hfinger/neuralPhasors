require 'torch';
require 'unsup';
require 'nn';
require 'gnuplot';
require 'Atan2';

n_hidden = 10
n_input = 100
input = torch.rand(100,2)

require 'unsup'x = torch.cmul(input[{{},1}],torch.cos(input[{{},2}]))
y = torch.cmul(input[{{},1}],torch.sin(input[{{},2}]))
input = torch.cat(x,y,2)


phase = nn.Sequential()
phase:add(nn.SplitTable(2))
phasepar1 = nn.ParallelTable()
linear4 = nn.Linear(n_input,n_hidden)
linear5 = nn.Linear(n_input,n_hidden)
phasepar1:add()
phasepar1:add()
phasepar2 = nn.ParallelTable()
phasepar2:add(nn.Reshape(n_hidden,1))
phasepar2:add(nn.Reshape(n_hidden,1))
phase:add(phasepar1)
phase:add(phasepar2)
phase:add(nn.JoinTable(2))
phase:add(nn.Atan2())


--phase = nn.Atan2()
--check output dimensions of gradient 
res = phase:forward(input)
back = phase:backward(input,res)
print(res)
print(back)

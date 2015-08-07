require 'Atan2'
require 'nn'

enc = nn.Sequential()
encpar = nn.ParallelTable()

--imaginary activation part
img = nn.Sequential()

par1 = nn.Parallel(2,2)

linear1 = nn.Linear(n_input,n_hidden)
linear2 = nn.Linear(n_input,n_hidden)
linear1:share(linear2,'weight')

cont1 = nn.Sequential()
cont1:add(linear1) -- w*x
cont1:add(nn.Reshape(10,1)) --reshape for concat
cont2 = nn.Sequential()
cont2:add(linear2) -- w*y
cont2:add(nn.Reshape(10,1)) --reshape for concat
par1:add(cont1)
par1:add(cont2)

img:add(par1)
img:add(nn.Replicate(2)) --duplicate for phase calculations
img:add(nn.SplitTable(1))

partable = nn.ParallelTable()

par2 = nn.Parallel(2,2)
par2:add(nn.Square())
par2:add(nn.Square())

partable:add(par2)
partable:add(nn.Atan2())

img:add(partable)
print(img:forward(input))
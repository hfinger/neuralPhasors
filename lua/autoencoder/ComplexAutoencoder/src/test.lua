require 'unsup';
require 'nn';
require 'gnuplot';
require 'Encoder';
require 'Decoder';
require 'image';
require 'optim';
require 'functions'
require 'randomkit'
require 'pl'
require 'mattorch'
require 'mattorch'
local lfs = require"lfs";
function catch()
    require 'cunn';
    require 'cutorch';
end
if pcall(catch) then print('a') else print('b') end
--require 'cutorch'
--require 'ComplexCrit';

test1 = whiteStripes(5)
test2 = whiteStripes(5)
test3 = whiteStripes(5)

mattorch.save('output.mat', {w1 = test1, w2 = test2, w3 = test3})
--print(torch.sqrt( torch.pow(autoencoder.output[1],2) + torch.pow(autoencoder.output[2],2)))











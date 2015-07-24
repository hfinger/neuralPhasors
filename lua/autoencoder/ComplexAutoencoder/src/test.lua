require 'nn'
require 'pl'
opt = lapp[[
   -f,--full          (default 5000)        use the full dataset or only n samples
   -o,--optimization  (default "SGD")       optimization: SGD | LBFGS 
   -r,--learningRate  (default 0.05)        learning rate, for SGD only
   -b,--batchSize     (default 10)          batch size
   -m,--momentum      (default 0)           momentum, for SGD only
   -d,--weightDecay   (default 5e-5)        weight decay for SGD
   -h,--hidden        (default 50)          number of hidden units
   -k,--kernel        (default 7)           kernelsize
   -p,--phase         (default false)       training with/without phase
   -n,--noise         (default 0.1)         noise level for denoising         
   --coefL1           (default 0)           L1 penalty on the weights
   --coefL2           (default 0)           L2 penalty on the weights
]]
print(opt)

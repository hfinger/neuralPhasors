require 'torch'
require 'eval'

local traindata, testdata = loadData(true)
testdata = createData(100, testdata, 5)
evaluate(false, true, false, false , false, false, testdata)


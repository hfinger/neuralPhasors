classdef TiledConv < handle
  %TILEDCONV Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
%     idsW2Wfull
    idsW2WfullSplitted
    idsW2Winv
    
    tileSizeX
    tileSizeY
    fIn
    fOut
    
    useGpu
  end
  
  methods
    
    function obj = TiledConv(tileSizeX,tileSizeY,fIn,fOut,useGpu)
      % the individual patch size will be [tileSizeX+1 tileSizeY+1]
      obj.setDims(tileSizeX,tileSizeY,fIn,fOut);
      if nargin<5
          obj.useGpu=false;
      else
          obj.useGpu=useGpu;
      end
    end
    
    function setDims(this,tileSizeX,tileSizeY,fIn,fOut)
      this.tileSizeX = tileSizeX;
      this.tileSizeY = tileSizeY;
      this.fIn = fIn;
      this.fOut = fOut;
      
%       %% dummyWfull
%       dummyWfull = zeros( 2*tileSizeX, 2*tileSizeY, fIn, tileSizeX, tileSizeY, fOut );
%       dummyWfull(:) = 1:numel(dummyWfull);
%       dummyW = zeros( 1+tileSizeX, 1+tileSizeY, fIn, tileSizeX, tileSizeY, fOut );
%       for dxout=1:tileSizeX
%         for dyout=1:tileSizeY
%           dummyW(:,:,:,dxout,dyout,:) = dummyWfull(dxout:dxout+this.tileSizeX,dyout:dyout+tileSizeY,:,dxout,dyout,:);
%         end
%       end
%       this.idsW2Wfull = dummyW(:);
      
      %% dummyWfullSplited
      dummyWfull = zeros( tileSizeX, tileSizeY, fIn, 2, 2, tileSizeX, tileSizeY, fOut );
      dummyWfull(:) = 1:numel(dummyWfull);
      dummyW = zeros( 1+tileSizeX, 1+tileSizeY, fIn, tileSizeX, tileSizeY, fOut );
      for dxout=1:tileSizeX
        for dyout=1:tileSizeY
          dummyW(1:end-dxout,     1:end-dyout,     :,dxout,dyout,:) = dummyWfull(dxout:end, dyout:end, :,1,1,dxout,dyout,:);
          dummyW(end-dxout+1:end, 1:end-dyout,     :,dxout,dyout,:) = dummyWfull(1:dxout,   dyout:end, :,2,1,dxout,dyout,:);
          dummyW(1:end-dxout,     end-dyout+1:end, :,dxout,dyout,:) = dummyWfull(dxout:end, 1:dyout,   :,1,2,dxout,dyout,:);
          dummyW(end-dxout+1:end, end-dyout+1:end, :,dxout,dyout,:) = dummyWfull(1:dxout,   1:dyout,   :,2,2,dxout,dyout,:);
        end
      end
      this.idsW2WfullSplitted = dummyW(:);
      this.idsW2Winv = [];
      
      if this.useGpu
        this.idsW2WfullSplitted = gpuArray(this.idsW2WfullSplitted);
        this.idsW2Winv = gpuArray(this.idsW2Winv);
      end
    end
    
    function WfullSplitted = W2WfullSplitted(this,W)
      % W has dimensions [1+tileSizeX, 1+tileSizeY, fIn, tileSizeX, tileSizeY, fOut]
      % Wfull has dimensions [tileSizeX tileSizeY fIn 2 2 tileSizeX tileSizeY fOut]
      if this.useGpu
        WfullSplitted = parallel.gpu.GPUArray.zeros( [this.tileSizeX, this.tileSizeY, this.fIn, 2, 2, this.tileSizeX, this.tileSizeY, this.fOut] , 'double');
      else
        WfullSplitted = zeros( [this.tileSizeX, this.tileSizeY, this.fIn, 2, 2, this.tileSizeX, this.tileSizeY, this.fOut] , 'double');
      end
      WfullSplitted(this.idsW2WfullSplitted) = W(:);
    end
    
    function W = WfullSplitted2W(this,WfullSplitted)
      % Wfull has dimensions [tileSizeX tileSizeY fIn 2 2 tileSizeX tileSizeY fOut]
      % W has dimensions [1+tileSizeX, 1+tileSizeY, fIn, tileSizeX, tileSizeY, fOut]
      if this.useGpu
        W = parallel.gpu.GPUArray.zeros( [1+this.tileSizeX, 1+this.tileSizeY, this.fIn, this.tileSizeX, this.tileSizeY, this.fOut] , 'double' );
      else
        W = zeros( [1+this.tileSizeX, 1+this.tileSizeY, this.fIn, this.tileSizeX, this.tileSizeY, this.fOut] , 'double' );
      end
      W(:) = WfullSplitted(this.idsW2WfullSplitted);
    end
    
    function out = convWfullSplitted(this,WfullSplitted,in,type)
      % WfullSplitted has dimensions [tileSizeX tileSizeY fIn 2 2 tileSizeX tileSizeY fOut]
      % in has dimensions [X*tileSizeX Y*tileSizeY fIn]
      % type is either 'valid' (standard) or 'full'
      % out has dimensions [X-1 Y-1 tileSizeX tileSizeY fOut] or [X+1 Y+1 tileSizeX tileSizeY fOut] for full convolution
      
      if nargin<4
        type = 'valid';
      end
      
      numTilesX = size(in,1) / this.tileSizeX;
      numTilesY = size(in,2) / this.tileSizeY;
      
      WfullSplitted = reshape(WfullSplitted,[this.tileSizeX*this.tileSizeY*this.fIn 2*2*this.tileSizeX*this.tileSizeY*this.fOut]);
      
      in = reshape(in,[this.tileSizeX numTilesX this.tileSizeY numTilesY this.fIn]);
      if strcmp(type,'full')
        s_in = size(in);
        s_in(end+1:5) = 1;
        if this.useGpu
            tmp = parallel.gpu.GPUArray.zeros(s_in+[0 2 0 2 0],'double');
        else
            tmp = zeros(s_in+[0 2 0 2 0],'double');
        end
        tmp(:,2:end-1,:,2:end-1,:) = in;
        in = permute(tmp,[2 4 1 3 5]);
        numTilesX = numTilesX+2;
        numTilesY = numTilesY+2;
      else
        in = permute(in,[2 4 1 3 5]);
      end
      in = reshape(in,[numTilesX*numTilesY this.tileSizeX*this.tileSizeY*this.fIn]);
      
      out = in * WfullSplitted;
      out = reshape(out,[numTilesX numTilesY 2 2 this.tileSizeX this.tileSizeY this.fOut]);
      out = out(1:end-1,1:end-1,1,1,:,:,:) + out(2:end,1:end-1,2,1,:,:,:) + out(1:end-1,2:end,1,2,:,:,:) + out(2:end,2:end,2,2,:,:,:);
      out = reshape(out,[numTilesX-1 numTilesY-1 this.tileSizeX this.tileSizeY this.fOut]);
    end
    
    function out = convW(this,W,in,type)
      % W has dimensions [tileSizeX+1 tileSizeY+1 fIn tileSizeX tileSizeY fOut]
      assert(size(W,1)==this.tileSizeX+1)
      assert(size(W,2)==this.tileSizeY+1)
      assert(size(W,3)==this.fIn)
      assert(size(W,4)==this.tileSizeX)
      assert(size(W,5)==this.tileSizeY)
      assert(size(W,6)==this.fOut)
      
      WfullSplitted = this.W2WfullSplitted(W);
      out = this.convWfullSplitted(WfullSplitted,in,type);
    end
    
    
    function Wfull = backpropErrorWFull(this,in,outErr)
      % calc backprop Error for a given outErr=dPsi/d_out and a given input in.
      % Or in other words compute all pairwise products of input and output
      % and sum over them such that it results in a matrix with size of Wfull.
      % 
      %
      % in has dimensions [X*tileSizeX Y*tileSizeY fIn]
      % outErr has dimensions [X-1 Y-1 tileSizeX tileSizeY fOut]
      
      if this.useGpu
          Wfull = parallel.gpu.GPUArray.zeros([this.tileSizeX this.tileSizeY this.fIn 2 2 this.tileSizeX this.tileSizeY this.fOut],'double');
      else
          Wfull = zeros([this.tileSizeX this.tileSizeY this.fIn 2 2 this.tileSizeX this.tileSizeY this.fOut],'double');
      end
      
      X = size(in,1) / this.tileSizeX;
      Y = size(in,2) / this.tileSizeY;
      %% check input dimensions:
      assert(mod(X,1)==0)
      assert(mod(Y,1)==0)
      assert(size(in,3)==this.fIn)
      
      in = reshape(in,[this.tileSizeX X this.tileSizeY Y this.fIn]);
      in = permute(in,[1 3 5 2 4]);
      
      %% check outErr dimensions:
      assert(size(outErr,1)==(X-1))
      assert(size(outErr,2)==(Y-1))
      assert(size(outErr,3)==this.tileSizeX)
      assert(size(outErr,4)==this.tileSizeY)
      assert(size(outErr,5)==this.fOut)
      
      outErr = reshape(outErr, [(X-1)*(Y-1) this.tileSizeX*this.tileSizeY*this.fOut]);
      for dx=0:1
        for dy=0:1
          % in has dimensions [tileSizeX tileSizeY fIn X Y]
          inPart = in(:,:,:,dx+1:end-1+dx,dy+1:end-1+dy);
          inPart = reshape(inPart,[ this.tileSizeX*this.tileSizeY*this.fIn (X-1)*(Y-1)]);
          try
            Wfull(:,:,:,dx+1,dy+1,:,:,:) = Wfull(:,:,:,dx+1,dy+1,:,:,:) + reshape(inPart * outErr, [this.tileSizeX this.tileSizeY this.fIn 1 1 this.tileSizeX this.tileSizeY this.fOut]);
          catch err
            disp(err)
          end
        end
      end
      
      
    end
    
    function W = backpropErrorW(this,in,outErr)
      % output pairwise backprop sum(in*out)
      %
      % in has dimensions [X*tileSizeX Y*tileSizeY fIn]
      % outErr has dimensions [X-1 Y-1 tileSizeX tileSizeY fOut]
      % W has dimensions [tileSizeX+1 tileSizeY+1 fIn tileSizeX tileSizeY fOut]
      WfullSplitted = this.backpropErrorWFull(in,outErr);
      W = this.WfullSplitted2W(WfullSplitted);
    end
    

    function Winv = invertW(this,W)
      % W has dimensions [tileSizeX+1 tileSizeY+1 fIn tileSizeX tileSizeY fOut]

      if isempty(this.idsW2Winv)
        %% Calculate indexes for faster conversion:
        dummyWinv = zeros( this.tileSizeX+1, this.tileSizeY+1, this.fOut, this.tileSizeX, this.tileSizeY, this.fIn );
        dummyWinv(:) = 1:numel(dummyWinv);
        dummyWinv = permute(dummyWinv,[4 5 6 7 8 1 2 3]);
        dummyWFull = zeros( this.tileSizeX, this.tileSizeY, this.fIn, 2, 2, this.tileSizeX, this.tileSizeY, this.fOut );
        for dxout=1:this.tileSizeX
          for dyout=1:this.tileSizeY
            dummyWFull(dxout,dyout,:, 2, 2, dxout:end, dyout:end,:) = dummyWinv(dxout,dyout,:,1,1,1:end-dxout,     1:end-dyout,     :);
            dummyWFull(dxout,dyout,:, 2, 1, dxout:end, 1:dyout,  :) = dummyWinv(dxout,dyout,:,1,1,1:end-dxout,     end-dyout+1:end, :);
            dummyWFull(dxout,dyout,:, 1, 2, 1:dxout,   dyout:end,:) = dummyWinv(dxout,dyout,:,1,1,end-dxout+1:end, 1:end-dyout,     :);
            dummyWFull(dxout,dyout,:, 1, 1, 1:dxout,   1:dyout,  :) = dummyWinv(dxout,dyout,:,1,1,end-dxout+1:end, end-dyout+1:end, :);
          end
        end
        
        dummyW = dummyWFull(this.idsW2WfullSplitted);
        
        this.idsW2Winv = dummyW(:);
        
        if this.useGpu
            this.idsW2Winv = gpuArray(this.idsW2Winv);
        end
      
      end
      
      %% Now put together the inverse just by indexing:
      if this.useGpu
        Winv = parallel.gpu.GPUArray.zeros( this.tileSizeX+1, this.tileSizeY+1, this.fOut, this.tileSizeX, this.tileSizeY, this.fIn, 'double' );
      else
        Winv = zeros( this.tileSizeX+1, this.tileSizeY+1, this.fOut, this.tileSizeX, this.tileSizeY, this.fIn, 'double' );
      end
      Winv(this.idsW2Winv) = W(:);
      
    end
    
  end
  
  methods (Static)
    
    function performance()
      tileSizeX = 10;
      tileSizeY = 10;
      fIn = 100;
      fOut = 5;
      numTilesX = 10;
      numTilesY = 10;
      
      tiledConv = TiledConv(tileSizeX,tileSizeY,fIn,fOut);
      tiledConvInv = TiledConv(tileSizeX,tileSizeY,fOut,fIn);
      
      for i=1:20
        disp(i)
        W = randn(tileSizeX+1,tileSizeY+1,fIn,tileSizeX,tileSizeY,fOut);
        in = randn(numTilesX*tileSizeX,numTilesY*tileSizeY,fIn);
        out = tiledConv.convW(W,in,'valid'); 
        Wbackprop = tiledConv.backpropErrorW(in,out);
        Winv = tiledConv.invertW(W);
        out = permute(out,[3 1 4 2 5]);
        out = reshape(out,[size(out,1)*size(out,2) size(out,3)*size(out,4) size(out,5)]);
        reconstruction = tiledConvInv.convW(Winv,out,'valid'); 
        
      end
      
      
    end
    
    function unittest()
      
      %% a tiled Conv with just one weight element 1 and one input element 1 should output just one element 1
      tiledConv = TiledConv(10,10,3,2);
      W = zeros(11,11,3,10,10,2); 
      W(3,3,2,5,5,1)=1;
      in = zeros(40,40,3); 
      in(17,17,2)=1;
      out = tiledConv.convW(W,in,'full'); 
      assert(sum(out(:)==1)==1) 
      
      %% check inversion:
      tiledConvInv = TiledConv(10,10,2,3);
      Winv = tiledConv.invertW(W);
      assert(numel(unique(tiledConv.idsW2Winv)) == numel(tiledConv.idsW2Winv));
      out = permute(out,[3 1 4 2 5]);
      out = reshape(out,[size(out,1)*size(out,2) size(out,3)*size(out,4) size(out,5)]);
      reconstr = tiledConvInv.convW(Winv,out,'valid'); 
      assert(sum(reconstr(:)==1)==1) 

      %% check if reconstruction is the same as the input
      reconstr = permute(reconstr,[3 1 4 2 5]);
      reconstr = reshape(reconstr,[size(reconstr,1)*size(reconstr,2) size(reconstr,3)*size(reconstr,4) size(reconstr,5)]);
      assert(isequal(reconstr,in));
      
      for i=1:1000
        TiledConv.unittestReconstr(randi(4,1),randi(4,1),randi(4,1),randi(4,1));
      end
      
    end
    
    function unittestReconstr(tileSizeX,tileSizeY,fIn,fOut)
%       disp([tileSizeX,tileSizeY,fIn,fOut])
      
      %% a tiled Conv with just one weight element 1 and one input element 1 should output just one element 1
      tiledConv = TiledConv(tileSizeX,tileSizeY,fIn,fOut);
      
      W = zeros(tileSizeX+1,tileSizeY+1,fIn,tileSizeX,tileSizeY,fOut);
      Wdim1=randi(size(W,1),1);
      Wdim2=randi(size(W,2),1);
      Wdim3=randi(size(W,3),1);
      
      numTilesX = randi(3,1);
      numTilesY = randi(3,1);
      in = zeros(numTilesX*tileSizeX,numTilesY*tileSizeY,fIn);
      inPosInTileX = randi(tileSizeX,1);
      inPosInTileY = randi(tileSizeY,1);
      inDim1 = (randi(numTilesX,1)-1) * tileSizeX + inPosInTileX;
      inDim2 = (randi(numTilesY,1)-1) * tileSizeY + inPosInTileY;
      inDim3 = Wdim3;
      in(inDim1,inDim2,inDim3) = 1;
      
      Wdim4 = -Wdim1 + inPosInTileX + 1;
      Wdim5 = -Wdim2 + inPosInTileY + 1;
      Wdim6=randi(size(W,6),1);
      
      Wdim4 = mod(Wdim4-1,size(W,4)) + 1;
      Wdim5 = mod(Wdim5-1,size(W,5)) + 1;
      W(Wdim1,Wdim2,Wdim3,Wdim4,Wdim5,Wdim6) = 1;
      
      out = tiledConv.convW(W,in,'full'); 
      
      %% check 
      assert(sum(out(:)==1)==1)
      
      %% check inversion:
      tiledConvInv = TiledConv(tileSizeX,tileSizeY,fOut,fIn);
      Winv = tiledConv.invertW(W);
      assert(numel(unique(tiledConv.idsW2Winv)) == numel(tiledConv.idsW2Winv));
      out = permute(out,[3 1 4 2 5]);
      out = reshape(out,[size(out,1)*size(out,2) size(out,3)*size(out,4) size(out,5)]);
      reconstr = tiledConvInv.convW(Winv,out,'valid'); 
      assert(sum(reconstr(:)==1)==1) 

      %% check if reconstruction is the same as the input
      reconstr = permute(reconstr,[3 1 4 2 5]);
      reconstr = reshape(reconstr,[size(reconstr,1)*size(reconstr,2) size(reconstr,3)*size(reconstr,4) size(reconstr,5)]);
      assert(isequal(reconstr,in));
      
      %% check backprop error / try to reconstruct W from in and out:
      
      in2 = zeros([size(in,1)+tileSizeX*2 size(in,2)+tileSizeY*2 size(in,3)]);
      in2(tileSizeX+1:end-tileSizeX,tileSizeY+1:end-tileSizeY,:) = in;
      out2 = tiledConv.convW(W,in2,'valid'); 
      Wbackprop = tiledConv.backpropErrorW(in2,out2);
      assert(isequal(W,Wbackprop))
    end
    
  end
  
end


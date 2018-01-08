function [ distance_matrix ] = calculateDistanceMatrix( voxel_coords, type, value )
    % This function calculates the distance matrix
    % INPUT:
    %       voxel_coords: data points with 3 dimensions
    %       type: 0 for t-nn calculation, 1 for calculation with threshold
    %       value: num of neighbors for type 0, distance for type 1

    if type == 1
        d=value;
        M=size(voxel_coords,1); N=3;
        X=voxel_coords(:,1);
        Y=voxel_coords(:,2);
        Z=voxel_coords(:,3);

        II=cell(M,1);
        JJ=cell(M,1);
        SS=cell(M,1);

        for jj=1:M
            sqdists= (X(jj)-X).^2 + (Y(jj)-Y).^2 + (Z(jj)-Z).^2;
            II{jj}=find( sqdists<=d).';
            JJ{jj}=II{jj};
            JJ{jj}(:)=jj;
            SS{jj}=sqdists(II{jj}).';
        end

        V_d=sqrt(sparse([II{:}],[JJ{:}],[SS{:}],M,M));
        distance_matrix = V_d;
        filename = ['../distances/' int2str(d) '_distance_matrix.mat'];
        save(filename, 'distance_matrix');
    elseif type == 0
        gen_nn_distance(voxel_coords, value, 10, 0);
    end
end

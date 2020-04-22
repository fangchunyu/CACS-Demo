function mat_rs = resize2D(mat, X, Y)
    % extract the original number of pixels from the size of the matrix
    [Nx_input, Ny_input] = size(mat);
    
    interp_mode = '*bicubic';
    
    % extract the desired number of pixels
    Nx_output = X;
    Ny_output = Y;

    % create normalised plaid grids of current discretisation
    [x_mat, y_mat] = ndgrid((0:Nx_input-1)/(Nx_input-1), (0:Ny_input-1)/(Ny_input-1));       
    %[x_mat, y_mat, z_mat] = ndgrid(gpuArray((0:Nx_input-1)/(Nx_input-1)), gpuArray((0:Ny_input-1)/(Ny_input-1)), gpuArray((0:Nz_input-1)/(Nz_input-1)));

    % create plaid grids of desired discretisation
    [x_mat_interp, y_mat_interp] = ndgrid((0:Nx_output-1)/(Nx_output-1), (0:Ny_output-1)/(Ny_output-1));
    %[x_mat_interp, y_mat_interp, z_mat_interp] = ndgrid(gpuArray((0:Nx_output-1)/(Nx_output-1)), gpuArray((0:Ny_output-1)/(Ny_output-1)), gpuArray((0:Nz_output-1)/(Nz_output-1)));
%{
    x_mat = gather(x_mat);
    y_mat = gather(y_mat);
    z_mat = gather(z_mat);
    x_mat_interp = gather(x_mat_interp);
    y_mat_interp = gather(y_mat_interp);
    z_mat_interp = gather(z_mat_interp);
    %}
    % compute interpolation; for a matrix indexed as [M, N, P], the
    % axis variables must be given in the order N, M, P
    %mat_rs = interp3_gpu(y_mat(1,:,1), x_mat(:,1,1), z_mat(1,1,:), mat, y_mat_interp, x_mat_interp, z_mat_interp);
    mat_rs = interp2(y_mat, x_mat, mat, y_mat_interp, x_mat_interp, interp_mode);        
    %matGPU = gpuArray(mat);
    %mat_rs = double(gather(interp3(y_mat, x_mat, z_mat, matGPU, y_mat_interp, x_mat_interp, z_mat_interp, interp_mode)));
end

function N_loads = calculate_N_loads(y1_n, y2_n, P, T, LN_Pi, Ene, state_rd, NumOfCom1, NumOfCom2, Cell_a, Cell_b, Cell_c)
    
    % Calculate the N_loads array based on the input parameters.
    P_ref_tot = 0.2 ; % bar
    T_ref = 300 ; % K

    % Supercell volume calculation
    Supercell = Cell_a * Cell_b * Cell_c;
    
    % Initialize N_loads array
    num_cases = length(y1_n); % Assuming y1_n, y2_n, P, and T are all vectors of same length
    N_loads = zeros(num_cases, 2);
    
    y1_n(y1_n == 0) = 1e-10;
    y2_n(y2_n == 0) = 1e-10;
    
    neg_idx = y2_n < 0; 
    y2_n(neg_idx) = 1e-10;
    y1_n(neg_idx) = 1; 

    % Reference molar ratio
    y1 = 0.5; 
    y2 = 0.5;  
    
    % Compute pressure and temperature ratios
    P_rat_1 = (P .* y1_n) / (P_ref_tot * y1);  % Element-wise multiplication (.*)
    P_rat_2 = (P .* y2_n) / (P_ref_tot * y2);  % Element-wise multiplication (.*)
    beta_new = 1 ./ T;  % Element-wise division (./)
    beta_ref = 1 / T_ref;
    
    % Initialize variables
    adder = zeros(num_cases, 2);  % To accumulate values for each case
   
    [i_grid, j_grid] = ndgrid(0:NumOfCom1, 0:NumOfCom2);
    
    i_grid = repmat(i_grid, 1, 1, num_cases);  % (NumOfCom1+1, NumOfCom2+1, num_cases)
    j_grid = repmat(j_grid, 1, 1, num_cases);
    
    LN_Pi_3D = repmat(LN_Pi, 1, 1, num_cases);   % (NumOfCom1+1, NumOfCom2+1, num_cases)
    Ene_3D = repmat(Ene, 1, 1, num_cases);       % (NumOfCom1+1, NumOfCom2+1, num_cases)
    
    P_rat_1_3D = reshape(P_rat_1, 1, 1, num_cases);
    P_rat_2_3D = reshape(P_rat_2, 1, 1, num_cases);
    T_3D = reshape(T, 1, 1, num_cases);
    beta_new_3D = reshape(beta_new, 1, 1, num_cases);
    
    LN_Pi_Rew = LN_Pi_3D + ...
                i_grid .* log(P_rat_1_3D) + ...
                j_grid .* log(P_rat_2_3D) + ...
                (i_grid + j_grid) .* log(T_ref ./ T_3D) - ...
                Ene_3D .* (beta_new_3D - beta_ref);
    
    % Compute expectation values (no need for loop over cases)
    exp_vals = exp(LN_Pi_Rew);  % Exponentiate the entire matrix at once
    
    valid_mask = state_rd > 0;
    
    % Apply the mask to exp_vals
    exp_vals_valid = exp_vals .* valid_mask;  
    
    % Compute weighted sums (element-wise multiplication and summation over the first two dimensions)
    adder(:,1) = squeeze(sum(sum(exp_vals_valid .* i_grid, 1), 2)); 
    adder(:,2) = squeeze(sum(sum(exp_vals_valid .* j_grid, 1), 2)); 

    % Compute nom (sum over the probability distribution)
    nom = squeeze(sum(sum(exp_vals_valid, 1), 2));

    % Compute final N_loads values
    N_loads(:, 1) = adder(:, 1) ./ nom / Supercell;
    N_loads(:, 2) = adder(:, 2) ./ nom / Supercell;

    N_loads(y1_n == 1e-10, 1) = 0;
    N_loads(y2_n == 1e-10, 2) = 0;

end
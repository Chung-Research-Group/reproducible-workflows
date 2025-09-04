function partial_loadings = IAST_func(num_components, Pressure, isotherm_params_array, gas_phase_mol_fraction, T, initial_guess)
    
    %% Parameter Unpacking and Decalarations
    size_of_mol_frac_vector = length(gas_phase_mol_fraction(:, 1));     % Length of mole fraction vector
    partial_pressures = Pressure .* gas_phase_mol_fraction;
    adsorbed_mole_fractions = zeros(size_of_mol_frac_vector, num_components);
    %
    %% Initialization of Guess for IAST Equations Solution
    % Solve the non-linear equations returened by the 'func' functions using fsolve method.
    % If Initial Guess has not been provided in input arguments, 
    % it will be guessed proportional to pure compoenent loading.
    
    if nargin < 6 || ~(any(sum(initial_guess, 2)))
        loading_array = Isotherm_functions(isotherm_params_array, partial_pressures, T, 1);
        loading_array(loading_array == 0) = 1e-07;
        initial_guess = loading_array ./ (sum(loading_array, 2));
    end
    

    tolerance = 1e-10;  % Use to avoid undefined values when the column is filled with only non-adsorbing gas
    guess = initial_guess(:, 1:num_components-1);
    

%     b = ones(num_components-1, 1);
%     jp = spdiags([b , b], 0:1, num_components-1, num_components);
%     jp = full(jp);
%     jp(end, :) = 1;
%     jp = sparse(jp);

%     options = optimset('Display','off', 'TolFun', 1e-10, 'MaxIter', 400, 'JacobPattern', jp);   % Options struc for fsolve
    options = optimset('TolFun', 1e-05, 'Display', 'off', 'Diagnostics', 'off', 'TolX', 1e-09, 'FinDiffRelStep', 1e-9, 'Algorithm','levenberg-marquardt');
    for k = 1:size_of_mol_frac_vector
        try
            if sum(loading_array(k, :), 2) - tolerance < 0
                adsorbed_mole_fractions(k, :) = zeros(1, num_components)+tolerance;

            % In case of single adsorbing component, the initial guess will give the correct loading
            elseif num_components == 1
                adsorbed_mole_fractions(k, :) = initial_guess(k, :);
            else

                
                function_handle = @(x)func(x, k);
                [solution, ~, exitflag] = fsolve(function_handle, guess(k, :), options);

                if exitflag == 0
                    error("IAST Error: Maximum Iterations for fsolve reached! Terminating the solution.");
                elseif exitflag < 0
                    error("IAST Error: Failed to solve the Equation!")
                end
                adsorbed_mole_fractions(k, :) = [solution, 1 - sum(solution, 2)];

                ads_mol_frac_unity_tol = 1e-10;    % Tolerance value to check that sum of adsorbed mole fraction is unity
                if ~(ismembertol(sum(adsorbed_mole_fractions(k, :)), 1, ads_mol_frac_unity_tol))
                    error("IAST Error: The sum of adsorbed mole fractions is not unity. Use better initial guess!")
                end
            end
        catch ME1
                error(strcat('IAST Error: Oh error occured!', ...
                            'refer below to see what is casuing the trouble. ', ' Error: ', ME1.message));    
        end
    end

    pressure_0 = partial_pressures ./ adsorbed_mole_fractions;
    loading_array = Isotherm_functions(isotherm_params_array, pressure_0, T, 1);
    inverse_loading = sum(adsorbed_mole_fractions./loading_array, 2);
    loading_total = 1 ./ inverse_loading;
    partial_loadings = (adsorbed_mole_fractions .* loading_total);

    function residual = func(adsorbed_MF, k)
        % Calculate the difference in spreading pressure for (N-1) components.
        % The spreading pressures are all calculated at the fictious pressure calculated based on the
        % adsorbed mole fraction and partial pressure of the components.
        % Fictitous Pressure = p_* = p_i / x_i
        
        s = 1e-05;       
        x = [(adsorbed_MF+s), (1-sum(adsorbed_MF)+s)];
        spreading_pressures = Isotherm_functions(isotherm_params_array, partial_pressures(k, :)./x, T(k, 1));
        
        spreading_pressures_diff = spreading_pressures(:, 1:num_components-1) - spreading_pressures(:, 2:num_components);
        if ~isreal(spreading_pressures_diff)
            residual = 1e05;
        else
            residual = spreading_pressures_diff;
        end
        if any(isnan(residual))
            ss =1;
        end
    end
end
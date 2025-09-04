function solution = Isotherm_functions(iso_params, P, T, calc_mode)
%% This functions calculates and returns:
% 1: Pure compoenent loading for the given pressure and temperature.
% Pure component reduced grand potential the given pressure and
% temperature.

    iso_type_flag = iso_params(1, :);
    R = 8.314;
    if (nargin<4)
        calc_mode = 0;
    end
    iso_params(6, :) = 1./iso_params(6, :);     % Transforming n_b because Sunghyun uses 1/n as power in his formula 
    iso_params(7, :) = 1./iso_params(7, :);     % Transforming n_d because Sunghyun uses 1/n as power in his formula
    %% Pure Loading
    if calc_mode == 1
        pure_loadings = zeros(size(P));
        for i=1:length(iso_params(1, :))    % Loop for all the components

            % DS-Langmuir-Freunclich Case
            b = iso_params(4, i) .* (exp(-iso_params(8, i)/R.*(1./T - 1/298))).^iso_params(6, i);
            d = iso_params(5, i) .* (exp(-iso_params(8, i)/R.*(1./T - 1/298))).^iso_params(7, i);    
            
            pure_loadings(:, i) = iso_params(2, i).*b.*P(:, i).^(iso_params(6, i)) ./ (1+b.*P(:, i).^(iso_params(6, i)))...
                                + iso_params(3, i).*d.*P(:, i).^(iso_params(7, i)) ./ (1+d.*P(:, i).^(iso_params(7, i)));
           
        end
        solution = pure_loadings;
    
    %% Reduced Potential Calculation    
    else
        reduced_grand_potential = zeros(size(P));
        for i=1:length(iso_params(1, :))    % Loop for all the components
        
            % DS-Langmuir-Freunclich Case             
            b = iso_params(4, i) .* exp(-iso_params(8, i)/R.*(1./T - 1/298)).^iso_params(6, i);
            d = iso_params(5, i) .* exp(-iso_params(8, i)/R.*(1./T - 1/298)).^iso_params(7, i);    
            
            reduced_grand_potential(:, i) = iso_params(2, i)/iso_params(6, i) .* log(1 + b.*P(:, i).^iso_params(6, i)) ...
                                          + iso_params(3, i)/iso_params(7, i) .* log(1 + d.*P(:, i).^iso_params(7, i));
        
        end
        solution = reduced_grand_potential; 
    end
end
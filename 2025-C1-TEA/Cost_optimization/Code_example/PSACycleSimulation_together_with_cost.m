function [ objectives, constraints ] = PSACycleSimulation_together_with_cost( x, type, N)

    % % Retrieve required process variables provided as inputs of the function
    % L           = process_vars(1)       ;   % Length of the column [m]   
    % P_0         = process_vars(2)       ;   % Adsorption pressure [Pa]
    % ndot_0      = process_vars(3)       ;   % Inlet molar flux [mol/s/m^2]
    % t_ads       = process_vars(4)       ;   % Time of adsorption step [s]
    % alpha       = process_vars(5)       ;   % Light product reflux ratio [-]
    % beta        = process_vars(6)       ;   % Heavy product reflux ratio [-]
    % P_I         = process_vars(7)       ;   % Intermediate pressure [Pa]
	% P_l         = process_vars(8)       ;   % Purge Pressure [Pa]
    % t_depres    = process_vars(9)       ;   % Time of depressurization [s]
    % t_pres      = process_vars(10)      ;   % Time of pressurization [s]

    process_variables = [x(9), x(1), x(1)*x(4)/8.314/298.15,  x(2),  x(3),  x(5), 1e4, x(6), x(7), x(8)] ;
	
	try
    [objectives, constraints] = PSACycle_together_with_cost(process_variables, [], type, N) ;
	catch
    %warning('Problem using function.  Assigning a value of 0 for objectives and constraints vialations');
    objectives  = 1e5 ;
	constraints(1) = 1 ;
	constraints(2) = 1 ;
	%constraints(3) = 1 ;
	end
	
end
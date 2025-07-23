function InputParams = ProcessInputParameters_together_with_cost(process_vars, N)
%InputParameters: Specify input parameters for PSA simulation
%   Number of finite volumes is always an input from the calling script/
%   function. In addition, x is used to allow varying specific parameters 
%   (times, pressures, heats of adsorptions, etc.) from the calling script/
%   function. It is not necessary for simple simulations.
    
%% State all input parameters for the simulation
    
    % Retrieve required process variables provided as inputs of the function
    L           = process_vars(1)       ;   % Length of the column [m]   
    P_0         = process_vars(2)       ;   % Adsorption pressure [Pa]
    ndot_0      = process_vars(3)       ;   % Inlet molar flux [mol/s/m^2]
    t_ads       = process_vars(4)       ;   % Time of adsorption step [s]
    alpha       = process_vars(5)       ;   % Light product reflux ratio [-]
    beta        = process_vars(6)       ;   % Heavy product reflux ratio [-]
    P_I         = process_vars(7)       ;   % Intermediate pressure [Pa]
	P_l         = process_vars(8)       ;   % Purge Pressure [Pa]
    
    % Operating bed parameters
    t_pres      = process_vars(10)      ;   % Ttime of pressurization step [s]
    t_CnCdepres = process_vars(9)       ;   % Time of depressurization step [s]
    t_CoCdepres = 70                    ;   % Maximum/time of depressurization step [s]
    t_LR        = t_ads                 ;   % Time of light reflux step [s]
    t_HR        = t_LR                  ;   % Time of heavy reflux step [s]
    tau         = 0.5                   ;   % Parameter used for determining speed of pressure change
    P_inlet     = 1.02                  ;   % Pressure of feed gas at the inlet of the adsorption step
    
    % Flue gas parameters and constants
    R          = 8.314                  ;   % Universal gas constant [J/mol/K : Pa*m^3/mol/K]
    T_0        = 298.15                 ;   % Feed temperature of flue gas [K]
    y_0        = 0.10                   ;   % Inlet gas CO2 mole fraction[-]
    Ctot_0     = P_0/R/T_0              ;   % Inlet total concentration [mol/m^3]
    v_0        = ndot_0/Ctot_0          ;   % Inlet velocity and scaling parameter [m/s]
    mu         = 1.61e-5                ;   % Viscosity of gas [Pa*s]
    epsilon    = 0.40                   ;   % Void fraction
    D_m        = 1.2995e-5              ;   % Molecular diffusivity [m^2/s]
    K_z        = 0.09                   ;   % Thermal conduction in gas phase [W/m/k]
    C_pg       = 32.8                   ;   % Specific heat of gas [J/mol/k]
    C_pa       = 32.8                   ;   % Specific heat of adsorbed phase [J/mol/k]
    MW_CO      = 0.028                  ;   % Molecular weight of CO [kg/mol]
    MW_CO2     = 0.044                  ;   % Molecular weight of CO2 [kg/mol]
	%feed_gas  = 'Constant Pressure'    ;   % Whether flue gas during the feed step has a constant pressure or velocity
    feed_gas   = 'Constant Velocity'    ;   % Whether flue gas during the feed step has a constant pressure or velocity
    
    % Adsorbent parameters
    ro_s        = 2162.5                ;   % Density of the adsorbent [kg/m^3]
    r_p         = 5e-3                  ;   % Radius of the pellets [m]
    C_ps        = 900                   ;   % Specific heat capacity of the adsorbent [J/kg/K]
    q_s         = 4.06                  ;   % Molar loading scaling factor [mol/kg]
    q_s0        = q_s*ro_s              ;   % Molar loading scaling factor [mol/m^3]
    k_CO_LDF    = 0.1500                ;   % Mass transfer coefficient for CO [1/s]
    k_CO2_LDF   = 0.1000                ;   % Mass transfer coefficient for CO2 [1/s]
    
    % Isotherm parameters

    q_s_b      = [0.41, 2.98]   ;   % Saturation loading on site b [mol/kg]
    q_s_d      = [3.79,  0]   ;   % Saturation loading on site d [mol/kg]
    b          = [2.61e-10,  2.61e-09]   ;   % Pre-exponential factor for site b [Pa-1]
    d          = [5.58e-11,  0]  ;   % Pre-exponential factor for site d [Pa-1]
    deltaU_b   = [-39964,  -17994]  ;   % Heat of adsorption for site b [J/mol]
    deltaU_d   = [-34895,  0]  ;   % Heat of adsorption for site d [J/mol]

    deltaU     = [-39129.19601, -18547.63072]    ; % Heat of adsorption for high (n_g)   
    
%% Distribute the values to the necessary variables
    Params     = zeros(39, 1) ;
    Params(1)  = N			  ;
    Params(2)  = deltaU(1)    ;
    Params(3)  = deltaU(2)    ;
    Params(4)  = ro_s		  ;
    Params(5)  = T_0		  ;
    Params(6)  = epsilon	  ;
    Params(7)  = r_p		  ;
    Params(8)  = mu			  ;
    Params(9)  = R			  ;
    Params(10) = v_0		  ;
    Params(11) = q_s0		  ;
    Params(12) = C_pg		  ;
    Params(13) = C_pa		  ;
    Params(14) = C_ps		  ;
    Params(15) = D_m		  ;
    Params(16) = K_z		  ;
    Params(17) = P_0		  ;
    Params(18) = L			  ;
    Params(19) = MW_CO		  ;
    Params(20) = MW_CO2		  ;
    Params(21) = k_CO_LDF	  ;
    Params(22) = k_CO2_LDF	  ;
    Params(23) = y_0		  ;
    Params(24) = tau		  ;
    Params(25) = P_l		  ;
    Params(26) = P_inlet	  ;
    Params(27) = 1			  ;   % Place for y at outlet of Adsorption = y at inlet of Light Reflux: y_LP
                                  % y_LR = 1 - No initial guess for inlet CO2 mole fraction in Light Reflux step
    Params(28) = 1			  ;   % Place for T at outlet of Adsorption = T at inlet of Light Reflux: T_LP
                                  % T_LR = 1 - No initial guess for inlet temperature in Light Reflux step
    Params(29) = 1			  ;   % Place for ndot at outlet of Adsorption = ndot at inlet of Light Reflux
                                  % ndot_LR = 1 - No initial guess for inlet ndotin Light Reflux step
    Params(30) = alpha    	  ;
    Params(31) = beta         ;
    Params(32) = P_I          ;
    Params(33) = y_0          ;   % Place for y at outlet of CnC depressurization = y at inlet of of Heavy Reflux: y_HP
                                  % y_HR = y_0 - Initial guess for inlet CO2 mole fraction in Heavy Reflux step
    Params(34) = T_0          ;   % Place for T at outlet of CnC depressurization = T at inlet of of Heavy Reflux: T_HP
                                  % T_HR = T_0 - Initial guess for inlet temperature in Heavy Reflux step
    Params(35) = ndot_0*beta  ; % 300/30   % Place for ndot at outlet of CnC depressurization = ndot at inlet of Heavy Reflux
                                  % ndot_HR = ndot_0*beta - Initial guess for inlet ndotin Heavy Reflux step
    Params(36) = 0.01    	  ;   % Place for y at outlet of Adsorption = y at inlet of CnC pressurization: y_LP
                                  % y_LR = 0.01 - Initial guess for inlet CO2 mole fraction in CnC pressurization step
    Params(37) = T_0    	  ;   % Place for T at outlet of Adsorption = T at inlet of CnC pressurization: T_LP
                                  % T_LR = T_0 - Initial guess for inlet temperature in CnC pressurization step
    Params(38) = ndot_0  	  ;   % Place for ndot at outlet of Adsorption = ndot at inlet of CnC pressurization 
                                  % NOTE: not used, seems not necessary. ndot_LR = ndot_0 - Initial guess for inlet ndot
                                  % in CnC pressurization step
    
    if strcmpi(feed_gas, 'Constant Pressure') == 1
        Params(end) = 1 ;
    elseif strcmpi(feed_gas, 'Constant Velocity') == 1
        Params(end) = 0 ;
    else
        error('Please specify whether inlet velocity or pressure is constant for the feed step')
    end
	
    Times          = [ t_pres; t_ads; t_CnCdepres; t_LR; t_CoCdepres; t_HR ] ;
 
    IsothermParams = [q_s_b, q_s_d, b, d, deltaU_b, deltaU_d, 0] ;
%   

%% Economic Parameters
    
    total_CO_flow                         = 4129.68      ;   % molar flow rate of CO in the flue gas [mol/s]
    electricity_cost                      = 0.07         ;   % Cost of electricity [$/kWh] (Turton 2018)
    CEPCI                                 = 798.8        ;   % CEPCI of present year (2024).
 
    EconomicParams    = zeros(3, 1)              ;
    EconomicParams(1) = total_CO_flow            ;
    EconomicParams(2) = electricity_cost         ;
    EconomicParams(3) = CEPCI                    ;
%   
%% Combine all lists into one variable of cells that can easily be passed
    InputParams{1} = Params          ;
    InputParams{2} = IsothermParams  ;
    InputParams{3} = Times           ;
    InputParams{4} = EconomicParams  ;
%   
end 

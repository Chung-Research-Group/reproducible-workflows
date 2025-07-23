function  [objectives, constraints, a, b, c, d, e, t1, t2, t3, t4, t5] = PSACycle_together_with_cost(vars, x0, type, N)
%Skarstrom: Simulate a 5-Step Modified Skarstrom PSA cycle
%   This function is able to simulate a 5-Step Modified Skarstrom PSA cycle
%   and provide the state variables and process objectives. Finite Volume 
%   method is used for calculating the derivatives due to the inherent mass
%   and energy conservation due to calculating the flux. An example of the 
%   FVM set up is shown below with the ordering of the volumes 
%   
%   Diagram of Column divided into sections shown below with corresponding
%   designation for inlet and outlet. N is a user defined value below
%   
%   Inlet                                                     Outlet
%      -------------------------------\    \----------------------
%       1   |   2   |   3   |   4   |           N|   N+1   |   N+2
%      --------------------------------\    \---------------------
%   
%   Input:
%   vars: Process variables which are in order: Length of column [m] 
%         adsoprtion pressure [Pa], inlet molar flux [mol/s/m^2], time of 
%         adsorption step [s], light product reflux ratio [-], heavy product
%         reflux ratio [-], and intermediate pressure. Not all of these 
%         variables are used for all cycles, so when the simulated cycle do
%         not have the step relating the "vars", assign these to zero. These
%         variables are inputed as they are desing variables and can be 
%         changed to be optimized.
%   
%   x0  : Initial profile of state variables in the column, it's not mandatory
%   N   : Number of Finite Volumes, it's not mandatory
%   
%   Output:
%   The output variables can be customized depend upon of the purpose. Here
%   are two cases:
%   1) Optimization: it is necessary to get the "objectives" variable: this
%      provides the purity and recovery of CO2. And also the constraints: 
%      Provides whether a constraint has been violated or not. There are 
%      some constraints that need to be satisfied. First, the recovery of 
%      the process must be over 90%. Second, the purity of the final product
%      must be greater than the entering stream.
%   
%      a-e: These are the state variables for the five steps: 
%      CoCPressurization, Adsorption, Heavy Reflux, CnCDepressurization and
%      Light Reflux.
%      The matrix is set up where each row contains all the values of the
%      state variables at a single moment in time, and each column contains
%      the values of a single state variable at a single location throughout
%      the step. As shown above in the graph, for each state variable there
%      are N volumes for which the value is known at the center of volume.
%      In addition, for calculating purity and recovery, the state variable
%      values at the two ends of the column are also provided. NOTE: since
%      there are no spatial derivatives in the solid loading equations, and
%      purity and recovery do not need these values, they are assumed to be
%      equal to the volume next to it. The order of the state variables are
%      as follows
%   
%      a-e(:, 1:N+2) are the dimensionless pressure. In order to retrieve the
%      true pressure, multiply the dimensionless pressure by the adsorption
%      pressure, P_0 or P_H.
%      a-e(:,  N+3:2*N+4) are the CO gas mole fraction
%      a-e(:, 2*N+5:3*N+6) are the dimensionless CO2 molar loadings. In order
%      to retrieve the true molar loading, multiply the dimensionless molar
%      loading by the molar loading scaling factor, q_s
%      a-e(:, 3*N+7:4*N+8) are the dimensionless N2 molar loadings. In order
%      to retrieve the true molar loading, multiply the dimensionless molar
%      loading by the molar loading scaling factor, q_s
%      a-e(:, 4*N+9:5*N+10) are the dimensionless column temperature. In order
%      to retrieve the true column temperature, multiply the column temperature
%      by the feed temperature, T_0
%   
%      t1-t5 are the dimensionless times for the five steps. In order to
%      retrieve the true time, multiply the dimensionless time by the ratio
%      of the length to the velocity scaling factor (L/v_0)
%   
%   2) Data collection intended to ANNs training: To get an output suitable
%      to train Artificial Neural Network of each step of the PSA cycle,
%      it is necessay to collect the initial and final states of variables
%      ans store it in adequate variables to be used by the ANN toolbox. 
%      These variables are: a_fin, b_fin, c_fin, d_fin, e_fin, a_in, b_in,
%      c_in, d_in, e_in.
%      Where fin means final state, and in means initial state.
%   
%   The following assumptions were made:
%   1) Ideal Gas Law is used to describe the gas phase
%   2) No concentration, pressure, or temperature gradient in the the radial
%      or azmuth directions
%   3) Linear Driving Force is used to describe the gas diffusion into the
%      adsorbent
%   4) Adsorbent properties, and void fraction are constant throughout the
%      column
%   5) Viscosity of the gas is independent of pressure
%   6) There is thermal equilibrium between the adsorbent and the gas phase
%   7) Ergun Equation is used to describe the pressure drop across the bed 
%   8) Column operates adiabatically, so any wall energy balance have to be
%      used
%   9) An axially dispersed plug flow model is used to represent bulk fluid
%      flow
%   
%% Initialize parameter for the simulation
%   Specify the parameters for the simulation
    
    % Initialize objectives and constraints output
    switch type
        case 'ProcessEvaluation'
            constraints = [0, 0];
        case 'EconomicEvaluation'
            constraints = [0, 0];
        otherwise
            error('Error. %s is not a recognizable type of operation.' , type) ;
    end
    
    % Input parameters
    InputParams     = ProcessInputParameters_together_with_cost(vars, N) ;
    Params          = InputParams{1}                  ;
    IsothermParams  = InputParams{2}                  ;
    Times           = InputParams{3}                  ;
    EconomicParams  = InputParams{4}                  ;
    
    % Retrieve process parameters
    N				= Params(1)	      ;
    ro_s			= Params(4)	      ;
    T_0				= Params(5)	      ;
    epsilon			= Params(6)	      ;
    r_p				= Params(7)	      ;
    mu				= Params(8)	      ;
    R				= Params(9)	      ;
    v_0				= Params(10)	  ;
    q_s0            = Params(11)/ro_s ;
    P_0				= Params(17)	  ;
    L				= Params(18)	  ;
    D               = L/3             ;
    MW_CO 			= Params(19)	  ;
    MW_CO2			= Params(20)	  ;
    y_0				= Params(23)	  ;
    ndot_0          = vars(3)	      ;
    P_l				= Params(25)	  ;
    P_inlet			= Params(26)	  ;
    alpha			= Params(30)	  ;
    beta 			= Params(31)	  ;
    y_HR            = Params(33)      ;
    T_HR            = Params(34)      ;
    ndot_HR         = Params(35)	  ;
    MaxTime_CoCPres = 300             ;

    % Call PSA cycle functions
    CoCPressurization_fxn   = @(t, x) FuncCoCPressurization(t, x, Params, IsothermParams)   ;
    Adsorption_fxn          = @(t, x) FuncAdsorption(t, x, Params, IsothermParams)          ;
    LightReflux_fxn         = @(t, x) FuncLightReflux(t, x, Params, IsothermParams)         ;
    HeavyReflux_fxn         = @(t, x) FuncHeavyReflux(t, x, Params, IsothermParams)         ;
    CnCDepressurization_fxn = @(t, x) FuncCnCDepressurization(t, x, Params, IsothermParams) ;
    
    % Call event function
    CoCPres_event_fxn = @(t, x) Eventfunction_CoCPres_step(t, x);

    % Retrieve times of PSA steps
    t_CoCPres   = Times(1) ;
    t_ads       = Times(2) ;
    t_HR        = Times(6) ;
    t_CnCDepres = Times(3) ;
    t_LR        = Times(4) ;
    
    % Dimensionless times 
    tau_CoCPres   = t_CoCPres*v_0/L   ;
    tau_ads       = t_ads*v_0/L       ;
    tau_HR        = t_HR*v_0/L        ;
    tau_CnCDepres = t_CnCDepres*v_0/L ;
    tau_LR        = t_LR*v_0/L        ;
    
    % Initialize the column (initial conditions)
    if nargin < 2 || isempty(x0)
        q                = Isotherm(y_0, P_l, 298.15, IsothermParams) ;
        x0               = zeros(5*N+10,1) ;
        x0(1:N+2)        = P_l/P_0         ;
        x0(N+3)          = y_0             ;
        x0(N+4:2*N+4)    = y_0             ;
        x0(2*N+5:3*N+6)  = q(1)/q_s0       ;
        x0(3*N+7:4*N+8)  = q(2)/q_s0       ;
        x0(4*N+9)        = 1               ;
        x0(4*N+10:5*N+10)= 298.15/T_0      ;
    end 
    
    opts1 = odeset( 'JPattern', JacPressurization(N), 'RelTol', 1e-6, 'Events', CoCPres_event_fxn)       ;
    opts2 = odeset( 'JPattern', JacAdsorption(N), 'RelTol', 1e-6)           ;
    opts3 = odeset( 'JPattern', JacAdsorption(N), 'RelTol', 1e-6)           ;
    opts4 = odeset( 'JPattern', Jac_CnCDepressurization(N), 'RelTol', 1e-6) ;
    opts5 = odeset( 'JPattern', Jac_LightReflux(N), 'RelTol', 1e-6)         ;
%   
%% Begin simulating PSA cycle. This is skipped if the first constraint is violated.
%   Run the simulation until the change in the temperature, gas mole fraction
%   and CO2 molar loading is less than 0.5%. For the molar loading and the
%   mole fraction, the simulation also stops if the absolute change in the
%   state variable is less than 5e-5 and 2.5e-4 respectively. 

if constraints(1) == 0
    
    for i=1:700
        
    %% 1. Simulate CoCPressurization step
        odesol_starttime = tic;
        [t1, a] = ode15s(CoCPressurization_fxn, [0 tau_CoCPres], x0, opts1) ;
        
        % Correct the output (clean up results from simulation)
        idx             = find(a(:, 1) < a(:, 2))         ;  % P_1  < P_2
        a(idx ,1)       = a(idx, 2)                       ;  % P_1  = P_2
        a(idx, N+3)     = a(idx, N+4)                     ;  % y_1  = y_2
        a(idx, 4*N+9)   = a(idx, 4*N+10)                  ;  % T_1  = T_2
        a(:, 2*N+5)     = a(:, 2*N+6)                     ;  % x1_1 = x1_2
        a(:, 3*N+7)     = a(:, 3*N+8)                     ;  % x2_1 = x2_2
        a(:, 3*N+6)     = a(:, 3*N+5)                     ;  % x1_N+2 = x1_N+1
        a(:, 4*N+8)     = a(:, 4*N+7)                     ;  % x2_N+2 = x2_N+1
        a(:, N+3:2*N+4) = max(min(a(:, N+3:2*N+4), 1), 0) ;  % 0 <= y => 1
        
        % Prepare initial conditions for Adsorption step
        x10         = a(end, :)' ;  % Final state of previous step is the
                                    % initial state for current step
        x10(1)      = P_inlet    ;  % BC z=0 P: P_1   = P_inlet
        x10(N+2)    = 1          ;  % BC z=1 P: P_N+2 = 1
        x10(N+3)    = y_0        ;  % BC z=0 y: y_1   = y_0
        x10(2*N+4)  = x10(2*N+3) ;  % BC z=1 y: y_N+2 = y_N+1
        x10(4*N+9)  = 1          ;  % BC z=0 T: T_1   = 1
        x10(5*N+10) = x10(5*N+9) ;  % BC z=1 T: T_N+2 = T_N+1
        
        % Initial conditions of states at first step of the PSA cycle
        statesIC = a(1, [2:N+1, N+4:2*N+3, 2*N+6:3*N+5, 3*N+8:4*N+7, 4*N+10:5*N+9]) ;
%       
    %% 2. Simulate Adsorption step
        [t2, b] = ode15s(Adsorption_fxn, [0 tau_ads], x10, opts2) ;
        
        % Correct the output (clean up results from simulation)
        idx             = find(b(:, N+1) < 1)             ;  % P_N+1 < 1 = P_N+2
        b(idx, N+2)     = b(idx, N+1)                     ;  % P_N+2 = P_N+1
        b(:, 2*N+5)     = b(:, 2*N+6)                     ;  % x1_1 = x1_2
        b(:, 3*N+7)     = b(:, 3*N+8)                     ;  % x2_1 = x2_2
        b(:, 3*N+6)     = b(:, 3*N+5)                     ;  % x1_N+2 = x1_N+1
        b(:, 4*N+8)     = b(:, 4*N+7)                     ;  % x2_N+2 = x2_N+1
        b(:, N+3:2*N+4) = max(min(b(:, N+3:2*N+4), 1), 0) ;  % 0 <= y => 1
        
        if Params(end) == 0
            %b = VelocityCorrection(b, ndot_0, 'HPEnd') ;
            b = velocitycleanup(b)                     ;
        end
        
        % Store final conditions for Adsorption step - all iterations
        % and the CO2 and total moles at the Front and End of the column
        [totalEnd, COEnd, TEnd]  = StreamCompositionCalculator(t2*L/v_0, b, 'LPEnd') ;
        
        % Add and update necessary parameters for Light Reflux step. These
        % are the composition and temperature going out of the adsorption 
        % (light product end), which are those for the inlet of the LR step
        y_LR        = COEnd/totalEnd ;
        T_LR        = TEnd            ;
        ndot_LR     = totalEnd/t_ads  ;
        Params(27)  = y_LR            ;
        Params(28)  = T_LR            ;
        Params(29)  = ndot_LR         ;
        
        % Call the function for Light Reflux step with updated parameters
        LightReflux_fxn = @(t, x) FuncLightReflux(t, x, Params, IsothermParams) ;
        
        % Prepare initial conditions for Heavy Reflux step
        x20         = b(end, :)' ;  % Final state of previous step is the
                                    % initial state for current step
        x20(1)      = P_inlet    ;  % BC z=0 P: P_1   = P_inlet
        x20(1)      = x20(2)     ;   
        x20(N+2)    = 1          ;  % BC z=1 P: P_N+2 = 1
        x20(N+3)    = y_HR       ;  % BC z=0 y: y_1   = y_HR
        x20(2*N+4)  = x20(2*N+3) ;  % BC z=1 y: y_N+2 = y_N+1
        x20(4*N+9)  = T_HR/T_0   ;  % BC z=0 T: T_1   = T_HR/T_0
        x20(5*N+10) = x20(5*N+9) ;  % BC z=1 T: T_N+2 = T_N+1
%       
    %% 3. Simulate Heavy Reflux step
        [t3, c] = ode15s(HeavyReflux_fxn, [0 tau_HR], x20, opts3) ;
        
        % Correct the output (clean up results from simulation)
        idx             = find(c(:, N+1) < 1)             ;  % P_N+1 < 1 = P_N+2
        c(idx, N+2)     = c(idx, N+1)                     ;  % P_N+2 = P_N+1
        c(:, 2*N+5)     = c(:, 2*N+6)                     ;  % x1_1 = x1_2
        c(:, 3*N+7)     = c(:, 3*N+8)                     ;  % x2_1 = x2_2
        c(:, 3*N+6)     = c(:, 3*N+5)                     ;  % x1_N+2 = x1_N+1
        c(:, 4*N+8)     = c(:, 4*N+7)                     ;  % x2_N+2 = x2_N+1
        c(:, N+3:2*N+4) = max(min(c(:, N+3:2*N+4), 1), 0) ;  % 0 <= y => 1
        
        if Params(end) == 0
            c = VelocityCorrection(c, ndot_HR, 'HPEnd') ;
            %c = velocitycleanup(c)                      ;
        end
        
        % Prepare initial conditions for CoCDepressurization step
        x30         = c(end,:)'   ;   % Final state of previous step is the
                                      % initial state for current step
        x30(1)      = x30(2)      ;   % BC z=0 P: P_1   = P_2
        x30(N+2)    = x30(N+1)    ;   % BC z=1 P: P_N+2 = P_N+1
        x30(N+3)    = x30(N+4)    ;   % BC z=0 y: y_1   = y_2
        x30(2*N+4)  = x30(2*N+3)  ;   % BC z=1 y: y_N+2 = y_N+1
        x30(4*N+9)  = x30(4*N+10) ;   % BC z=0 T: T_1   = T_2
        x30(5*N+10) = x30(5*N+9)  ;   % BC z=1 T: T_N+2 = T_N+1
%       
    %% 4. Simulate CnCDepressurization step
        [t4, d] = ode15s(CnCDepressurization_fxn, [0 tau_CnCDepres], x30, opts4) ;
        
        % correct the output (clean up results from simulation)
        idx             = find(d(:, 2) < d(:, 1))         ;  % P_2  < P_1
        d(idx ,1)       = d(idx, 2)                       ;  % P_1  = P_2
        d(:, 2*N+5)     = d(:, 2*N+6)                     ;  % x1_1 = x1_2
        d(:, 3*N+7)     = d(:, 3*N+8)                     ;  % x2_1 = x2_2
        d(:, 3*N+6)     = d(:, 3*N+5)                     ;  % x1_N+2 = x1_N+1
        d(:, 4*N+8)     = d(:, 4*N+7)                     ;  % x2_N+2 = x2_N+1
        d(:, N+3:2*N+4) = max(min(d(:, N+3:2*N+4), 1), 0) ;  % 0 <= y => 1
        
        % Prepare initial donditions for Light Reflux step
        x40         = d(end,:)'    ;  % Final state of previous step is the
                                      % initial state for durrent step
        x40(1)      = P_l/P_0      ;  % Bd z=0 P: P_1   = P_l/P_0
        %x40(N+2)    = 2*P_l/P_0    ;  % Bd z=1 P: P_N+2 = 2*P_l/P_0 % NOTE: be aware of this alpha here, on the other dode is just 2
        x40(N+3)    = x40(N+4)     ;  % Bd z=0 y: y_1   = y_2
        x40(2*N+4)  = y_LR         ;  % Bd z=1 y: y_N+2 = y_LR
        x40(4*N+9)  = x40(4*N+10)  ;  % Bd z=0 T: T_1   = T_2
        x40(5*N+10) = T_LR/T_0     ;  % Bd z=1 T: T_N+2 = T_LR/T_0
%       
    %% 5. Simulate Light Reflux step
        [t5, e] = ode15s(LightReflux_fxn, [0 tau_LR], x40, opts5) ;
        
        % Correct the output (clean up results from simulation)
        idx             = find(e(:, 2) < e(:, 1))         ;  % P_2  < P_1
        e(idx ,1)       = e(idx, 2)                       ;  % P_1  = P_2
        e(:, 2*N+5)     = e(:, 2*N+6)                     ;  % x1_1 = x1_2
        e(:, 3*N+7)     = e(:, 3*N+8)                     ;  % x2_1 = x2_2
        e(:, 3*N+6)     = e(:, 3*N+5)                     ;  % x1_N+2 = x1_N+1
        e(:, 4*N+8)     = e(:, 4*N+7)                     ;  % x2_N+2 = x2_N+1
        e(:, N+3:2*N+4) = max(min(e(:, N+3:2*N+4), 1), 0) ;  % 0 <= y => 1
        
        e = VelocityCorrection(e, ndot_LR*alpha, 'LPEnd') ;
        %e = velocitycleanup(e)                            ;
        
        % Store final conditions for Light Reflux step - all iterations
        % and the CO2 and total moles at the Front and End of the column
        [totalFront, COFront, TFront]  = StreamCompositionCalculator(t5*L/v_0, e, 'HPEnd') ;
        
        % Calculate necessary parameters for Heavy Reflux step
        y_HR       = COFront/totalFront    ;
        T_HR       = TFront                 ;
        ndot_HR    = totalFront.*beta/t_HR  ;
        Params(33) = y_HR                   ;
        Params(34) = T_HR                   ;
        Params(35) = ndot_HR                ;
        
        HeavyReflux_fxn = @(t, x) FuncHeavyReflux(t, x, Params, IsothermParams) ;
        
        % Prepare initial conditions for CoCPressurization step
        x0         = e(end, :)' ;  % Final state of previous step is the
                                   % initial state for current step
        x0(1)      = x0(2)      ;  % BC z=0 P: P_1   = P_2
        x0(N+2)    = x0(N+1)    ;  % BC z=1 P: P_N+2 = P_N+1
        x0(N+3)    = y_0        ;  % BC z=0 y: y_1   = y_0
        x0(2*N+4)  = x0(2*N+3)  ;  % BC z=1 y: y_N+2 = y_N+1
        x0(4*N+9)  = 1          ;  % BC z=0 T: T_1   = 1
        x0(5*N+10) = x0(5*N+9)  ;  % BC z=1 T: T_N+2 = T_N+1
        
        % Final conditions of states at lat step of the PSA cycle
        statesFC = e(end, [2:N+1, N+4:2*N+3, 2*N+6:3*N+5, 3*N+8:4*N+7, 4*N+10:5*N+9]) ;
%       
    %% Check CCS condition
        % CCS condition of states
        CSS_states = norm(statesIC-statesFC) ;
        % Mass balance condition
        [~, ~, massBalance] = ProcessEvaluation(a, b, c, d, e, t1, t2, t3, t4, t5) ;
        
        % Check if CCS has been acheived or not
        if CSS_states <= 1e-3 && abs(massBalance-1) <= 0.005
            break
        end      
    end
    
    %% Process and Economic evaluation
    
    [purity, recovery, MB] = ProcessEvaluation(a, b, c, d, e, t1, t2, t3, t4, t5) ;
    
    %% Column Scheduling
    % by Subraveti et al. 2021
    cycle_time = t_CoCPres + t_ads + t_HR + t_CnCDepres + t_LR       ; % s
    feed_time = t_CoCPres + t_ads ; % s
    
    % minimum number of columns per train
    N_col = ceil(cycle_time/feed_time) ;
    
    % duration of an idle time (included after evacuation (CnCDepres) step 
    % so that the bed profiles are least affected

    %idle_time = N_col*feed_time - cycle_time ; % s
    
    % Number of compressor required for feeding and heavy reflux = 1
    N_comp_feed = 1 ;
    N_comp_hr = 1 ;
       
    % minimum number of evacuation (CnCDepres) and light reflux vacuum pumps
    N_pump_evac = ceil(t_CnCDepres/feed_time) ; 
    N_pump_lr = ceil(t_LR/feed_time) ;
    
    %% Process Scale-up
    % by Subraveti et al. 2021
    total_CO_flow = EconomicParams(1) ; % mol/s
    
    % Calculate the amount of flue gas that is fed during the cycle
    [~, n_CO_pres, ~] = StreamCompositionCalculator(t1*L/v_0, a, 'HPEnd') ;
    [~, n_CO_ads, ~]  = StreamCompositionCalculator(t2*L/v_0, b, 'HPEnd') ;
    CO_fed = n_CO_pres + n_CO_ads ;  % mols/m^2
    
    % Calculate column area
    r_in = D/2                                      ; % m
    Col_Area = pi()*r_in.^2                         ; % m^2

    % Average molar flow rate of the feed to each train
    feed_flow_rate_avg = CO_fed*Col_Area/feed_time ; % mol/s

    % Number of parallel trains
    M_train = ceil(total_CO_flow/feed_flow_rate_avg) ; 
    
    %% Cost Calculation
    % by Turton (2018)
    CEPCI = EconomicParams(3)  ;  % 2024

    %% CAPEX
    % Calculate column cost ($)
    V = Col_Area*L ; % volume of adsorption column (m^3)
    P_H = P_0/1e5 ; % adsorption pressure (bar)
    Cbmc_col = ColumnCost(V, P_H, CEPCI) ;

    % Calculate compressor cost for feeding ($)
    % Asuume feed pressure = 1 bar
    Max_pres_comp_power = MaxCompressorPower(a, 1e5) ;
    Max_ads_comp_power  = MaxCompressorPower(b, 1e5) ;

    Max_comp_power_feed = max(Max_pres_comp_power,Max_ads_comp_power) ; % kW

    [Cbmc_comp_feed, Cbmc_drive_feed] = CompressorCost(Max_comp_power_feed, CEPCI) ;
    
    % Calculate compressor cost for heavy reflux ($)
    Max_hr_comp_power = MaxCompressorPower(c, 1e5) ;
    [Cbmc_comp_hr, Cbmc_drive_hr] = CompressorCost(Max_hr_comp_power, CEPCI) ;

    % Calculate vacuump pump cost for evacuation (CnCDepres) step ($)
    Max_evac_pump_power = MaxVacuumPumPPower(d, 1e5, 'HPEnd') ; % kW
    Cbmc_evac_pump = VacuumPumpCost(Max_evac_pump_power, CEPCI) ;
    
    % Calculate vacuump pump cost for light reflux step ($)
    Max_lr_pump_power = MaxVacuumPumPPower(e, 1e5, 'HPEnd') ; % kW
    Cbmc_lr_pump = VacuumPumpCost(Max_lr_pump_power, CEPCI) ;

    % Calculate initial adsorbent cost ($)
    %adsorbent_price = 2 ; % $/kg
    %initial_adsorbent_cost = adsorbent_price*V*(1-epsilon)*ro_s ;

    % Calculate total bare modular cost (TBMC) ($)
    TBMC = N_col*M_train*Cbmc_col+N_comp_feed*M_train*(Cbmc_comp_feed+Cbmc_drive_feed)+...
             N_comp_hr*M_train*(Cbmc_comp_hr+Cbmc_drive_hr)+...
                N_pump_evac*M_train*Cbmc_evac_pump+N_pump_lr*M_train*Cbmc_lr_pump ;

    % Calculate total modular cost (TMC) ($)
    % Project contingency cost and engineering cost are taken as 15 and 3%
    % of TBMC, respectively.
    TMC = 1.18*TBMC  ;

    % Calultae Grass roots cost or CAPEX ($)
    % A total of 35% of TBMC is assumed as cost of site development,
    % auxiliary buildings, and off-site facilities and utilities.
    CAPEX = TMC+0.50*TBMC    ;

    % Calculate equivalent annual cost (EAC) ($/yr)
    % by Subraveti et al. 2021
    discount_rate = 0.08  ;
    lifetime = 25 ; % yr
    
    EAC = CAPEX*discount_rate/(1-(1/((1+discount_rate)^(lifetime)))) ; % $/yr
    
    %% OPEX 
    % Calculate variable OPEX
    % Calculate cost of electricity ($/yr)
    C_elec = EconomicParams(2)         ; % $/kWh
    
    % Calculate the energy required during the Pressurization step
    % Asuume feed pressure = 1 bar
    E_pres  = CompressionEnergy(t1*L/v_0, a, 1e5)   ; % kWh
    
    % Calculate the energy required during the Feed step
    % Asuume feed pressure = 1 bar
    E_feed  = CompressionEnergy(t2*L/v_0, b, 1e5)   ; % kWh
    
    % Calculate the energy required during the Heavy Reflux Step
    % Asuume feed pressure = 1 bar
    E_HR    = CompressionEnergy(t3*L/v_0, c, 1e5)   ; % kWh
    
    % Calculate the energy required during the Counter Current Depressurization step
    E_evac = VacuumEnergy(t4*L/v_0, d, 1e5, 'HPEnd')        ; % kWh
    
    % Calculate the energy required during the Light Reflux Step
    E_lr  = VacuumEnergy(t5*L/v_0, e, 1e5, 'HPEnd')        ; % kWh
    
    % Calculate the total energy required
    energy_per_cycle = E_pres + E_feed + E_HR + E_evac + E_lr ; % [kWh per cycle]
    
    % Calculate the CO recovered during the cycle [ton CO per cycle and mol/cycle]
    [~, n_CO_CnCD, ~]   = StreamCompositionCalculator(t4*L/v_0, d, 'HPEnd') ;
    [~, n_CO_LR, ~]     = StreamCompositionCalculator(t5*L/v_0, e,  'HPEnd') ;
    CO_recovered_cycle  = (n_CO_CnCD+(1-beta)*n_CO_LR)*r_in^2*pi()*MW_CO/1e3 ;        
    
    %Calculate the energy requirments (kWh/tonne)
    energy_requirments = energy_per_cycle./CO_recovered_cycle             ;

    % Calculate the cost of electricity ($/yr)

    total_CO_mass_flow = total_CO_flow*28/1000/1000*3600*24*365 ; % tonne/yr

    Cost_of_electricity = energy_requirments*recovery*total_CO_mass_flow*C_elec  ;

    % Calculate adsorbent replacement cost ($/yr)
    % It is assumed that adsorbent is completely replaced every 1.5 yrs

    adsorbent_price = 2 ; % $/kg
    Cost_of_adsorbent_replacing = N_col*M_train*V*(1-epsilon)*ro_s*adsorbent_price/1.5 ;

    % Calculate Fixed OPEX
    
    % Calculate labor cost ($/yr)
    % It is assumed that 10 operators including technicians are required
    % for the capture unit at a labor rate of $34.50/h
    labor_rate = 34.50 ; % $/h
    number_of_operators = 10 ;
    labor_cost = labor_rate*number_of_operators*365*24 ;

    % Calculate supervisory cost ($/yr)
    supervisory_cost = 0.25*labor_cost ;

    % Calculate maintenance cost ($/yr)
    % Assume 5 %/yr
    maintenance_cost = 0.1*TMC ;

    % Calculate operating supplies ($/yr)
    operating_supplies_cost = 0.2*maintenance_cost ;

    % Calculate admin overhead cost ($/yr)
    admin_overhead_cost = 0.15*(labor_cost+supervisory_cost+maintenance_cost) ;

    % Calculate plant overhead cost ($/yr)
    plant_overhead_cost = 0.70*(labor_cost+supervisory_cost+maintenance_cost) ;

    % Calculate operating cost, OPEX ($/yr)
    OPEX = Cost_of_electricity+Cost_of_adsorbent_replacing+...
        labor_cost+supervisory_cost+maintenance_cost+...
        operating_supplies_cost+admin_overhead_cost+plant_overhead_cost ;

    %% Capture_cost ($/yr)
    % Calculate total annual cost (TAC) ($/yr)
    TAC = EAC + OPEX ;

    % Calculate capture cost ($/tonne)
    CO_capture_cost = TAC/(recovery*total_CO_mass_flow)  ;
 
    %% Compile objectives and constraint violations
    %con = recovery/MB - 0.9       ;
    %if con < 0
    %    constraints(2) = abs(con) ;
    %end
    
    switch type
        case 'ProcessEvaluation'
    
            objectives  = [0, 0]    ;

            objectives(1) = (1-purity)^2       ;
            objectives(2) = (1-recovery/MB)^2  ;
            
            con = purity - y_0            ;
            if con < 0
                constraints(2) = abs(con) ;
                %constraints(2) = 0        ;
            end
            
        case 'EconomicEvaluation'

            objectives = CO_capture_cost  ;
            
            con = purity - 0.99          ;
            if con < 0
                constraints(2) = abs(con) ;
            end
    end
%   
end 
%   
%% Complementary Functions
   
    function [n_tot, n_CO, Temp] = StreamCompositionCalculator(time, state_vars, ProductEnd)
    %MoleStreamCalculator: Calculate the composition of streams
    %   Calculate the composition (moles of CO and total moles), and 
    %   temperature of a stream at any end of the column end. 
    %   
    %   Input:
    %   time      : Dimensional time vector supplied by the ODE solver
    %   state_vars: Dimensionless state variable matrix supplied for the
    %               step from the ODE solver
    %   ProductEnd: End of the column where the composition will be 
    %               calculated. OPTIONS: HPEnd and LPEnd. HPEnd stands for
    %               heavy product end, which is at the bottom of the column,
    %               and LPEnd stands for light product end, which is at the
    %               top of the column.
    %   
    %   Output:
    %   n_tot     : Average total number of moles at the desired end of the
    %               requested step 
    %   n_CO      : Average number of moles of CO at the desired end of the
    %               requested step 
    %   Temp      : Average temperature at the desired end of the requested
    %               step  
    %   
    %% Check number of inputs
        % If no value is given for the ProductEnd, this is set up by default
        % as HPEnd, which is the heavy product end for any step.
        if nargin < 3
            ProductEnd = 'HPEnd';
        end
    %   
    %%  
        % Differential section length of the column
        dz = L/N ;
        
        % Dimensionalize all variables at the two ends of the columns
        % Collect pressure, temperature and mole fraction at the column end
        % of interest.
        if strcmpi(ProductEnd, 'HPEnd') == 1
        
            P = state_vars(:, 1:2)*P_0     ;
            y = state_vars(:, N+3)         ;
            T = state_vars(:, 4*N+9)*T_0   ;
        
            % Calculate the density of the gas [kg/m^3]
            ro_g = (y*MW_CO + (1-y)*MW_CO2).*P(:, 1)/R./T ;
        
            % calculate concentrations [mol/m^3]
            C_tot = P(:, 1)/R./T ;
            C_CO = C_tot.*y     ;
        
        elseif strcmpi(ProductEnd, 'LPEnd') == 1
        
            P = state_vars(:, N+1:N+2)*P_0 ;
            y = state_vars(:, 2*N+4)       ;
            T = state_vars(:, 5*N+10)*T_0  ;
        
            % calculate the density of the gas [kg/m^3]
            ro_g = (y*MW_CO + (1-y)*MW_CO2).*P(:, 2)/R./T ;
        
            % calculate concentrations [mol/m^3]
            C_tot = P(:, 2)/R./T ;
            C_CO = C_tot.*y     ;
        
        else
            error('Please specify in which end of the column the composition will be calculated. OPTIONS: HPEnd and LPEnd')
        end
        
        % Calculate the pressure gradient at the edges [Pa/m]
        dPdz = 2*(P(:, 2)-P(:, 1))/dz ;
        
        % calculate superficial velocity using ergun equation [m/s]
        viscous_term =  150*mu*(1-epsilon)^2/4/r_p^2/epsilon^3             ;
        kinetic_term = (1.75*(1-epsilon)/2/r_p/epsilon^3) * ro_g           ;
        v            = -sign(dPdz).*(-viscous_term+(abs(viscous_term^2+...  
                        4*kinetic_term.*abs(dPdz))).^(.5))/2./kinetic_term ;
        
        % calculate molar fluxes [mol/m^2/s]
        ndot_tot = abs(v.*C_tot) ;
        ndot_CO = abs(v.*C_CO) ;
        
        % calculate total moles per column area [mol/m^2]
        n_tot = trapz(time, ndot_tot) ;
        n_CO = trapz(time, ndot_CO) ;
        
        % calculate the average temperature of the gas [K]. Only important
        % if the emissions are being used in another step in the cycle (e.g.
        % light reflux, light product pressurization)
        energy_flux_tot = ndot_tot.*T                  ;
        energy_tot      = trapz(time, energy_flux_tot) ;
        Temp            = energy_tot/n_tot         ;
    %   
    end 
    
    function [purity, recovery, mass_balance] = ProcessEvaluation(varargin)
    %ProcessEvaluation: Calculate the process objectives
    %   Calculate the purity and the recovery of the heavy product and the
    %   overall mass balance of the heavy product for the cycle.
    %   
    %   Input:
    %   a-e         : The dimensionless state variables of each step for the
    %                 cycle
    %   t1-t5       : The time steps for each step of the cycle
    %   
    %   Output:
    %   purity      : Purity of the heavy product
    %   recovery    : Recovery of the heavy product
    %   mass_balance: Mass balance of the heavy product (ensure all heavy
    %                 product that enters, leaves. No accumulation)
    %   
    %%  
        step  = cell(nargin/2,1) ;
        tau   = cell(nargin/2,1) ;
        for st = 1:nargin/2
            step{st} = varargin{st}          ;
            tau{st}  = varargin{nargin/2+st} ;
        end
    %   
    %% Calculate total moles going in and out of column [mols/m^2]
        
        [~, n_CO_CoCPres_HPEnd, ~]                       = StreamCompositionCalculator(tau{1}, step{1}, 'HPEnd') ;
        [~, n_CO_ads_HPEnd, ~]                           = StreamCompositionCalculator(tau{2}, step{2}, 'HPEnd') ;
        [~, n_CO_ads_LPEnd, ~]                           = StreamCompositionCalculator(tau{2}, step{2}, 'LPEnd') ;
        [~, n_CO_HR_LPEnd, ~]                            = StreamCompositionCalculator(tau{3}, step{3}, 'LPEnd') ;
        [~, n_CO_HR_HPEnd, ~]                            = StreamCompositionCalculator(tau{3}, step{3}, 'HPEnd') ;
        [n_tot_CnCDepres_HPEnd, n_CO_CnCDepres_HPEnd, ~] = StreamCompositionCalculator(tau{4}, step{4}, 'HPEnd') ;
        [~, n_CO_LR_LPEnd, ~]                            = StreamCompositionCalculator(tau{5}, step{5}, 'LPEnd') ;
        [n_tot_LR_HPEnd, n_CO_LR_HPEnd, ~]               = StreamCompositionCalculator(tau{5}, step{5}, 'HPEnd') ;
	%   
    %% Calculate Purity, recovery and mass balance of the column
        
        purity       = (n_CO_CnCDepres_HPEnd+(1-beta)*n_CO_LR_HPEnd)/(n_tot_CnCDepres_HPEnd+(1-beta)*n_tot_LR_HPEnd) ;
        recovery     = (n_CO_CnCDepres_HPEnd+(1-beta)*n_CO_LR_HPEnd)/(n_CO_CoCPres_HPEnd+n_CO_ads_HPEnd)           ;
        
        mass_balance = (n_CO_CnCDepres_HPEnd+n_CO_ads_LPEnd+n_CO_HR_LPEnd+n_CO_LR_HPEnd)/... 
                       (n_CO_CoCPres_HPEnd+n_CO_ads_HPEnd+n_CO_HR_HPEnd+n_CO_LR_LPEnd)                             ;
    %   
    end 
    
    function Maxpower = MaxCompressorPower(state_vars, Patm)
    %MaxCompressorPower: Calculate the maximum compression power from entrance
    %   Calculate the maximum power required by the compressor to increase
    %   the pressure to the desired value. 
    %   
    %   Input:
    %   state_vars: Dimensionless state variable matrix supplied for the
    %               step from the ODE solver
    %   Patm      : Atmospheric pressure [Pa]. This value can also be changed
    %               if the flue gas is at a higher pressure. This value is
    %               the cutoff for the compressor. Above this value, energy
    %               is required. Below this value, no energy is required
    %   
    %   Output:
    %   MaxPower  : Maximum compressor power [kW]  
    %   
    %% Check number of inputs
    
        % Differential section length of the column
        dz = L/N  ;
        
        % Compresor parameters
        adiabatic_index       = 1.4  ;
        compressor_efficiency = 0.72 ;
        
        % Calculate the pressure gradient at the edges [Pa/m]
        P = state_vars(:, 1:2)*P_0     ;
        y = state_vars(:, N+3)         ;
        T = state_vars(:, 4*N+9)*T_0   ;
        
        % Calculate the density of the gas [kg/m^3]
        ro_g = (y*MW_CO + (1-y)*MW_CO2).*P(:, 1)/R./T ;
        
        dPdz = 2*(P(:, 2)-P(:, 1))/dz ;  % [Pa/m]
        
        % calculate superficial velocity using ergun equation [m/s]
        viscous_term =  150*mu*(1-epsilon)^2/4/r_p^2/epsilon^3             ;
        kinetic_term = (1.75*(1-epsilon)/2/r_p/epsilon^3) * ro_g           ;
        v            = -sign(dPdz).*(-viscous_term+(abs(viscous_term^2+...  
                        4*kinetic_term.*abs(dPdz))).^(.5))/2./kinetic_term ;
        
        % Calculate the compression ratio along with the impact it has on the
        % energy requirments
        ratio_term = ((P(:, 1)/Patm).^((adiabatic_index-1)/adiabatic_index)-1) ;
        ratio_term = max(ratio_term, 0)                                        ;
        
        %Calculate the total energy required by the compressor
        integral_term = abs(v.*P(:, 1).*ratio_term) ;
        
        power = integral_term.*((adiabatic_index)./(adiabatic_index-1))./compressor_efficiency*pi()*r_in.^2 ;  % W
        
        Maxpower = max(power)/1000 ; % W to kW
    end
    
    function Maxpower = MaxVacuumPumPPower(state_vars, Patm, ProductEnd)
    %MaxVacuumPumPPower: Calculate the maximum power at both ends of the column
    %   Calculate the maximum power required by the vacuum pump to operate 
    %   at the desired pressure. 
    %   
    %   Input:
    %   state_vars: Dimensionless state variable matrix supplied for the
    %               step from the ODE solver
    %   ProductEnd: End of the column where the composition will be 
    %               calculated. OPTIONS: HPEnd and LPEnd. HPEnd stands for
    %               heavy product end, which is at the bottom of the column,
    %               and LPEnd stands for light product end, which is at the
    %               top of the column.
    %   Patm      : Atmospheric pressure [Pa].This value is the cutoff for 
    %               the vacuum. Above this value, energy is not required.
    %               Below this value, energy is required
    %   
    %   Output:
    %   Maxpower    : Maximum vacuum pump power [kW]  
    %   
    %%  
        % Differential section length of the column
        dz = L/N ;
        
        % Vacuum parameters
        adiabatic_index   = 1.4  ;
        vacuum_efficiency = 0.72 ;
        
        % Dimensionalize all variables at the two ends of the columns
        % Collect pressure, temperature and mole fraction at the column end
        % of interest.
        if strcmpi(ProductEnd, 'HPEnd') == 1
        
            P = state_vars(:, 1:2)*P_0     ;
            y = state_vars(:, N+3)         ;
            T = state_vars(:, 4*N+9)*T_0   ;
        
            % Calculate the density of the gas [kg/m^3]
            ro_g = (y*MW_CO + (1-y)*MW_CO2).*P(:, 1)/R./T ;
            
            P_out = P(:, 1) ;
        
        elseif strcmpi(ProductEnd, 'LPEnd') == 1
        
            P = state_vars(:, N+1:N+2)*P_0 ;
            y = state_vars(:, 2*N+4)       ;
            T = state_vars(:, 5*N+10)*T_0  ;
        
            % calculate the density of the gas [kg/m^3]
            ro_g = (y*MW_CO + (1-y)*MW_CO2).*P(:, 2)/R./T ;
            
            P_out = P(:, 2) ;
        
        else
            error('Please specify in which end of the column the composition will be calculated. OPTIONS: HPEnd and LPEnd')
        end
        
        % Calculate the pressure gradient at the edges [Pa/m]
        dPdz = 2*(P(:, 2)-P(:, 1))/dz ;
        
        % calculate superficial velocity using ergun equation [m/s]
        viscous_term =  150*mu*(1-epsilon)^2/4/r_p^2/epsilon^3             ;
        kinetic_term = (1.75*(1-epsilon)/2/r_p/epsilon^3) * ro_g           ;
        v            = -sign(dPdz).*(-viscous_term+(abs(viscous_term^2+...  
                        4*kinetic_term.*abs(dPdz))).^(.5))/2./kinetic_term ;
        
        % Calculate the compression ratio along with the impact it has on the
        % energy requirments
        ratio_term = ((Patm./P_out).^((adiabatic_index-1)/adiabatic_index)-1) ;
        ratio_term = max(ratio_term, 0)                                       ;
        
        integral_term = abs(v.*P_out.*ratio_term) ;
        
        %Calculate the total energy required by the vacuum pump  
        power=integral_term.*((adiabatic_index)./(adiabatic_index-1))./vacuum_efficiency*pi()*r_in.^2;  % W
        
        Maxpower = max(power)/1000 ; % W to kW
    %   
    end 

    function energy = CompressionEnergy(time, state_vars, Patm)
    %CompressionEnergy: Calculate the  compression energy from entrance
    %   Calculate the total energy required by the compressor to increase
    %   the pressure to the desired value. Intended to be used with PSA
    %   dimensionless simulations. 
    %   
    %   Input:
    %   time      : Dimensional time vector supplied by the ODE solver
    %   state_vars: Dimensionless state variable matrix supplied for the
    %               step from the ODE solver
    %   Patm      : Atmospheric pressure [Pa]. This value can also be changed
    %               if the flue gas is at a higher pressure. This value is
    %               the cutoff for the compressor. Above this value, energy
    %               is required. Below his value, no energy is required
    %   
    %   Output:
    %   energy    : Energy Requirments [kWh]  
    %   
    %% Check number of inputs
    
        % Differential section length of the column
        dz = L/N  ;
        
        % Compresor parameters
        adiabatic_index       = 1.4  ;
        compressor_efficiency = 0.72 ;
        
        % Calculate the pressure gradient at the edges [Pa/m]
        P = state_vars(:, 1:2)*P_0     ;
        y = state_vars(:, N+3)         ;
        T = state_vars(:, 4*N+9)*T_0   ;
        
        % Calculate the density of the gas [kg/m^3]
        ro_g = (y*MW_CO + (1-y)*MW_CO2).*P(:, 1)/R./T ;
        
        dPdz = 2*(P(:, 2)-P(:, 1))/dz ;
        
        % calculate superficial velocity using ergun equation [m/s]
        viscous_term =  150*mu*(1-epsilon)^2/4/r_p^2/epsilon^3             ;
        kinetic_term = (1.75*(1-epsilon)/2/r_p/epsilon^3) * ro_g           ;
        v            = -sign(dPdz).*(-viscous_term+(abs(viscous_term^2+...  
                        4*kinetic_term.*abs(dPdz))).^(.5))/2./kinetic_term ;
        
        % Calculate the compression ratio along with the impact it has on the
        % energy requirments
        ratio_term = ((P(:, 1)/Patm).^((adiabatic_index-1)/adiabatic_index)-1) ;
        ratio_term = max(ratio_term, 0)                                        ;
        
        %Calculate the total energy required by the compressor
        integral_term = abs(v.*P(:, 1).*ratio_term) ;
        
        energy = trapz(time, integral_term).*((adiabatic_index)./(adiabatic_index-1))./compressor_efficiency*pi()*r_in.^2 ;
        
        energy = energy/3.6e6 ;
    end
    
    function energy = VacuumEnergy(time, state_vars, Patm, ProductEnd)
    %VacuumEnergy: Calculate the vacuum energy at both ends of the column
    %   Calculate the total energy required by the vacuum pump to operate 
    %   at the desired pressure. Intended to be used with PSA dimensionless 
    %   simulations. 
    %   
    %   Input:
    %   time      : Dimensional time vector supplied by the ODE solver
    %   state_vars: Dimensionless state variable matrix supplied for the
    %               step from the ODE solver
    %   ProductEnd: End of the column where the composition will be 
    %               calculated. OPTIONS: HPEnd and LPEnd. HPEnd stands for
    %               heavy product end, which is at the bottom of the column,
    %               and LPEnd stands for light product end, which is at the
    %               top of the column.
    %   Patm      : Atmospheric pressure [Pa].This value is the cutoff for 
    %               the vacuum. Above this value, energy is not required.
    %               Below his value, energy is required
    %   
    %   Output:
    %   energy    : Energy Requirments [kWh]
    %   
    %% Check number of inputs
        % If no value is given for the ProductEnd, this is set up by default
        % as HPEnd, which is the heavy product end for any step.
        if nargin < 4
            ProductEnd = 'HPEnd';
        end
    %   
    %%  
        % Differential section length of the column
        dz = L/N ;
        
        % Vacuum parameters
        adiabatic_index   = 1.4  ;
        vacuum_efficiency = 0.72 ;
        
        % Dimensionalize all variables at the two ends of the columns
        % Collect pressure, temperature and mole fraction at the column end
        % of interest.
        if strcmpi(ProductEnd, 'HPEnd') == 1
        
            P = state_vars(:, 1:2)*P_0     ;
            y = state_vars(:, N+3)         ;
            T = state_vars(:, 4*N+9)*T_0   ;
        
            % Calculate the density of the gas [kg/m^3]
            ro_g = (y*MW_CO + (1-y)*MW_CO2).*P(:, 1)/R./T ;
            
            P_out = P(:, 1) ;
        
        elseif strcmpi(ProductEnd, 'LPEnd') == 1
        
            P = state_vars(:, N+1:N+2)*P_0 ;
            y = state_vars(:, 2*N+4)       ;
            T = state_vars(:, 5*N+10)*T_0  ;
        
            % calculate the density of the gas [kg/m^3]
            ro_g = (y*MW_CO + (1-y)*MW_CO2).*P(:, 2)/R./T ;
            
            P_out = P(:, 2) ;
        
        else
            error('Please specify in which end of the column the composition will be calculated. OPTIONS: HPEnd and LPEnd')
        end
        
        % Calculate the pressure gradient at the edges [Pa/m]
        dPdz = 2*(P(:, 2)-P(:, 1))/dz ;
        
        % calculate superficial velocity using ergun equation [m/s]
        viscous_term =  150*mu*(1-epsilon)^2/4/r_p^2/epsilon^3             ;
        kinetic_term = (1.75*(1-epsilon)/2/r_p/epsilon^3) * ro_g           ;
        v            = -sign(dPdz).*(-viscous_term+(abs(viscous_term^2+...  
                        4*kinetic_term.*abs(dPdz))).^(.5))/2./kinetic_term ;
        
        % Calculate the compression ratio along with the impact it has on the
        % energy requirments
        ratio_term = ((Patm./P_out).^((adiabatic_index-1)/adiabatic_index)-1) ;
        ratio_term = max(ratio_term, 0)                                       ;
        
        integral_term = abs(v.*P_out.*ratio_term) ;
        
        %Calculate the total energy required by the compressor  
        energy=trapz(time, integral_term).*((adiabatic_index)./(adiabatic_index-1))./vacuum_efficiency*pi()*r_in.^2;
        
        energy = energy/3.6e6 ;
    %   
    end 
    
    function x_new = VelocityCorrection(x, n_hr, CorrectionEnd)
    %% Check number of inputs
        % If no value is given for the CorrectionEnd, this is set up by
        % default as HPEnd, which is the heavy product end for any step.
        if nargin < 3
            CorrectionEnd = 'HPEnd';
        end
    %   
    %%  
        x_new = x   ;
        
        % Differential section length of the column
        dz    = L/N ;
        
        % Dimensionalize all variables at the two ends of the columns
        if strcmpi(CorrectionEnd, 'HPEnd') == 1
        
            T = x(:, 4*N+9)*T_0 ;
            y = x(:, N+3)       ;
            P = x(:, 2)*P_0     ;
        
        elseif strcmpi(CorrectionEnd, 'LPEnd') == 1
        
            T = x(:, 5*N+10)*T_0 ;
            y = x(:, 2*N+4)      ;
            P = x(:, N+1)*P_0    ;
        
        else
            error('Please specify in which end of the column the velocity correction will be calculated. OPTIONS: HPEnd and LPEnd')
        end
        
        MW = MW_CO2+(MW_CO-MW_CO2)*y ;
        
        a_1   = 150*mu*(1-epsilon)^2*dz/2/4/r_p^2/epsilon^3/R./T    ;
        a_2_1 = 1.75*(1-epsilon)/2/r_p/epsilon/epsilon/epsilon*dz/2 ;
        a_2   = a_2_1/R./T*n_hr.*MW                                 ;
        
        a_a =  a_1+a_2  ;
        b_b =  P./T/R   ;
        c_c = -n_hr     ;
        
        vh = (-b_b+sqrt(b_b.^2-4.*a_a.*c_c))/2./a_a ;
        
        a_p = a_1.*T*R       ;
        b_p = a_2_1.*MW/R./T ;
        
        % Correction
        if strcmpi(CorrectionEnd, 'HPEnd') == 1
        
            x_new(:, 1)   = ((a_p.*vh+P)./(1-b_p.*vh.*vh))./P_0 ;
        
        elseif strcmpi(CorrectionEnd, 'LPEnd') == 1
        
            x_new(:, N+2) = ((a_p.*vh+P)./(1-b_p.*vh.*vh))./P_0 ;
        
        end
	%   
    end 
    
    function x_new = velocitycleanup(x)

        x_new=x;
        numb1=150*mu*(1-epsilon)^2/4/r_p^2/epsilon^2;
        ro_gent=x(:,2).*P_0/R/T_0;
        numb2_ent=(ro_gent.*(MW_CO2+(MW_CO-MW_CO2)*x(:,N+3)).*(1.75*(1-epsilon)/2/r_p/epsilon));


        x_new(:,1)=(numb1*v_0+numb2_ent*v_0*v_0)*L/P_0/2/N + x(:,2);
    end 

    function C_col = ColumnCost(V, p_H, CEPCI)
    %Column_Cost: Calculate the bare modular costs for adsorption column
    %   Calculate the bare modular costs (BMC) for adsorption column by Turton
    %   (2018)
    %   
    %   Input:
    %   V         : Volume of an adsorption column (m^3)
    %   p_H       : adsorption pressure (bar)
    %   CEPCI     : CEPCI of present year (2024).
    %   
    %   Output:
    %   C_col     : Bare modular costs for adsorption column [$]
    %%  
        % Calculate column purchase cost for a base condition 
        Input_col = log10(V) ;
        Exp_term_col = 3.4974+0.4485*Input_col +0.1074*Input_col*Input_col  ;
        Cb_col = 10^(Exp_term_col) ; % $

        % Calculate pressur factor
        Fp_left_numerator = p_H*2*r_in ;
        Fp_left_denominator = 2*(850-0.6*p_H) ;
        
        Fp_left = ((Fp_left_numerator/Fp_left_denominator)+0.00315)/0.0063 ;
        
        Fp = max(Fp_left, 1.25) ;
        
        % Calculate cost factor to account for indirect cost (freight,
        % overhead, and engineering cost) plus materials of construction
        % and operating pressure.
        
        % Calculate material factor
        FM = 1 ; % for carbon steel

        % Calculate bare modele factor (cost factor)
        FBM = 2.25+1.82*FM*Fp ;
        
        % CEPCI (base)
        CEPCI_base = 397 ; % 2001
        
        % Calculate BMC for column
        C_col = Cb_col*CEPCI/CEPCI_base*FBM ;
    end 

    function [C_comp, C_drive] = CompressorCost(W, CEPCI)
    %CompressorCost: Calculate the bare modular costs for compressor
    %   Calculate the bare modular costs (BMC) for compressor by Turton
    %   (2018)
    %   
    %   Input:
    %   W         : Maximum shaft power of compressor (kW)
    %   CEPCI     : CEPCI of present year (2024).
    %   
    %   Output:
    %   C_comp     : Bare modular costs for compressor [$]
    %   C_drive    : Bare modular costs for drive [$]
    %%  
        % Calculate column purchase cost for a base condition 
        compressor_efficiency = 0.72 ;
        Input_comp = log10(W*compressor_efficiency) ; % fluid power (kW)
        Input_drive = log10(W) ;  % compressor shaft work (kW)
        
        Exp_term_comp = 2.2897+1.3604*Input_comp-0.1027*Input_comp*Input_comp ;
        Exp_term_drive = 1.956+1.7142*Input_drive-0.2282*Input_drive*Input_drive ;
        
        Cb_comp = 10^(Exp_term_comp) ;
        Cb_drive = 10^(Exp_term_drive) ;
        
        % CEPCI (base)
        CEPCI_base = 397 ; % 2001
        
        % Calculate BMC for compressor and drive
        % FBM for compressor = 2.7 for carbon steel
        % FBM for drive = 1.5 for carbon steel
        C_comp = 2.7*CEPCI/CEPCI_base*Cb_comp ; % $
        C_drive = 1.5*CEPCI/CEPCI_base*Cb_drive ; % $
    end 

    function C_pump = VacuumPumpCost(W, CEPCI)
    %Column_Cost: Calculate the bare modular costs for vacuum pumps
    %   Calculate the bare modular costs (BMC) for vacuum pumps by Turton
    %   (2018)
    %   
    %   Input:
    %   W         : Maximum shaft power of vacuum pump (kW)
    %   CEPCI     : CEPCI of present year (2024).
    %   
    %   Output:
    %   C_pump     : Bare modular costs for vacuum pumps [$]
    %%  
        % Calculate column purchase cost for a base condition 
        Input_vacuum_pump = log10(W) ;
        Exp_term_vacuum_pump = 3.3892+0.0536*Input_vacuum_pump+0.1536*Input_vacuum_pump*Input_vacuum_pump ;
        Cb_pump = 10^(Exp_term_vacuum_pump) ;

        % Calculate pressur factor
        Fp = 1 ; % for carbon steel
        
        % Calculate material factor
        FM = 1.6 ; % for carbon steel

        % Calculate bare modele factor (cost factor)
        FBM = 1.89+1.35*FM*Fp ;

        % CEPCI (base)
        CEPCI_base = 397 ; % 2001
        
        % Calculate BMC for column
        C_pump = FBM*CEPCI/CEPCI_base*Cb_pump ;
    end 
%   
%% Jacobian patterns functions
    
    function J_pres      = JacPressurization(N)
    %JacPressurization: Calculates a Jacobian pattern for the step
    %   This function calculates a jacobian pattern for the pressurization
    %   step. This is to be used with the ode solver to decrease the 
    %   computational time
    %   
    %   Input:
    %       N     : Number of finite volumes used in the column
    %   Output:
    %       J_pres: The sparse Jacobian pattern scheme
    %   
    %% Create individual segments
        % Four band Jacobian scheme for advection terms
        B4 = ones(N+2, 4)                      ;
        A4 = full(spdiags(B4, -2:1, N+2, N+2)) ;
        
        % One band Jacobian scheme for adsorption/desoprtion term
        B1           = ones(N+2, 1)                   ;
        A1           = full(spdiags(B1, 0, N+2, N+2)) ;
        A1(1, 1)     = 0                              ;
        A1(N+2, N+2) = 0                              ;
        
        % Zero band Jacobian Term
        A0 = zeros(N+2) ;
	%   
    %% Create Overall Jacobian based on individual segments
        J_pres = [ A4, A4, A1, A1, A4;... 
                   A4, A4, A1, A1, A4;... 
                   A1, A1, A1, A0, A1;... 
                   A1, A1, A0, A1, A1;... 
                   A4, A1, A1, A1, A4 ]  ;
	%   
    %% Modify based on boundary conditions
        % Pressure Inlet
        J_pres(1, :) = 0 ;
        J_pres(1, 1) = 1 ;
        
        % Pressure Outlet
        J_pres(N+2, :) = J_pres(N+1, :) ;
        J_pres(:, N+2) = 0              ;
        
        % Mole Fraction Inlet
        J_pres(N+3, :) = 0 ;
        J_pres(:, N+3) = 0 ;
        
        % Mole Fraction Outlet
        J_pres(2*N+4, :) = J_pres(2*N+3, :) ;
        J_pres(:, 2*N+4) = 0                ;
        
        % Temperature Inlet
        J_pres(4*N+9, :) = 0 ;
        J_pres(:, 4*N+9) = 0 ;
        
        % Temperature Outlet
        J_pres(5*N+10, :) = J_pres(5*N+9, :) ;
        J_pres(:, 5*N+10) = 0                ;
        
        J_pres=sparse(J_pres);
	%   
    end 
    
    function J_ads       = JacAdsorption(N)
    %JacAdsorption: Calculates a Jacobian pattern for the step
    %   This function calculates a jacobian pattern for the pressurization
    %   step. This is to be used with the ode solver to decrease the 
    %   computational time
    %   
    %   Input:
    %       N    : Number of finite volumes used in the column
    %   Output:
    %       J_ads: The sparse Jacobian pattern scheme
    %   
    %% Create individual segments
        % Four band Jacobian scheme for advection terms
        B4 = ones(N+2, 4)                      ;
        A4 = full(spdiags(B4, -2:1, N+2, N+2)) ;
        
        % One band Jacobian scheme for adsorption/ desoprtion term
        B1           = ones(N+2, 1)                   ;
        A1           = full(spdiags(B1, 0, N+2, N+2)) ;
        A1(1, 1)     = 0                              ;
        A1(N+2, N+2) = 0                              ;
        
        % Zero band Jacobian Term
        A0 = zeros(N+2) ;
	%   
    %% Create Overall Jacobian based on individual segments
        J_ads = [ A4, A4, A1, A1, A4;... 
                  A4, A4, A1, A1, A4;... 
                  A1, A1, A1, A0, A1;... 
                  A1, A1, A0, A1, A1;... 
                  A4, A1, A1, A1, A4 ]  ;
	%   
    %% Modify based on boundary conditions
        % Pressure Inlet
        J_ads(1, :) = 0 ;
        J_ads(:, 1) = 0 ;
        
        % Pressure Outlet
        J_ads(N+2, :) = 0 ;
        J_ads(:, N+2) = 0 ;
        
        % Mole Fraction Inlet
        J_ads(N+3, :) = 0 ;
        J_ads(:, N+3) = 0 ;
        
        % Mole Fraction Outlet
        J_ads(2*N+4, :) = J_ads(2*N+3, :) ;
        J_ads(:, 2*N+4) = 0               ;
        
        % Temperature Inlet
        J_ads(4*N+9, :) = 0 ;
        J_ads(:, 4*N+9) = 0 ;
        
        % Temperature Outlet
        J_ads(5*N+10, :) = J_ads(5*N+9, :) ;
        J_ads(:, 5*N+10) = 0               ;
        
        J_ads = sparse(J_ads) ;
	%   
    end 
    
    function J_CnCdepres = Jac_CnCDepressurization(N)
    %Jac_CnCDepressurization: Calculates a Jacobian pattern for the step
    %   This function calculates a jacobian pattern for the pressurization
    %   step. This is to be used with the ode solver to decrease the 
    %   computational time
    %   
    %   Input:
    %       N          : Number of finite volumes used in the column
    %   Output:
    %       J_CnCdepres: The sparse Jacobian pattern scheme
    %   
    %% Create individual segments
        % Four band Jacobian scheme for advection terms
        B4         = ones(N+2, 4)                      ;
        A4         = full(spdiags(B4, -1:2, N+2, N+2)) ;
        A4(1, :)   = A4(2, :)                          ;
        A4(N+2, :) = A4(N+1, :)                        ;
        
        % One band Jacobian scheme for adsorption/ desoprtion term
        B1           = ones(N+2, 1)                   ;
        A1           = full(spdiags(B1, 0, N+2, N+2)) ;
        A1(1, 1)     = 0                              ;
        A1(N+2, N+2) = 0                              ;
        A1(1, 2)     = 1                              ;
        A1(N+2, N+1) = 1                              ;
        
       % Zero band Jacobian Term
        A0 = zeros(N+2) ;
	%   
    %% Create Overall Jacobian based on individual segments
        J_CnCdepres = [ A4, A4, A1, A1, A4;... 
                        A4, A4, A1, A1, A4;... 
                        A1, A1, A1, A0, A1;... 
                        A1, A1, A0, A1, A1;... 
                        A4, A1, A1, A1, A4 ]  ;
	%   
    %% Modify based on boundary conditions
        % Pressure Inlet
        J_CnCdepres(1, :) = 0 ;
        J_CnCdepres(1, 1) = 1 ;
        
        % Pressure Outlet
        J_CnCdepres(N+2, :) = J_CnCdepres(N+1, :) ;
        
        % Mole Fraction Inlet
        J_CnCdepres(N+3, :) = J_CnCdepres(N+4, :) ;
        
        % Mole Fraction Outlet
        J_CnCdepres(2*N+4) = J_CnCdepres(2*N+3) ;
        
        % Molar Loading
        J_CnCdepres(2*N+5, :)       = 0 ;
        J_CnCdepres(3*N+6:3*N+7, :) = 0 ;
        J_CnCdepres(4*N+8, :)       = 0 ;
        
        % Temperature Inlet
        J_CnCdepres(4*N+9, :) = J_CnCdepres(4*N+10, :) ;
        
        %Temperature Outlet
        J_CnCdepres(5*N+10, :) = J_CnCdepres(5*N+9) ;
        
        J_CnCdepres = sparse(J_CnCdepres) ;
	%   
    end 
    
    function J_LR        = Jac_LightReflux(N)
    %Jac_LightReflux: Calculates a Jacobian pattern for the step
    %   This function calculates a jacobian pattern for the pressurization
    %   step. This is to be used with the ode solver to decrease the 
    %   computational time
    %   
    %   Input:
    %       N   : Number of finite volumes used in the column
    %   Output:
    %       J_LR: The sparse Jacobian pattern scheme
    %   
    %% Create individual segments
        % Four band Jacobian scheme for advection terms
        B4         = ones(N+2, 4)                      ;
        A4         = full(spdiags(B4, -1:2, N+2, N+2)) ;
        A4(1, :)   = A4(2, :)                          ;
        A4(N+2, :) = A4(N+1, :)                        ;
        
        % One band Jacobian scheme for adsorption/ desoprtion term
        B1           = ones(N+2, 1)                   ;
        A1           = full(spdiags(B1, 0, N+2, N+2)) ;
        A1(1, 1)     = 0                              ;
        A1(N+2, N+2) = 0                              ;
        A1(1, 2)     = 1                              ;
        A1(N+2, N+1) = 1                              ;
        
        % Zero band Jacobian Term
        A0 = zeros(N+2) ;
	%   
    %% Create Overall Jacobian based on individual segments
        J_LR = [ A4, A1, A1, A1, A4;... 
                 A4, A4, A1, A1, A4;... 
                 A1, A1, A1, A0, A1;... 
                 A1, A1, A0, A1, A1;... 
                 A4, A1, A1, A1, A4 ]  ;
	%   
    %% Modify based on boundary conditions
        % Pressure Inlet
        J_LR(1, :) = 0 ;
        J_LR(1, 1) = 1 ;
        
        % Pressure Outlet
        J_LR(N+2, :) = J_LR(N+1, :) ;
        
        % Mole Fraction Inlet
        J_LR(N+3, :) = J_LR(N+4, :) ;
        
        % Mole Fraction Outlet
        J_LR(2*N+4) = J_LR(2*N+3) ;
        
        % Molar Loading
        J_LR(2*N+5, :)       = 0 ;
        J_LR(3*N+6:3*N+7, :) = 0 ;
        J_LR(4*N+8, :)       = 0 ;
        
        % Temperature Inlet
        J_LR(4*N+9, :) = J_LR(4*N+10, :) ;
        
        % Temperature Outlet
        J_LR(5*N+10, :) = J_LR(5*N+9) ;
        
        J_LR = sparse(J_LR) ;
	%   
    end 
%   
%% Pres Termination Function
    
    function [value, isterminal, direction] = Eventfunction_CoCPres_step(t, x)
    %Eventfunction_CoCPres_step: Function used with ode15s to stop the step
    %   Termination function for the CoC Pressurization step to stop the
    %   solution once the maximum time is reached

    %%
        if (toc(odesol_starttime)-MaxTime_CoCPres) >=0
            value = 0;
            error("Maximum time for ode evaluation for CoC Pressurization step reached. Terminatin the solution")
        else
            value = 1;
        end
        isterminal = 1;
        direction = 0;
    end
end 
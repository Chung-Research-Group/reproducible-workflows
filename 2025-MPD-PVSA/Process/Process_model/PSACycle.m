function  [objectives, constraints, a, b, c, d, e, t1, t2, t3, t4, t5] = PSACycle(vars, x0, type, N)
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
    InputParams     = ProcessInputParameters(vars, N) ;
    Params          = InputParams{1}                  ;
    MPDParams       = InputParams{2}                  ;
    Times           = InputParams{3}                  ;
    LN_Pi           = InputParams{4}                  ;
    Ene             = InputParams{5}                  ;
    state_rd        = InputParams{6}                  ;
    %EconomicParams  = InputParams{4}                  ;
    
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
    %D               = L/3             ;
    MW_CO			= Params(19)	  ;
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
    Maxpurgetime    = 3600            ;
    Maxpresstime    = 300             ;

    NumOfCom1       = MPDParams(1)    ;
    NumOfCom2       = MPDParams(2)    ;
    Cell_a          = MPDParams(3)    ;
    Cell_b          = MPDParams(4)    ;
    Cell_c          = MPDParams(5)    ;
    Unit            = MPDParams(6)    ;

    % Call PSA cycle functions
    CoCPressurization_fxn   = @(t, x) FuncCoCPressurization(t, x, Params, MPDParams, LN_Pi, Ene, state_rd)   ;
    Adsorption_fxn          = @(t, x) FuncAdsorption(t, x, Params, MPDParams, LN_Pi, Ene, state_rd)          ;
    LightReflux_fxn         = @(t, x) FuncLightReflux(t, x, Params, MPDParams, LN_Pi, Ene, state_rd)         ;
    HeavyReflux_fxn         = @(t, x) FuncHeavyReflux(t, x, Params, MPDParams, LN_Pi, Ene, state_rd)         ;
    CnCDepressurization_fxn = @(t, x) FuncCnCDepressurization(t, x, Params, MPDParams, LN_Pi, Ene, state_rd) ;

    % Call event function
    myeventfunc_CoCPress_step            = @(t, x) eventfunction_CoCpress_step(t, x);
    myeventfunc_Adsorption_step          = @(t, x) eventfunction_adsorption_step(t, x);
    myeventfunc_heavy_purge_step         = @(t, x) eventfunction_heavy_purge_step(t, x);
    myeventfunc_CnCDepress_step          = @(t, x) eventfunction_CnCdepress_step(t, x);
    myeventfunc_light_purge_step         = @(t, x) eventfunction_light_purge_step(t, x);
    
    % Call event function
    %CoCPres_event_fxn = @(t, x) Eventfunction_CoCPres_step(t, x);

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
        q                = calculate_N_loads(y_0, (1-y_0), P_l/1e5, 298.15, LN_Pi, Ene, state_rd, NumOfCom1, NumOfCom2, Cell_a, Cell_b, Cell_c) ;
        x0               = zeros(5*N+10,1) ;
        x0(1:N+2)        = P_l/P_0         ;
        x0(N+3)          = y_0             ;
        x0(N+4:2*N+4)    = y_0             ;
        x0(2*N+5:3*N+6)  = q(1)*Unit/q_s0  ;
        x0(3*N+7:4*N+8)  = q(2)*Unit/q_s0  ;
        x0(4*N+9)        = 1               ;
        x0(4*N+10:5*N+10)= 298.15/T_0      ;
    end 
    
    opts1 = odeset('Events', myeventfunc_CoCPress_step, 'RelTol',1e-06, 'AbsTol',1e-09) ;
    opts2 = odeset('Events', myeventfunc_Adsorption_step, 'RelTol',1e-06, 'AbsTol',1e-09) ;
    opts3 = odeset('Events', myeventfunc_heavy_purge_step, 'RelTol',1e-06, 'AbsTol',1e-09) ;
    opts4 = odeset('Events', myeventfunc_CnCDepress_step, 'RelTol',1e-06, 'AbsTol',1e-09) ;
    opts5 = odeset('Events', myeventfunc_light_purge_step, 'RelTol',1e-06, 'AbsTol',1e-09) ;

    %opts1 = odeset( 'JPattern', JacPressurization(N), 'RelTol', 1e-6, 'Events', CoCPres_event_fxn)       ;
    %opts2 = odeset( 'JPattern', JacAdsorption(N), 'RelTol', 1e-6)           ;
    %opts3 = odeset( 'JPattern', JacAdsorption(N), 'RelTol', 1e-6)           ;
    %opts4 = odeset( 'JPattern', Jac_CnCDepressurization(N), 'RelTol', 1e-6) ;
    %opts5 = odeset( 'JPattern', Jac_LightReflux(N), 'RelTol', 1e-6)         ;
%   
%% Begin simulating PSA cycle. This is skipped if the first constraint is violated.
%   Run the simulation until the change in the temperature, gas mole fraction
%   and CO2 molar loading is less than 0.5%. For the molar loading and the
%   mole fraction, the simulation also stops if the absolute change in the
%   state variable is less than 5e-5 and 2.5e-4 respectively. 

if constraints(1) == 0
    
    for i=1:250
        
    %% 1. Simulate CoCPressurization step
        press_time_count = tic;
        [t1, a] = ode15s(CoCPressurization_fxn, [0 tau_CoCPres], x0, opts1) ;
        [~] = toc(press_time_count);
        
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
        adsorption_time_count = tic;
        [t2, b] = ode15s(Adsorption_fxn, [0 tau_ads], x10, opts2) ;
        [~] = toc(adsorption_time_count);
        
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
        LightReflux_fxn = @(t, x) FuncLightReflux(t, x, Params, MPDParams, LN_Pi, Ene, state_rd) ;
        
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
        heavy_ref_time_count = tic;
        [t3, c] = ode15s(HeavyReflux_fxn, [0 tau_HR], x20, opts3) ;
        [~] = toc(heavy_ref_time_count);
        
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
        depress_time_count = tic;
        [t4, d] = ode15s(CnCDepressurization_fxn, [0 tau_CnCDepres], x30, opts4) ;
        [~] = toc(depress_time_count);
        
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
        light_ref_time_count = tic;
        [t5, e] = ode15s(LightReflux_fxn, [0 tau_LR], x40, opts5) ;
        [~] = toc(light_ref_time_count);
        
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
        
        HeavyReflux_fxn = @(t, x) FuncHeavyReflux(t, x, Params, MPDParams, LN_Pi, Ene, state_rd) ;
        
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
        if i >=10
            if CSS_states <= 1e-2 && abs(massBalance-1) <= 0.01
                break
            end      
        end

        if i == 250
            error('Computation did not complete within 700 iterations.');
        end

    end
    
    %% Process and Economic evaluation
    
    [purity, recovery, MB] = ProcessEvaluation(a, b, c, d, e, t1, t2, t3, t4, t5) ;
    

    %% Compile objectives and constraint violations
    con = recovery/MB - 0.9       ;
    if con < 0
        constraints(2) = abs(con) ;
    end
    
    switch type
        case 'ProcessEvaluation'
    
            objectives  = [0, 0]    ;

            objectives(1) = (1-purity)^2       ;
            objectives(2) = (1-recovery/MB)^2  ;
            
            con = purity - y_0            ;
            if con < 0
                constraints(3) = abs(con) ;
                %constraints(2) = 0        ;
            end
            
        case 'EconomicEvaluation'

            objectives = CO_capture_cost  ;
            
            con = purity - 0.99            ;
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
        
        [n_tot_CoCPres_HPEnd, n_CO2_CoCPres_HPEnd, ]                       = StreamCompositionCalculator(tau{1}, step{1}, 'HPEnd') ;
        [n_tot_ads_HPEnd, n_CO2_ads_HPEnd, ~]                              = StreamCompositionCalculator(tau{2}, step{2}, 'HPEnd') ;
        [n_tot_ads_LPEnd, n_CO2_ads_LPEnd, ~]                              = StreamCompositionCalculator(tau{2}, step{2}, 'LPEnd') ;
        [n_tot_HR_LPEnd, n_CO2_HR_LPEnd, ~]                                = StreamCompositionCalculator(tau{3}, step{3}, 'LPEnd') ;
        [n_tot_HR_HPEnd, n_CO2_HR_HPEnd, ~]                                = StreamCompositionCalculator(tau{3}, step{3}, 'HPEnd') ;
        [n_tot_CnCDepress_HPEnd, n_CO2_CnCDepres_HPEnd, ~]                 = StreamCompositionCalculator(tau{4}, step{4}, 'HPEnd') ;
        [n_tot_LR_LPEnd, n_CO2_LR_LPEnd, ~]                                = StreamCompositionCalculator(tau{5}, step{5}, 'LPEnd') ;
        [n_tot_LR_HPEnd, n_CO2_LR_HPEnd, ~]                                = StreamCompositionCalculator(tau{5}, step{5}, 'HPEnd') ;
	%   
    %% Calculate Purity, recovery and mass balance of the column

        purity        = ((n_tot_ads_LPEnd-n_CO2_ads_LPEnd)*(1-alpha)+(n_tot_HR_LPEnd-n_CO2_HR_LPEnd))/(n_tot_ads_LPEnd*(1-alpha)+n_tot_HR_LPEnd) ;
        recovery      = ((n_tot_ads_LPEnd-n_CO2_ads_LPEnd)*(1-alpha)+(n_tot_HR_LPEnd-n_CO2_HR_LPEnd))/((n_tot_ads_HPEnd-n_CO2_ads_HPEnd)+(n_tot_CoCPres_HPEnd-n_CO2_CoCPres_HPEnd)) ;

        comp_2_out = (n_tot_CnCDepress_HPEnd - n_CO2_CnCDepres_HPEnd)...
                     + (n_tot_ads_LPEnd - n_CO2_ads_LPEnd)...
                     + (n_tot_HR_LPEnd - n_CO2_HR_LPEnd)...
                     + (n_tot_LR_HPEnd - n_CO2_LR_HPEnd);

        comp_2_in = (n_tot_CoCPres_HPEnd - n_CO2_CoCPres_HPEnd)...
                     + (n_tot_ads_HPEnd - n_CO2_ads_HPEnd)...
                     + (n_tot_HR_HPEnd - n_CO2_HR_HPEnd)...
                     + (n_tot_LR_LPEnd - n_CO2_LR_LPEnd);

        mass_balance = comp_2_out/comp_2_in;
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
%   
%% PRESSURIZATION STEP EVENY FUNCTION
    function [value, isterminal, direction] = eventfunction_CoCpress_step(~, ~)
    % This function will stop the solution and throw an exception if the ode solution
    % takes more than the time specified.

        if toc(press_time_count) >= Maxpresstime
            value = 0;
            error("Maximum time for ode evaluation for Pressurization step reached. Terminating the solution!");
        else
            value = 1;
        end
        isterminal = 1;
        direction = 1;
    end
%
%% ADSORPTION STEP EVENT FUNCTION
    function [value, isterminal, direction] = eventfunction_adsorption_step(~, ~)
    % This function will stop the solution and throw an exception if the ode solution
    % takes more than the time specified.
        
        if (toc(adsorption_time_count)-Maxpurgetime) >= 0
            value = 0;
            error("Maximum time for ode evaluation for Adsorption step reached. Terminating the solution!")
        else
            value = 1;
        end
        isterminal = 1;
        direction = 0;
    end
%
%% HEAVY REFLUX EVENT FUNCTION
    function [value, isterminal, direction] = eventfunction_heavy_purge_step(~, ~)
    % This function will stop the solution and throw an exception if the ode solution
    % takes more than the time specified.
        
        if (toc(heavy_ref_time_count)-Maxpurgetime) >= 0
            value = 0;
            error("Maximum time for ode evaluation for Heavy reflux Purge_step reached. Terminating the solution!")
        else
            value = 1;
        end
        isterminal = 1;
        direction = 0;
    end
%
%% DEPRESSURIZATION STEP EVENT FUNCTION
    function [value, isterminal, direction] = eventfunction_CnCdepress_step(~, ~)
    % This function will stop the solution and throw an exception if the ode solution
    % takes more than the time specified.

        if toc(depress_time_count) >= Maxpresstime
            value = 0;
            error("Maximum time for ode evaluation for Depressurization step reached. Terminating the solution!");
        else
            value = 1;
        end
        isterminal = 1;
        direction = 1;
    end
%% LIGHT REFLUX FUNCTION
    function [value, isterminal, direction] = eventfunction_light_purge_step(~, ~)
    % This function will stop the solution and throw an exception if the ode solution
    % takes more than the time specified.
        if (toc(light_ref_time_count)-Maxpurgetime) >= 0
            value = 0;
            error("Maximum time for ode evaluation for Light Purge_step reached. Terminating the solution!")
        else
            value = 1;
        end
        isterminal = 1;
        direction = 0;
    end
end  
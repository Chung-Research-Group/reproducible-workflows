function q=Isotherm(y, P, T, isotherm_parameters)
%hypotheticalisotherm- Calculate the molar loadings of a two component mixture
%   Calculate the total loading of a two component mixture loading. It is
%   assumed that the isotherm is able to be represented by a dual site
%   langmuir freundlich competitive isotherm.
%   The isotherm can be in terms of Partial pressure [Pa]. 
%   Inputs:
%   y: mole fraction of component one [-]. It is assumed the mole fraction
%   of component 2 is 1-y
%   P: Total pressure of the gas [Pa]
%   T: Temperature of the gas [K]
%   isothermparams: Parameters for the isotherm.

    R=8.314; % J/mol/K
    
    q_sb_1=isotherm_parameters(1); % mol/kg
    q_sd_1=isotherm_parameters(2); % mol/kg
    b_1=isotherm_parameters(3); % 1/Pa
    d_1=isotherm_parameters(4); % 1/Pa
    n_b_1=isotherm_parameters(5);  % -
    n_d_1=isotherm_parameters(6);  % -
    delta_H_1=isotherm_parameters(7) ; % J/mol

    q_sb_2=isotherm_parameters(8); % mol/kg
    q_sd_2=isotherm_parameters(9); % mol/kg
    b_2=isotherm_parameters(10); % 1/Pa
    d_2=isotherm_parameters(11); % 1/Pa
    n_b_2=isotherm_parameters(12);  % -
    n_d_2=isotherm_parameters(13);  % -
    delta_H_2=isotherm_parameters(14) ; % J/mol

    P_1=y.*P.*exp(-delta_H_1/R.*(1./T-1/298));
    P_2=(1-y).*P.*exp(-delta_H_2/R.*(1./T-1/298));

    input_1=P_1;
    input_2=P_2;

    B_1=b_1.*input_1.^(1./n_b_1) ;
    B_2=b_2.*input_2.^(1./n_b_2) ;

    D_1=d_1.*input_1.^(1./n_d_1) ;
    D_2=d_2.*input_2.^(1./n_d_2) ;
    
    q1=q_sb_1.*B_1./(1+B_1+B_2) + q_sd_1.*D_1./(1+D_1+D_2) ;
    q2=q_sb_2.*B_2./(1+B_1+B_2) + q_sd_2.*D_2./(1+D_1+D_2) ;
	
	q = [q1, q2] ;

end
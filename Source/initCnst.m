function cnst = initCnst()
    cnst = struct(...    %
     'c',[], ...      %   speed of light in vacuum [m/s]
     'e',[], ...      %   electron charge [C]
     'eps0',[], ...   %   permittivity of vacuum [F/m]
     'mp',[]) ;       %   proton mass [kg]

    % various math
    cnst.pi       = 3.141592653589793;
    cnst.twopi    = 6.283185307179586;
    cnst.ln10     = 2.302585092994046;
    
    % various physical
    cnst.c        = 2.9979E+8;      %   speed of light in vacuum [m/s]
    cnst.cumps    = 299.79;         %     ... ditto [um/ps]
    cnst.e        = 1.6022E-19;     %   electron charge [C]
    cnst.mp       = 1.6726E-27;     %   proton mass [kg]
    cnst.eps0     = 8.8542E-12;     %   permittivity of vacuum
                                    %   [F/m]
    
    % various plasma
    cnst.vTe1eV = 0.419;    % VTe = [4.19e7*sqrt(Te)] in um/ps with Te
                            % in eV
    cnst.wpe         = 5.64e4;   % wpe = 5.64e4*sqrt(ne) in rad/s
    cnst.lamDebye    = 7.43e2;   % lamDebye = 7.43e2*sqrt(Te/ne) in
                                 % cm
    
    % default laser params
    cnst.lambda0 = 0.351e-6; % m
    cnst.omega0 = (cnst.c)*(cnst.twopi)/(cnst.lambda0); % 1/secs
    

function delta = chromatinOde(t,x,time,v, tf)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% delta = chromatinOde(t,x,ode_options,v)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Full chromatin ODE model - stiff system, designed to be solved using ode15 (see this function's help file to set
% options).  
% In phase 1 (v.PHASE == 1),simulation is run w/o stimulus, until convergence to initialize steady-state,
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% Set the perisistant variable that discrete-delay rxns refer to 

% persistent DELAY;
% if isempty(t)
%     % The sim function calls the ode function once without arguments to reset persistent variables
%     sz = 5000; % Size of total delay memory - if DELAYS.idx exceeds this, function will throw error
%     DELAY.t = zeros(sz,1);
%     DELAY.E0 = zeros(sz,1); % increase if getting errors
%     DELAY.idx     = 1; % index (starts at 1)
%     return;
% end


% Slice parameters, get previous concentrations
[v.PARAMS, v.SPECIES] = chromatinInitialize();
p = v.PARAMS;
delta = zeros(size(x));

tf = interp1(time, tf, t, 'linear');

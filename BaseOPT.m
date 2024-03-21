% $Date: 2012-10-01 10:51:29 +0200 (Mon, 01 Oct 2012) $
% $Revision: 9881 $

%% Construct faulted grid, and populate with layered permeability field  
mrstModule add ad-core ad-blackoil ad-props ad-fi optimization

% Define model ------------------------------------------------------------
nx = 20; ny = 20; nz = 4;
G = cartGrid([nx ny nz], [500, 500, 8]*meter);
G = computeGeometry(G);

c = G.cells.centroids;
[I, J, K] = gridLogicalIndices(G);

rock.perm  = max(30*sin(c(:,2)/25+.5*sin(c(:,1)/25))-10, .01)*200*milli*darcy;
rock.poro  = repmat(0.2, [G.cells.num, 1]);



% Wells and initial rates -------------------------------------------------
radius = .1;
simTime = 1600*day;
pv = poreVolume(G, rock);
injRate = 1*sum(pv)/simTime;
W = [];

% Rate controlled well
W = addWell(W, G, rock, 43:nx*ny:nx*ny*nz, ...
             'Type', 'rate', 'Val', injRate/2,... 
             'Comp_i', [1,0], ....
             'Name', 'I1'); 
             %'Type', 'bhp' , 'Val', 6000*psia);  
             
 W = addWell(W, G, rock, 58:nx*ny:nx*ny*nz, ...
             'Type', 'rate', 'Val', injRate/2,... 
             'Comp_i', [1,0], ....
             'Name', 'I2');
             %'Type', 'bhp' , 'Val', 6000*psia);    
             
W = addWell(W, G, rock, 358:nx*ny:nx*ny*nz, ...
             'Type', 'lrat', 'Val', -injRate/2,... 
             'Comp_i', [0,1], ....
             'Name', 'I3');
             %'Type', 'bhp' , 'Val', 6000*psia);  
     
W = addWell(W, G, rock, 343:nx*ny:nx*ny*nz, ...
             'Type', 'lrat', 'Val', -injRate/2,... 
             'Comp_i', [0,1], ....
             'Name', 'I4');
             %'Type', 'bhp' , 'Val', 6000*psia);   
figure,            
plotCellData(G, (rock.perm(:,1)));
plotWell(G, W), view(2)

%% Define base schedule
% We set up 4 control-steps each 160 days. We explicitly set shorter time
% steps during start.
ts = { [1 2 5 7 10 15 20 20 20 20 20 20]'*day, ...
                      repmat(160/7, 7, 1)*day, ...
                      repmat(160/7, 7, 1)*day, ...
                      repmat(160/7, 7, 1)*day, ...
                      repmat(160/7, 7, 1)*day, ...
                      repmat(160/7, 7, 1)*day, ...
                      repmat(160/7, 7, 1)*day};
       
numCnt = numel(ts);
[schedule.control(1:numCnt).W] = deal(W);
schedule.step.control = rldecode((1:7)', cellfun(@numel, ts));
schedule.step.val     = vertcat(ts{:});

%% Set fluid properties (slightly compressible oil)
fluid = initSimpleADIFluid('mu',    [.3, 5, 0]*centi*poise, ...
                           'rho',   [1000, 700, 0]*kilogram/meter^3, ...
                           'n',     [2, 2, 0]);
c = 1e-5/barsa;
p_ref = 300*barsa;
fluid.bO = @(p) exp((p - p_ref)*c);

%% Initialize model and run base-case simulation
model  = TwoPhaseOilWaterModel(G, rock, fluid);
state0 = initResSol(G, p_ref, [0 1]);
schedule_base = schedule;
[wellSols_base, states_base] = simulateScheduleAD(state0, model, schedule_base);

Time = cumsum(schedule.step.val);
plotWellSols(wellSols_base, Time);

%% Set prices ($/stb), discount rate (/year) and plot evolutiuon of NPV 
npvopts     = {'OilPrice',             100.0 , ...
               'WaterProductionCost',   7.0 , ...
               'WaterInjectionCost',    7.0 , ...
               'DiscountFactor',        0.05 };
v_base  = NPVOW(G, wellSols_base, schedule_base, npvopts{:});
v_base  = cell2mat(v_base);
t_base  = cumsum(schedule_base.step.val);
figure, plot(convertTo(t_base,day), cumsum(v_base), '-o', 'LineWidth', 2);
title('Base run evolution NPV'), ylabel('NPV, USD'), xlabel('days')

%% Set up box limits for scaling and define function evaluation handle
li = [  10, 300]/day;  % Injector limits  
lp = [-300, -10]/day;  % Producer limits 
scaling.boxLims = [li;li;lp;lp];  % control scaling  
scaling.obj     = sum(v_base);    % objective scaling    
% Get initial scaled controls 
u_base = schedule2control(schedule_base, scaling);
% Define objective function with above options
obj = @(wellSols, schedule, varargin)NPVOW(G, wellSols, schedule, varargin{:}, npvopts{:});
% Get function handle for objective evaluation
f = @(u)evalObjective(u, obj, state0, model, schedule, scaling);

%% Define linear equality and inequality constraints, and run optimization
% Constraints are applied to all control steps such for each step i, we
% enforce 
%       linEq.A*u_i   = linEq.b
%       linIneq.A*u_i = linIneq.b
%
% All rates should add to zero to preserve reservoir pressure (I1, I2, P1, P2)
linEq = struct('A', [1 1 1 1], 'b', 0);
% We also impose a total water injection constraint <= 500/day
linIneq = struct('A', [1 1 0 0], 'b', 500/day);  
% Constraints must be scaled!
linEqS   = setupConstraints(linEq,   schedule, scaling);
linIneqS = setupConstraints(linIneq, schedule, scaling);
% Run optimization with default options
[v, u_opt, history] = unitBoxBFGS(u_base, f, 'linEq', linEqS, 'linIneq', linIneqS);

%% Evaluate evolution of NPV for optimal schedule and compare with base
schedule_opt = control2schedule(u_opt, schedule, scaling);
[wellSols_opt, states_opt] = simulateScheduleAD(state0, model, schedule_opt);

%Time = cumsum(schedule.step.val);
plotWellSols(wellSols_opt, Time);
%%
v_opt  = NPVOW(G, wellSols_opt, schedule_opt, npvopts{:});
v_opt  = cell2mat(v_opt);
t_opt  = cumsum(schedule_opt.step.val);
figure, plot(convertTo([t_base, t_opt],day), cumsum([v_base, v_opt]), '-o', 'LineWidth', 2);
legend('Base', 'Optimal', 'Location', 'se'), xlabel('Time, days'), ylabel('NPV, USD')


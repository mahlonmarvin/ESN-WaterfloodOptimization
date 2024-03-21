%% Construct faulted grid, and populate with layered permeability field 
%%
grdecl = makeModel3([20, 20, 4], [500, 500, 8]*meter);
G = processGRDECL(grdecl);
G = computeGeometry(G);

% Set up permeability based on K-indices
[I, J, K] = gridLogicalIndices(G);

px       = 100*milli*darcy*ones(G.cells.num,1);
px(K==2) = 500*milli*darcy;
px(K==3) = 50*milli*darcy;

% Introduce anisotropy by setting K_x = 5*K_z.
perm = [px, px, 0.1*px];
rock = makeRock(G, perm, 0.2);

%% Set fluid properties (slightly compressible oil)
fluid = initSimpleADIFluid('mu',    [.3, 5, 0]*centi*poise, ...
                           'rho',   [1000, 700, 0]*kilogram/meter^3, ...
                           'n',     [2, 2, 0]);
c = 1e-5/barsa;
p_ref = 300*barsa;
fluid.bO = @(p) exp((p - p_ref)*c);

%% Define wells 
% We define two vertical producers and two vertical injectors. Rates
% correspond to one reservoir pore volume injected during 640 days
simTime = 2555*day;
pv = poreVolume(G, rock);
injRate = 1*sum(pv)/simTime;

% Place wells
[nx, ny] = deal(G.cartDims(1), G.cartDims(2));
ppos = [nx-2, 4; 4   , ny-4];
ipos = [4   , 5; nx-3, ny-3];
offset = 5;
W = [];
u1 = injRate/2;
u2 = injRate/2;
u3 = -u1;
u4 = -u2;

% Injectors            
W = verticalWell(W, G, rock, ipos(1,1), ipos(1,2), [], 'sign', 1,...
                'Name', 'I1', 'comp_i', [1 0], 'Val', injRate/2, 'Type', 'rate');

W = verticalWell(W, G, rock, ipos(2,1), ipos(2,2), [], 'sign', 1, ...
                'Name', 'I2', 'comp_i', [1 0], 'Val', injRate/2, 'Type', 'rate');
% Producers
W = verticalWell(W, G, rock, ppos(1,1), ppos(1,2), [], 'sign', -1, ...
                'Name', 'P1', 'comp_i', [0 1], 'Val', -injRate/2, 'Type', 'lrat');
W = verticalWell(W, G, rock, ppos(2,1), ppos(2,2), [], 'sign', -1, ...
                'Name', 'P2', 'comp_i', [0 1], 'Val', -injRate/2, 'Type', 'lrat');
figure,            
plotCellData(G, log(rock.perm(:,1)));
plotWell(G, W), view([1 1 1])
  
  % Create model-object of class TwoPhaseOilWaterModel
  model  = TwoPhaseOilWaterModel(G, rock, fluid);
  
  
  % Set up 10 control-steps each 2555 days (7 years)
time = 2555*day;           % period of simulation
n = 2555;
dt = time/n;
timesteps = repmat(dt, n, 1);
schedule = simpleSchedule(timesteps, 'W', W);

% model simulation
state0 = initResSol(G, 4500*barsa, [0.15, 0.85]);
[wellSols, states]         = simulateScheduleAD(state0, model, schedule);


% Plot Bottom hole Pressure for both wells
%li = [0 800]/day;
%lp = [100 800]*psia;
%plotSchedules(schedule, 'singleplot', true, 'boxconst', [li;li;lp;lp])
%%
plotWellSols(wellSols, time)

%%
% Compute NPV
d   = 0.05;    % yearly discount factor
ro  = 100;      % oil revenue/price ($/stb)
rwp =  7;      % water production handling costs ($/stb)
rwi =  7;      % water injection cost ($/stb) 

npvopts = {'OilPrice',             ro , ...
           'WaterProductionCost', rwp , ...
           'WaterInjectionCost',  rwi , ...
           'DiscountFactor',        d};
vals     = cell2mat(NPVOW(G, wellSols, schedule, npvopts{:}));

% Plot evolution of NPV and indicate peak value:
npv = cumsum(vals);
figure,  plot(time, cumsum(vals), '--b', 'LineWidth', 2);
title('Evolution of NPV [Naira]')
legend('base')
xlabel('Time')
ylabel('NPV')

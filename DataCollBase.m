%% Construct faulted grid, and populate with layered permeability field
%%
% Define model ------------------------------------------------------------
nx = 20; ny = 20; nz = 4;
G = cartGrid([nx ny nz]);
G = computeGeometry(G);

c = G.cells.centroids;
[I, J, K] = gridLogicalIndices(G);

rock.perm  = max(30*sin(c(:,2)/25+.5*sin(c(:,1)/25))-10, .01)*200*milli*darcy;
rock.poro  = repmat(0.2, [G.cells.num, 1]);


%% Set fluid properties (slightly compressible oil)
fluid = initSimpleADIFluid('mu',    [.3, 5, 0]*centi*poise, ...
                           'rho',   [1000, 700, 0]*kilogram/meter^3, ...
                           'n',     [2, 2, 0]);
c = 1e-5/barsa;
p_ref = 300*barsa;
fluid.bO = @(p) exp((p - p_ref)*c);

% Wells and initial rates -------------------------------------------------
radius = .1;
simTime = 2555*day;
pv = poreVolume(G, rock);
injRate = 1*sum(pv)/simTime;
W = [];
% Injectors (lower-left and upper-right)
ci(1) = 1;
ci(2) = G.cells.num;
for k  = 1:2
    W = addWell(W, G, rock, ci(k), 'Type' , 'rate', ...
                                   'Val'  , injRate, ...
                                   'Name' , sprintf('I%d', k), ...
                                   'comp_i', [1 0], ...
                                   'Sign' , 1);
end
% Producers (upper-left and -right)
cp(1) = G.cartDims(1);
cp(2) = 1 + (G.cartDims(2)-1)*G.cartDims(1);
for k  = 1:2
    W = addWell(W, G, rock, cp(k), 'Type', 'rate', ...
                                   'Val' , -injRate, ...
                                   'Name', sprintf('P%d', k), ...
                                   'comp_i', [0 1], ...
                                   'Sign', -1);
end
figure,            
plotCellData(G, (rock.perm(:,1)));
plotWell(G, W), view(2)
  
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

% Compute NPV
d   = 0;    % yearly discount factor
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

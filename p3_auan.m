%% 2D ISING MODEL
% Assumes square lattice

clear

% Params
lattice = 2.^[3 4 5]; J = 1;
% M_therm = 1000; M = 1200;
M_therm = 1000; M = 10;

% Number of runs taken as the grand average
grand_average = 16;

% For readability
MC = 1;

% Temperature interval struct
dT = 0.1; dT_crit = 0.01;
T = struct('min1',1.7,'max1',2-dT, ...
  'crit1', 2, 'crit2', 2.5-dT_crit, ...
  'min2', 2.5, 'max2', 3.3);

% Flip from high to low T
% due to cheaper therm phase
kBTJ = fliplr( ...
    [T.min1:dT:T.max1 ...
    T.crit1:dT_crit:T.crit2 ...
    T.min2:dT:T.max2])/J;


% Possible energy states of nearest nbrs
E_pos = -4:2:4;

% Alloc
order = zeros(length(kBTJ), length(lattice));
susc = zeros(length(kBTJ), length(lattice));

% ***************** Unused *********************
% err_MAT = zeros(length(kBTJ), length(lattice));
% mctau_MAT = zeros(length(kBTJ), length(lattice));
% c_MAT= zeros(length(kBTJ), length(lattice));

% Constant loop variable is a criterion of parpool
% Total steps in T
T_steps = length(kBTJ);
n_lattices = length(lattice);

% Parallel for loop
parfor ga = 1:grand_average
  for sz = 1:n_lattices
    % Lattice size for reuse
    N = lattice(sz);

    % Simulation info
    disp(['Lattice size N = ' num2str(N) ' on T = [' num2str(T.min1) ',' num2str(T.max2) ']']);
    disp([num2str(M*N*N) ' sweeps over ' num2str(MC*N*N) ' MC steps']);
    disp('********* Simulation started *********');
    p1 = datetime(now,'ConvertFrom','datenum');
    disp(p1);

    % Random spin matrix
    S = lattice_init(N);

    % Temperature loop
    for i = 1:T_steps
      t_iter1 = datetime(now,'ConvertFrom','datenum');
      % Populate Boltzmann-Gibbs factors.
      % Easily modified to a Glauber step.
      h = min(1,exp(-2*E_pos./kBTJ(i)));

      % Thermalization phase
      for ij = 1:M_therm*N*N
        S = metropolis(S,N,h);
      end

      % Init vars
      c = 0; rm = 0; rm2 = 0;
      rm1 = abs(sum(sum(S)));

      % MC phase
      for im = 1:M*N*N
        for ij = 1:MC*N*N
          S = metropolis(S,N,h);
        end

        % Accumulate measurements
        rm0 = abs(sum(sum(S)));
        rm=rm+rm0;
        rm2=rm2+rm0*rm0;
        c=c+rm0*rm1;
        rm1=rm0;

      end
      % Averages
      rm=rm/(M*N*N);
      rm2=rm2/(M*N*N)-rm*rm;

      % ***************** Unused *********************
      % c=(c/M-rm*rm*N*N)/rm2*N*N;
      % if (c~=1.0)
      %   tau=c/(1-c);
      % end

      % Get iteration time
      t_iter2 = datetime(now,'ConvertFrom','datenum');
      disp(['T = ' num2str(kBTJ(i)) ' [' num2str(i) '/' num2str(T_steps) '] of run ' num2str(ga) ' took ' char(t_iter2-t_iter1)])
      %err=sqrt(rm2*(2*tau+1)/M);

      order(i,sz,ga) = rm;
      susc(i,sz,ga) = N*rm2/kBTJ(i);

      % ***************** Unused *********************
      % err_MAT(i,sz) = err;
      % mctau_MAT(i,sz) = MC*tau;
      % c_MAT(i,sz) = c;

    end
    p2 = datetime(now,'ConvertFrom','datenum');
    disp('********* Simulation finished *********');
    disp(p2-p1);

  end
end

% Avg of 3D tensor
order_avg = mean(order, 3);
susc_avg = mean(susc, 3);
save('order.mat', 'order');
save('susc.mat', 'susc');

subplot(2,1,1)
plot(kBTJ,(order_avg./lattice.^2), '.-')
hold on
legend
title('Order')
subplot(2,1,2)
legend
plot(kBTJ,susc_avg./lattice.^2, '.-')
hold on
title('Susceptibilty')
legend

function nbrs = site_nbrs(S, site)
  % Get nbrs given site
  N = size(S,1);
  % Linear idx to matrix idx
  [i,j] = ind2sub([N N], site);
  % Get nbrs state
  nbrs = S(mod(i-2, N)+1,j) ...
    + S(mod(i, N)+1,j) ...
    + S(i,mod(j, N)+1) ...
    + S(i,mod(j-2,N)+1);
end

function S = lattice_init(N)
  S = (-1).^randi(2,N,N);
end

function S = metropolis(S,N,h)
  % Random site in S
  site = randi(N*N);
  dE = S(site)*site_nbrs(S,site);
  % Metropolis step
  if dE <= 0 || rand() < h(((dE+4)/2+1))
    S(site) = -S(site);
  end
end

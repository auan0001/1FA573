%% 2D ISING MODEL
% Assumes square lattice

clear

% Params
% lattice = 2.^[3 4 5]; J = 1;
lattice = 2.^[4]; J = 1;
% M = 100
M_therm = 1000; M = 10000;

% Number of runs taken as the grand average
grand_avg = 64;

% For readability
MC = 1;

% Temperature interval struct
dT = 0.05; dT_crit = 0.005;
T = struct('min1',1.7,'max1',2.05-dT, ...
  'crit1', 2.05, 'crit2', 2.7-dT_crit, ...
  'min2', 2.7, 'max2', 3);

% Flip from high to low T
% due to cheaper therm phase
kBTJ = fliplr( ...
    [T.min1:dT:T.max1 ...
    T.crit1:dT_crit:T.crit2 ...
    T.min2:dT:T.max2])/J;


% Possible energy states of nearest nbrs
E_pos = -4:2:4;

% Alloc
% order = zeros(length(kBTJ), length(lattice), grand_avg);
% susc = zeros(length(kBTJ), length(lattice), grand_avg);
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
parfor ga = 1:grand_avg
  for sz = 1:n_lattices
    % Lattice size for reuse
    N = lattice(sz);

    % Simulation info
    % disp(['Lattice size N = ' num2str(N) ' on T = [' num2str(T.min1) ',' num2str(T.max2) ']']);
    disp([num2str(M*N*N) ' sweeps over ' num2str(MC*N*N) ' MC steps']);
    % disp('********* Simulation started *********');
    p1 = datetime(now,'ConvertFrom','datenum');
    disp(p1);

    % Random spin matrix
    S = lattice_init(N);

    % Temperature loop
    for i = 1:T_steps
      t_iter1 = datetime(now,'ConvertFrom','datenum');
      % Populate Boltzmann-Gibbs factors
      % for Metropolis-Hastings choice
      h = min(1,exp(-2*E_pos./kBTJ(i)));
      % Glauber Heat-Bath choice
      % g = 1/(1 + exp(2*E_pos./kBTJ(i)));

      % Thermalization phase
      for ij = 1:M_therm*N*N
        S = metropolis(S,N,h);
      end

      % Init vars
      c = 0; m = 0; m2 = 0; E = 0; E2 = 0; E0 = 0;
      m1 = abs(sum(sum(S)));
      nbrs = lattice_nbrs(S);
      E1 = lattice_E(S, nbrs);

      % MC phase
      for im = 1:M*N
        for ij = 1:MC*N
          S = metropolis(S,N,h);
        end

        % Accumulate measurements
        E0 = lattice_E(S, nbrs);
        m0 = abs(sum(sum(S)));
        m = m + m0;
        E = E + E0;
        m2 = m2 + m0*m0;
        E2 = E2 + E0*E0;

        % c=c+m0*m1;
        % m1=m0;

      end
      % Averages
      m = m/(M*N);
      E = E/(M*N);
      m2 = m2/(M*N)-m*m;
      E2 = E2/(M*N)-E*E;

      % ***************** Unused *********************
      % c=(c/M-m*m*N*N)/m2*N*N;
      % if (c~=1.0)
      %   tau=c/(1-c);
      % end

      % Get iteration time
      t_iter2 = datetime(now,'ConvertFrom','datenum');
      disp(['T = ' num2str(kBTJ(i)) ' [' num2str(i) '/' num2str(T_steps) '] of run ' num2str(ga) ' took ' char(t_iter2-t_iter1)])
      % err=sqrt(rm2*(2*tau+1)/M);

      % Put measurement matrices as pages in rank 3 tensor
      order(i,sz,ga) = m;
      susc(i,sz,ga) = m2/kBTJ(i);
      energy(i,sz,ga) = E;
      heat(i,sz,ga) = E2/(kBTJ(i)^2);

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

% Avg of tensor in dim 3
order_avg = mean(order, 3);
susc_avg = mean(susc, 3);
heat_avg = mean(heat, 3);

save('order.mat', 'order_avg');
save('susc.mat', 'susc_avg');
save('heat.mat', 'heat_avg');
save('temp.mat', 'kBTJ');

subplot(3,1,1)
hold on
plot(kBTJ,(order_avg./lattice.^2), '-')
legend
title('Order')
subplot(3,1,2)
hold on
legend
plot(kBTJ,susc_avg, '-')
title('Susceptibilty')
legend
subplot(3,1,3)
hold on
plot(kBTJ,heat_avg./lattice.^2, '-')
legend

function E = lattice_E(S, nbrs)
  % Energy Hamiltonian of the whole matrix
  E = -2*abs(sum(sum(S.*nbrs)));
end

function nbrs = lattice_nbrs(S)
  % Circular shift to compute neighbours
  % Similar to MPI_Cart_shift or np.roll
  nbrs = circshift(S,[0 1]) ...
    + circshift(S, [0 -1]) ...
    + circshift (S,[1 0]) ...
    + circshift(S, [-1 0]);
end

function nbrs = site_nbrs(S, N, site)
  % Get nbrs given site
  % N = size(S,1);
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
  dE = S(site)*site_nbrs(S,N,site);
  % Metropolis step
  if dE <= 0 || rand() < h(((dE+4)/2+1))
    S(site) = -S(site);
  end
end

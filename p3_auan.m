clear
% Lattice and sample size
N = 8;
N_samples = 1000;

% Init spin matrix
S = lattice_init(N);

% Length of phases
therm_phase = N_samples;
sweep_phase = N_samples;

% kBT/J is our unit (from hints)
J = 1;
unit_min = 1; 
unit_max = 5; 
kBT_J = linspace(unit_min,unit_max,N_samples)./J; 

% TODO: ising function as driver function
for i = 1:length(kBT_J)
  % Thermalization phase for each k_B*T/J
  for j = 1:therm_phase
    % nbrs = lattice_nbrs(S);
    S = metropolis(S, lattice_nbrs(S), kBT_J(i));
  end

  % Sweep phase for each k_B*T/J
  for j = 1:sweep_phase
    nbrs = lattice_nbrs(S);
    S = metropolis(S, nbrs, kBT_J(i));
    % E(i) = lattice_E(S, nbrs);
    M(i) = lattice_M(S);
  end
end

win_min = 50;
win_step = 1;
win_max = 200;
win_tol = 1e-2;

% Interested in magnitude of M this time
% Probably need to be normalized in some way
M = abs(M);

% TODO: Vectorize and test
% Find where the absolute difference of the sums flattens out
% To obtain a good rolling mean
% Small window size <-> cheap, unsmooth, unrealistic
% Too large window size <-> expensive, smooth, unrealistic
% Optimal window size <-> cheap, smooth, realistic
win_sz = win_min:win_step:win_max;
for sz = 1:length(win_sz)
  M_smooth = movmean(M, win_sz(sz));
  sum_abs_diff(sz) = sum(abs(M_smooth - M));
end

% Find 'optimal' i.e cheap window size
win_opt = find(diff(sum_abs_diff) <= win_tol, 1, 'last');
% If we cannot find a good size, use the cheapest
if isempty(win_opt)
  % Use smallest size
  win_opt = min(win_sz);
  disp(['WARNING! Did not find optimal window size'])
  disp(['Using fallback size: ' num2str(win_opt)])
else
  disp(['Window size: ' num2str(win_sz(win_opt))])
end

% Moving mean of M
M_mov_mean = movmean(M, win_opt);

figure;
% plot(kBT_J,M,'.')
hold on
grid on
title(['Lattice size N =' num2str(N)])
plot(kBT_J,M_mov_mean,'-', 'Linewidth',2)
xlabel('{k_bT}/J')
% ylabel('<E>')
ylabel('M [magnitude]')
hold off

function S = metropolis(S, nbrs, kBT_J)
  % Assumes square spin lattice
  N = size(S,1);
  % Random lattice site
  site = randi(N*N);
  % Compute all nbrs once every iter
  % Compute dE and check spin flip
  dE = 2*S(site)*nbrs(site); 
  spin_flip_prob = exp(-dE/kBT_J);
  if dE <= 0 || rand() <= spin_flip_prob
    S(site) = -S(site); 
  end
end

function spin = lattice_init(N)
  % Assumes square spin lattice
  spin = (-1).^randi(2,N,N);
end

function E = lattice_E(S, nbrs)
  % Energy Hamiltonian of the whole matrix
  E = -mean(S.*nbrs,'all');
end

function M = lattice_M(S)
  % Energy Hamiltonian of the whole matrix
  M = sum(sum(S));
end

function nbrs = lattice_nbrs(S)
  % Circular shift to compute neighbours
  % Similar to MPI_Cart_shift or np.roll
  nbrs = circshift(S,[0 1]) ...
    + circshift(S, [0 -1]) ...
    + circshift (S,[1 0]) ...
    + circshift(S, [-1 0]);
end

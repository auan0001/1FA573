clear
% Lattice and sample size
N = [8 16 32];
N_samples = 200;

% Init spin matrix

% Length of phases
therm_phase = 1000;

% kBT/J is our unit (from hints)
J = 1;
unit_min = 0.1; 
unit_max = 5; 

% From high to low temp gives efficient therm phase
kBT_J = linspace(unit_max,unit_min,100)./J; 

for sz = 1:length(N)
  S = lattice_init(N(sz));
  M = zeros(1,length(kBT_J));
  tic
  % TODO: ising function as driver function
  for i = 1:length(kBT_J)
    M_MC = 0;
    % Thermalization phase for each k_B*T/J
    for j = 1:therm_phase
      S = metropolis(S, kBT_J(i));
    end

    % Sweep phase for each k_B*T/J
    for j = 1:N_samples
      for k = 1:N(sz)*N_samples
        S = metropolis(S, kBT_J(i));
        % Accumulate magnetization
        M_MC = M_MC + lattice_M(S);
      end
       % Compute average after every MC step
      M(i) = M_MC/(N(sz)*N_samples^2);
    end
  end
  toc

figure(1);
plot(kBT_J,M,'.')
hold on
grid on
title(['Lattice sizes N =' num2str(N)])
% plot(kBT_J,M_mov_mean,'-', 'Linewidth',2)
xlabel('{k_bT}/J')
% ylabel('<E>')
ylabel('M [magnitude]')
end
legend('N=8','N=16','N=32')
hold off

function S = metropolis(S, kBT_J)
  % Assumes square spin lattice
  N = size(S,1);
  % Random lattice site
  site = randi(N*N);
  % Compute dE using nbrs and check spin flip
  dE = 2*S(site)*site_nbrs(S,site);
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
  % Absolute magnetization
  M = abs(mean(S, 'all'));
end

function nbrs = lattice_nbrs(S)
  % Circular shift to compute neighbours
  % Similar to MPI_Cart_shift or np.roll
  nbrs = circshift(S,[0 1]) ...
    + circshift(S, [0 -1]) ...
    + circshift (S,[1 0]) ...
    + circshift(S, [-1 0]);
end

function nbrs = site_nbrs(S, site)
  % Get nbrs given site
  N = size(S,1);
  [i,j] = ind2sub([N N], site);
  n = mod(i-2, N)+1;
  s = mod(i, N)+1;
  w  = mod(j-2,N) + 1;
  e = mod(j, N)+1;
  nbrs = S(n,j)+S(s,j)+S(i,e)+S(i,e);
end

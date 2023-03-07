clear

N = 32;

% Init spin matrix
S = lattice_init(N);

% Temperature and spin interaction
J = 1;

% Length of thermalisation phase
therm_phase = 300;
sweep_phase = 300;

% kBT/J
kBT_J = linspace(0,3,50)./J; 


for i = 1:length(kBT_J)
  for j = 1:therm_phase
    S = metropolis_hastings(S,kBT_J(i));
  end

  for j = 1:sweep_phase
    S = metropolis(S,kBT_J(i));
    nbrs = lattice_nbrs(S);
    E(i) = lattice_energy(S, nbrs);
  end
end

plot(kBT_J,E,'*')

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
  % Energy mean of whole matrix
  E = mean(S.*nbrs,'all'));

function nbrs = lattice_nbrs(S)
  % Circular shift to compute neighbours
  % Similar to MPI_Cart_shift
  nbrs = circshift(S,[0 1]) ...
    + circshift(S, [0 -1]) ...
    + circshift (S,[1 0]) ...
    + circshift(S, [-1 0]);
end

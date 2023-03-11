clear

% lattice = 2.^[4 4 5];
lattice = 2.^[3 4 5];

M_therm = 1000;
% M_therm = 1;
% M = 8192;
M = 100;
MC = 1;

T_max = 3.5;
T_min = 1;
T_steps = 100;

J = 1;
kBTJ = linspace(T_max,T_min,T_steps)/J;
E_pos = -4:2:4;

rm_MAT = zeros(length(kBTJ), length(lattice));
rm2_MAT = zeros(length(kBTJ), length(lattice));
err_MAT = zeros(length(kBTJ), length(lattice));
mctau_MAT = zeros(length(kBTJ), length(lattice));
c_MAT= zeros(length(kBTJ), length(lattice));

for sz = 1:length(lattice)
  N = lattice(sz);

  disp(['Lattice size N = ' num2str(N) ' on T = [' num2str(T_max) ',' num2str(T_min) ']']);
  disp([num2str(M*N*N) ' sweeps over ' num2str(MC*N*N) ' MC steps']);
  disp('********* Simulation started *********');
  p1 = datetime(now,'ConvertFrom','datenum');
  disp(p1);

  S = lattice_init(N);
  for i = 1:length(kBTJ)
    % Populate Boltzmann-Gibbs factors
    h = min(1,exp(-2*E_pos./kBTJ(i)));
    % Thermalization phase
    for ij = 1:M_therm*N*N
      site = randi(N*N);
      dE = S(site)*site_nbrs(S,site);
      if dE <= 0 || rand() < h(((dE+4)/2+1))
        S(site) = -S(site);
      end
    end
    % Init vars
    c = 0;
    rm = 0;
    rm2 = 0;
    rm1 = abs(sum(sum(S)));
    % Measuring phase
    for im = 1:M*N*N
      for ij = 1:MC*N*N
        site = randi(N*N);
        dE = S(site)*site_nbrs(S,site);
        if dE <= 0 || rand() < h(((dE+4)/2+1))
          S(site) = -S(site);
        end
      end
      rm0 = abs(sum(sum(S)));
      rm=rm+rm0;
      rm2=rm2+rm0*rm0;
      c=c+rm0*rm1;
      rm1=rm0;
    end
    rm=rm/(M*N*N);
    rm2=rm2/(M*N*N)-rm*rm;
    c=(c/M-rm*rm*N*N)/rm2*N*N;
    if (c~=1.0)
      tau=c/(1-c);
    end
    disp(['T = ' num2str(kBTJ(i)) ' [' num2str(i) '/' num2str(T_steps) ']'])
    % disp('Writing')
    err=sqrt(rm2*(2*tau+1)/M);
    T_MAT(i,sz) = kBTJ(i);
    rm_MAT(i,sz) = rm;
    rm2_MAT(i,sz) = rm2/kBTJ(i);
    err_MAT(i,sz) = err;
    mctau_MAT(i,sz) = MC*tau;
    c_MAT(i,sz) = c;
  end
  p2 = datetime(now,'ConvertFrom','datenum');
  disp('********* Simulation finished *********');
  disp(p2-p1);
  save('T.mat', 'T_MAT');
  save('rm.mat', 'rm_MAT');
  save('rm2.mat', 'rm2_MAT');
end

subplot(2,1,1)
plot(kBTJ,(rm_MAT./lattice.^2), '.-')
legend
title('Order')
subplot(2,1,2)
legend
plot(kBTJ,rm2_MAT./lattice.^2, '.-')
title('Susceptibilty')
legend

function nbrs = site_nbrs(S, site)
  % Get nbrs given site
  N = size(S,1);
  [i,j] = ind2sub([N N], site);
  nbrs = S(mod(i-2, N)+1,j) ...
    + S(mod(i, N)+1,j) ...
    + S(i,mod(j, N)+1) ...
    + S(i,mod(j-2,N)+1);
end

function spin = lattice_init(N)
  % Assumes square spin lattice
  spin = (-1).^randi(2,N,N);
end

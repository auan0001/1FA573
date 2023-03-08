% FINAL PROJECT IN COMPUTATIONAL PHYSICS
% PROBLEM 3.

clear all 
close all 

N_samples = 1000; 
L = 8; 

% kBT/J is our unit (from hints)
J = 1;
unit_min = 1; 
unit_max = 4; 
kBT_J = linspace(unit_min,unit_max,N_samples)./J; 

test_results = zeros(N_samples, 1); 
tic
for i = 1:N_samples
    test_results(i) = ising(L, kBT_J(i)); 
end 

toc 

hold on 
plot(kBT_J, test_results, '.')
title('Lattice sizes N = 8x8', 'Fontsize',20)
xlabel('{k_bT}/J', 'Fontsize',20)
ylabel('M [magnitude]', 'Fontsize',20)
grid on
hold off

% function returns the average magnetization of an L-dim Ising system 
function M = ising(L, kBT_J)

%r = r_generate(kBT_J); % generates an r-matrix
S = (-1).^randi(2, L); % generating a pseudo-random matrix with values 1 and 2
N_thermal =  2000; % number of thermal equalization steps 
mag = 0; % pre-allocation
N_runs = 1e2; % number of sampling/ averaging runs 
N_sims = 100; % number of runs between each sampling 

% performing thermal equalization of the matrix 
for i = 1:N_thermal
    S = metropolis(S, L, kBT_J); 
end 

for i = 1:N_runs 
    
    for j = 1:N_sims
        S = metropolis(S, L, kBT_J);
        mag = mag + sum(S, 'all'); 
    end 
    
end 

M = abs(mag/(N_runs*N_sims))/L^2; 

end 

% implementation of the metropolis algorithm 
function S_new = metropolis(S, L, kBT_J)

% "typewriter algorithm" that goes through S and flips spins 
for i = 1:L
    for j = 1:L
        %S(i,j) = flip(S(i,j), spin_sum(S, i, j, L), r);
        S(i,j) = flip_alternative(S(i,j), spin_sum(S, i, j, L), kBT_J); 
    end 
end
S_new = S; 
end

% function returns the neighbour-sums of a point (i,j)
function tot = spin_sum(S, i, j, L)

% we have to do a bunch of tests to impose harmonic BC's 
if i == 1  % if we are at the top 
	upp = S(L,j); 
else
    upp = S(i-1,j); 
end 

if i == L  % if we're at the bottomn 
    low = S(1,j); 
else
    low = S(i+1,j); 
end 

if j == 1 % if we're at LHS
    left = S(i,L); 
else
    left = S(i,j-1); 
end 

if j == L % if we're at RHS
    right = S(i,1); 
else
    right = S(i,j+1); 
end 

tot = upp + low + left + right; 

end 

% function returns new spin-state of a point. Biplab's algorithm 
function new_spin = flip_alternative(state, spin_sum, kBT_J)

J = 1; % same as above 

dE = 2*J*state*spin_sum; % change in energy 

if dE < 0
    new_spin = state*(-1); 
else 
    r = rand(); %generating a random integer 
    
    if r < exp(-dE/kBT_J)
        new_spin = state*(-1);
    else 
       new_spin = state; 
    end 
end 
end 

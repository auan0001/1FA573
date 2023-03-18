#include <NumCpp.hpp>
#include <iostream>

// Macros for readability
#define SITE s[0],s[1]
#define SZ N*N
#define n1 (s[0]+1)%N, s[1]
#define n2 s[0], (s[1]+1)%N
#define n3 (s[0]-1)%N, s[1]
#define n4 s[0], (s[1]-1)%N
#define idx(dH) (dH+4)/2
#define avgnorm(x) x/((double)M*N*N)
#define norm(x) x/((double)N*N)

// Global constants
const int M = 2000; 
const int therm = 1000;
const int MC = 100;

// Metropolis-Hastings
void metropolis(nc::NdArray<int>& S, int& N, nc::NdArray<double>& h);

// Hamiltonian
double H(nc::NdArray<int>& S, int& N);

int main (int argc, char *argv[])
{
  std::ofstream out;
  std::ifstream in;
  out.open("run.txt");
  std::cout << "GO!" << std::endl;

  struct {
    const double max = 2.3;
    const double min = 2.2;
    const double step = 15;
  } T;

  nc::random::seed(666);

  // Measurements
  double m0 = 0;
  double m = 0;
  double E = 0;
  double E2 = 0;
  double m2 = 0;
  double m4 = 0;
  double chi = 0;
  double CB = 0;
  double U = 0;
  int E0 = 0;

  // Lattice sizes
  // auto lattices = nc::NdArray<int>{8,16,32};
  // auto lattices = nc::NdArray<int>{4,8,16};
  auto lattices = nc::NdArray<int>{8,16,32};

  // To compute Boltzmann-Gibbs factors
  auto h = nc::NdArray<double>(1,5);
  auto one = nc::ones<double>(1,5);

  // Temperature from high to low
  auto kBTJ = nc::fliplr(nc::linspace(T.min, T.max, T.step));

  // Set of possible energies in nbrhood
  auto E_pos = nc::linspace<double>(-4.0, 4.0, 5.0);

  // order lambda
  auto order = [](auto S, auto N){return nc::abs(nc::sum<int>(S).item())/((double)SZ);};

  // Temperature loop
  for (auto& N : lattices) {
    // Init random matrix
    auto S = nc::random::choice<int>({-1,1},SZ).reshape(N,N);
    for (auto& T : kBTJ) {
      // Populate G-B
      h = nc::minimum(one, nc::exp(-2.0*E_pos/T));


      // Thermalization
      for (size_t i = 0; i < therm*SZ; i++) {
          metropolis(S, N, h);
      }

      // Reset measurements
      m = 0; m2 = 0; chi = 0; E = 0; E2 = 0;

      // Simulation and measurements
      for (size_t i = 0; i < M; i++) {
        for (size_t j = 0; j < MC*SZ; j++) {
          metropolis(S, N, h);
        }

        // Order param
        m0 = order(S,N);
        // Energy
        // E0 = H(S,N);

        // Accumulate
        m += m0;
        m2 += m0*m0;
        m4 += m0*m0*m0*m0;
        // E += E0;
        // E2 += E0*E0;
      }
      // Averages over M samples
      m = m/M;
      m2 = m2/M;
      m4 = m4/M;
      // E = E/M;
      E2 = E2/M;

      chi = SZ*(m2-m*m)/(T); // Susceptibilty
      // CB = (E2-E*E)/(T*T); // Specific heat 
      U = 1-m4/(3*(m2*m2)); // Binder's parameter

      // Write to file
      out << m << ' ' << chi << ' ' << CB << ' ' << U << ' ' << T << "\n";
    }
  }
  return 0;
}

// Metropolis-Hastings
void metropolis(nc::NdArray<int>& S, int& N, nc::NdArray<double>& h) {
  nc::NdArray<int>s = nc::random::randInt({1,2},N);
  double dH = S(SITE)*(S(n1)+S(n2)+S(n3)+S(n4));
  if (dH <= 0 || nc::random::rand<double>() < h[idx(dH)]) 
    S(SITE) = -S(SITE);
}

// Hamiltonian
double H(nc::NdArray<int>& S, int& N) {
  auto nbrs = nc::roll(S, 1, nc::Axis::COL)
    + nc::roll(S, -1, nc::Axis::COL)
    + nc::roll(S, 1, nc::Axis::ROW)
    + nc::roll(S, -1, nc::Axis::ROW);
  return -nc::abs(nc::sum(nc::matmul(S,nbrs))).item()/((double)SZ);
}

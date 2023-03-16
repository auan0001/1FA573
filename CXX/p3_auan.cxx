#include <NumCpp.hpp>
#include <iostream>

// Indexing macros for readability
#define SITE s1,s2
#define N1 (s1+1)%N, s2
#define N2 s1, (s2+1)%N
#define N3 (s1-1)%N, s2
#define N4 s1, (s2-1)%N
#define idx(dH) (dH+4)/2

int main (int argc, char *argv[])
{
  std::ofstream out;

  std::cout << "GO!" << std::endl;

  struct
  {
    const double max = 5;
    const double min = 0.1;
    const double step = 40;

  } T;

  // Energy change in Hamiltonian
  double dH = 0;

  nc::random::seed(666);
  // const int N = 32;
  const int M = 1000;
  const int therm = 1000;
  const int MC = 1;

  // Random sites
  int s1;
  int s2;


  // Lattice sizes
  // auto lattices = nc::NdArray<int>{8,16,32};
  auto lattices = nc::NdArray<int>{8};

  // Boltzmann-Gibbs factors
  auto h = nc::NdArray<double>(1,5);
  auto one = nc::ones<double>(1,5);

  // Temperature array
  auto kBTJ = nc::fliplr(nc::linspace(T.min, T.max, T.step));

  // Set of possible nbr energies
  auto E_pos = nc::linspace<double>(-4.0, 4.0, 5.0);
  for (auto& N : lattices) {
    out.open("lattice");
    auto S = nc::random::choice<int>({-1,1},N*N).reshape(N,N);
    for (auto& T : kBTJ) {
      // Populate G-B
      h = nc::minimum(one, nc::exp(-2.0*E_pos/T));

      // Reset measurements
      double m = 0; int m0 = 0;

      // Thermalization
      for (size_t i = 0; i < therm; i++) {
        // Init matrix
        s1 = nc::random::randInt(N); s2 = nc::random::randInt(N);
        dH = S(s1,s2)*(S(N1)+S(N2)+S(N3)+S(N4));
        if (dH <= 0 || nc::random::rand<double>() < h[idx(dH)]) {
          S(SITE)=-S(SITE);
        }
      }

      // Simulation and measurements
      for (size_t i = 0; i < M*N*N; i++) {
        for (size_t j = 0; j < MC*N*N; j++) {
          s1 = nc::random::randInt(N); s2 = nc::random::randInt(N);
          dH = S(s1,s2)*(S(N1)+S(N2)+S(N3)+S(N4));
          if (dH <= 0 || nc::random::rand<double>() < h[idx(dH)]) {
            S(SITE)=-S(SITE);
          }
        }

        // Reduce and accumulate
        m0 = nc::abs(nc::sum<int>(S).item());
        m += m0/((double)N*N);
      }
      // Normalize and write to file
      out << m/((double)N*N*M) << ",";
    }
  }
  return 0;
}

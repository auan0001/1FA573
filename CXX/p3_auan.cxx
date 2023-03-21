#include <NumCpp.hpp>
#include <iostream>

// Macros
#define SITE s[0],s[1]
#define SZ N*N
#define n1 (s[0]+1)%N, s[1]
#define n2 s[0], (s[1]+1)%N
#define n3 (s[0]-1)%N, s[1]
#define n4 s[0], (s[1]-1)%N
#define idx(dH) (dH+4)/2
#define avgnorm(x) x/((double)M*N*N)
#define norm(x) x/((double)N*N)

// Columns
#define TEMP 0
#define ORDER 1
#define CHI 2
#define CB 3
#define U 4

// Global constants
const int M = 200000; 
const int therm = 1000;
const int MC = 1;

// Metropolis-Hastings
void metropolis(nc::NdArray<int>& S, int& N, nc::NdArray<double>& h);

// Hamiltonian
double H(nc::NdArray<int>& S, int& N);

// Print data columns
void tofile(nc::NdArray<double>& measurement, std::string file);

int main (int argc, char *argv[]) {
  if (argc < 3 || argc >= 4) {
    std::cout << "USAGE: <lattice dim> <file>" << std::endl;
    return EXIT_FAILURE;
  }

  // Input args
  int N = atoi(argv[1]);
  std::string file = argv[2];

  // Temperature steps
  struct {
    const double max = 5;
    const double min = 0.1;
    const double step = 80;
  } T;

  nc::random::seed(1337);

  // Measurements
  double m0 = 0;
  double m = 0;
  double E0 = 0;
  double E = 0;
  double E2 = 0;
  double m2 = 0;
  double m4 = 0;
  double chi = 0;
  double cb = 0;
  double u = 0;
  int n_measures = 5;

  // To compute Boltzmann-Gibbs factors
  auto h = nc::NdArray<double>(1,5);
  auto one = nc::ones<double>(1,5);

  // Temperature from high to low
  auto kBTJ = nc::fliplr(nc::linspace(T.min, T.max, T.step));

  // Set of possible energies in nbrhood
  auto E_pos = nc::linspace<double>(-4.0, 4.0, 5.0);

  // Order parameter lambda
  auto order = [](auto S, auto N){
    return nc::abs(nc::sum<int>(S).item())/((double)SZ);
  };

  // Init random matrix
  auto S = nc::random::choice<int>({-1,1},SZ).reshape(N,N);

  // Measurement matrix
  auto measurement = nc::NdArray<double>(kBTJ.shape().cols,n_measures);

  for (size_t T = 0; T < kBTJ.shape().cols; T++) {
    // Populate G-B factors
    h = nc::minimum(one, nc::exp(-2.0*E_pos/kBTJ[T]));

    // Thermalization
    for (size_t i = 0; i < therm*SZ; i++) {
      metropolis(S, N, h);
    }

    // Reset measurements
    m = 0; m2 = 0; m4 = 0; E = 0; E2 = 0;

    // Simulation and measurements
    for (size_t i = 0; i < M; i++) {
      for (size_t j = 0; j < MC*SZ; j++) {
        metropolis(S, N, h);
      }

      // Order param
      m0 = order(S,N);
      // Energy
      E0 = H(S,N);

      // Accumulate
      m += m0;
      m2 += m0*m0;
      m4 += m0*m0*m0*m0;
      E += E0;
      E2 += E0*E0;
    }
    // Average over M samples
    m = m/M;
    m2 = m2/M;
    m4 = m4/M;
    E = E/M;
    E2 = E2/M;

    // Add to measurements
    measurement(T,TEMP) = kBTJ[T];
    measurement(T,ORDER) = m;
    measurement(T,CHI) = SZ*(m2-m*m)/(kBTJ[T]);
    measurement(T,CB) = (E2-E*E)/(kBTJ[T]*kBTJ[T]);
    measurement(T,U) = 1-m4/(3*(m2*m2));

  }
  // Print columns
  tofile(measurement, file);
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

// Write columns to file
void tofile(nc::NdArray<double>& measurement, std::string file) {
  std::ofstream out;
  out.open("/home/auan/1FA573/CXX/Data/"+file);
  // Header
  out << "temp " << "order " << "chi " << "cb " << "u" << std::endl;
  for (size_t i = 0; i < measurement.shape().rows; i++) {
    for (size_t j = 0; j < measurement.shape().cols; j++) {
      out << measurement(i,j) << ' ';
    }
    out << std::endl;
  }
}

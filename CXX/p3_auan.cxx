#include <NumCpp.hpp>
#include <iostream>

// Indexing, normalizing and averages
#define SITE s[0],s[1]
#define SZ N*N
#define n1 (s[0]+1)%N, s[1]
#define n2 s[0], (s[1]+1)%N
#define n3 (s[0]-1)%N, s[1]
#define n4 s[0], (s[1]-1)%N
#define idx(dH) (dH+4)/2
#define avgnorm(x) x/((double)M*N*N)
#define norm(x) x/((double)N*N)

// Metropolis-Hastings
void metropolis(nc::NdArray<int>& S, int& N, nc::NdArray<double>& h);

// Heat-bath
void heatbath(nc::NdArray<int>& S, int& N, nc::NdArray<double>& g);

// Hamiltonian
double H(nc::NdArray<int>& S, int& N);

// Print data columns
void tofile(nc::NdArray<double>& measurements, std::string file);

int main (int argc, char *argv[]) {
  if (argc < 3 || argc >= 4) {
    std::cout << "USAGE: <lattice dim> <file>" << std::endl;
    return EXIT_FAILURE;
  }

  // Input args
  int N = atoi(argv[1]);
  std::string file = argv[2];

  // Data columns
  enum DATA_COLS {TEMP, ORDER, CHI, CB, U};

  // Global constants
  const int M = 200000; 
  const int therm = 1000;
  const int MC = 1;

  // Temperature steps
  struct {
    const double max = 5;
    const double min = 0.1;
    const double step = 80;
  } T;

  nc::random::seed(1337);

  // Measurements
  double m0 = 0,
         m = 0,
         E0 = 0,
         E = 0,
         E2 = 0,
         m2 = 0,
         m4 = 0,
         chi = 0,
         cb = 0,
         u = 0;
  int n_meas = 5;

  // To compute Boltzmann-Gibbs factors
  auto h = nc::NdArray<double>(1,5);
  auto g = nc::NdArray<double>(1,5);
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
  auto measurements = nc::NdArray<double>(kBTJ.shape().cols,n_meas);

  for (size_t T = 0; T < kBTJ.shape().cols; T++) {
    // Populate look-up arrays
    h = nc::minimum(one, nc::exp(-2.0*E_pos/kBTJ[T]));
    // g = 1.0/(1.0 + nc::exp(-2.0*E_pos/kBTJ[T]));

    // Thermalization
    for (size_t i = 0; i < therm*SZ; i++) {
      metropolis(S, N, h);
      // heatbath(S, N, g);
    }

    // Reset measurements
    m = m2 = m4 = E = E2 = 0;

    // Simulation and measurements
    for (size_t i = 0; i < M; i++) {
      for (size_t j = 0; j < MC*SZ; j++) {
        metropolis(S, N, h);
        // heatbath(S, N, g);
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
    measurements(T,TEMP) = kBTJ[T];
    measurements(T,ORDER) = m;
    measurements(T,CHI) = (m2-m*m)/(kBTJ[T]);
    measurements(T,CB) = (E2-E*E)/(kBTJ[T]*kBTJ[T]);
    measurements(T,U) = 1-m4/(3*(m2*m2));

  }
  // Print columns
  tofile(measurements, file);
  return 0;
}

// Metropolis-Hastings
void metropolis(nc::NdArray<int>& S, int& N, nc::NdArray<double>& h) {
  nc::NdArray<int> s = nc::random::randInt({1,2},N);
  double dH = S(SITE)*(S(n1)+S(n2)+S(n3)+S(n4));
  if (dH <= 0 || nc::random::rand<double>() < h[idx(dH)]) 
    S(SITE) = -S(SITE);
}

// Heat-bath
void heatbath(nc::NdArray<int>& S, int& N, nc::NdArray<double>& g) {
  nc::NdArray<int> s = nc::random::randInt({1,2},N);
  double sj = S(n1)+S(n2)+S(n3)+S(n4);
  S(SITE) = (nc::random::rand<double>() < g[idx(sj)]) ? 1: -1;
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
void tofile(nc::NdArray<double>& measurements, std::string file) {
  std::ofstream out;
  out.open("/home/auan/1FA573/CXX/Data/"+file);
  // Header
  out << "temp " << "order " << "chi " << "cb " << "u" << std::endl;
  // Data
  for (size_t i = 0; i < measurements.shape().rows; i++) {
    for (size_t j = 0; j < measurements.shape().cols; j++) {
      out << measurements(i,j) << ' ';
    }
    out << std::endl;
  }
}

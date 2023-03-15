#include <NumCpp.hpp>
#include <NumCpp/Core/Constants.hpp>
#include <NumCpp/Core/Types.hpp>
#include <NumCpp/Functions/append.hpp>
#include <NumCpp/Functions/arange.hpp>
#include <NumCpp/Functions/divide.hpp>
#include <NumCpp/Functions/exp.hpp>
#include <NumCpp/Functions/min.hpp>
#include <NumCpp/Functions/minimum.hpp>
#include <NumCpp/Functions/ones.hpp>
#include <NumCpp/Functions/sum.hpp>
#include <NumCpp/NdArray/NdArrayCore.hpp>
#include <NumCpp/Random/choice.hpp>
#include <NumCpp/Random/randInt.hpp>
#include <iostream>
#include <algorithm>
#include <iterator>
#include <ostream>
int main (int argc, char *argv[])
{
  struct
  {
    double max = 4;
    double min = 0.1;
    double step = 10;
    
  } T;

  nc::random::seed(666);
  const int N = 8;
  // Init matrix
  auto S = nc::random::choice<int>({-1,1},N*N).reshape(N,N);
  auto lattices = nc::NdArray<int>{8,16,32};
  auto h = nc::NdArray<double>{1,5};
  auto ones = nc::NdArray<double>{1.0,1.0,1.0,1.0,1.0};
  auto kBTJ = nc::linspace(T.min, T.max, T.step);
  auto E_pos = nc::linspace<double>(-4.0, 4.0, 5.0);


  std::cout << -2.0*nc::mean(S) << std::endl;
  std::cout << S << std::endl;

  for (auto& T : kBTJ) {
    auto s1 = nc::random::randInt(N);
    auto s2 = nc::random::randInt(N);

    h = nc::minimum(ones, nc::exp(-2.0*E_pos*T));
    std::cout << h << std::endl;
  }


  return 0;
}

#include <Utils.d/MathUtils.h>

long Factorial(int n)
{
  long f = 1;
  for(int i=2; i<=n; ++i) f *= i;
  return f;
}

long double DFactorial(int n)
{
  long double f = 1.0;
  for(int i=2; i<=n; ++i) f *= ((long double) i);
  return f;
}

int Combination(int n, int r)
{
  return Factorial(n)/(Factorial(r)*Factorial(n-r));
}

double DCombination(int n, int r)
{
  return DFactorial(n)/(DFactorial(r)*DFactorial(n-r));
}


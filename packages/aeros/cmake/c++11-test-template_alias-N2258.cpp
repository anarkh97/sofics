template <typename A, typename B>
class X { };
 
template <typename F>
using Xint = X<F, int>;

int main()
{
   Xint<float> x;
   return 0;
}

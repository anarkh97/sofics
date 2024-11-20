#include <algorithm>
#include <vector>
#include <iostream>

int main()
{
    std::vector<int> v;
    v.push_back(7);

    struct DivisibleBy
    {
        const int d;
        DivisibleBy(int n) : d(n) {}
        bool operator()(int n) const { return n % d == 0; }
    };
 
    if (std::all_of(v.begin(), v.end(), DivisibleBy(7))) {
        std::cout << "All numbers are divisible by 7\n";
    }

    return 0;
}

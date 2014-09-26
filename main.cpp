#include <iostream>
#include <random>

#include "quadratic_fit.h"

double qfunc(double x)
{
    const double a = 1.23;
    const double b = -9.87;
    const double c = 1e-2;

    return a * std::pow(x, 2) + b * x + c;
}

int main()
{
    std::minstd_rand gen;
    std::uniform_real_distribution<> dis(-1, 1);
    quadratic_fit<double> qfit(8);

    for (std::size_t i = 0; i < 8; ++i) {
        double x = dis(gen);
        double y = qfunc(x);
        qfit.add(x, y);

        std::cout << "Point " << i << ": (" << x << ", " << y << ")\n";
    }

    auto c = qfit.compute();

    std::cout << "a = " << c[0] << std::endl;
    std::cout << "b = " << c[1] << std::endl;
    std::cout << "c = " << c[2] << std::endl;

    return 0;
}

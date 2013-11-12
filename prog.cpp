#include <cmath>
#include <iostream>

bool VERBOSE_BISECTION = false;
bool VERBOSE_CURVE = false;

const double e1 = 1.0;
const double m1 = 1.0;
const double l1 = 3.5;

const double e2 = 5.0;
const double m2 = 1.0;
const double l2 = 1.5;

// поперечный размер, см
const double b = 3.0;

// скорость света, см/с
const double c = 3e10;
const double pi = 4 * std::atan(1);

double e_condition(double u1, double u2)
{
    return u1 * std::tan(u1 * l1) / e1 + u2 * std::tan(u2 * l2) / e2;
}

double m_condition(double u1, double u2)
{
    return m1 * std::tan(u1 * l1) / u1 + m2 * std::tan(u2 * l2) / u2;
}

template <typename F>
double bisection(F f, double left, double right, double precision=5e-3)
{
    // Метод бисекции поиска корня (трансцендентного) уравнения
    double center = (left + right) / 2;
    if (f(left) * f(right) > 0)
    {
        if (VERBOSE_BISECTION) {std::cout << "no solution" << std::endl;}
        return 0;
    }
    if ((right - left) < precision)
    {
        if (VERBOSE_BISECTION) {std::cout << "accurately" << std::endl;}
        return center;
    }
    if (f(center) * f(left) < 0)
    {
        if (VERBOSE_BISECTION) {std::cout << "left" << std::endl;}
        return bisection(f, left, center);
    }
    if (VERBOSE_BISECTION) {std::cout << "right" << std::endl;}
    return bisection(f, center, right);
}

template <typename F>
std::vector<std::vector<double>> curve(F condition, int n1, int n2,
        double precision=5e-3)
{
    // границы области, в которой ищется решение
    double left = n1 * pi / 2.0 / l1;
    double right = (n1 + 1) * pi / 2.0 / l1;
    double bottom = n2 * pi / 2.0 / l2;
    double top = (n2 + 1) * pi / 2.0 / l2;

    const int count = (right - left) / precision - 2;

    std::vector<double> U1;
    std::vector<double> U2;
    std::vector<double>::iterator i1;
    double u1, u2;

    for (int i = 0; i < count; i++) {
        U1.push_back(left + (i + 1) * precision);
    }

    i1 = U1.begin();
    int index;
    while (i1 != U1.end()) {
        u2 = bisection([i1, condition](double u2) -> double
                {
                    return condition(*i1, u2);
                },
                bottom + precision, top - precision);
        if (u2)
        {
            if (VERBOSE_CURVE) {std::cout << *i1 << ", " << u2 << std::endl;}
            U2.push_back(u2);
            ++i1;
        }
        else
        {
            index = i1 - U1.begin();
            U1.erase(i1);
            std::vector<double>(U1).swap(U1);
            i1 = U1.begin() + index;
        }
    }

    std::vector<std::vector<double>>result = {U1, U2};
    return result;
}


int main()
{
    auto x = curve(e_condition, 0, 1);
    std::cout << x[0].size() << ", " << x[1].size() << std::endl;
    return 0;
}

#include <mgl2/mgl.h>
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


std::vector<std::vector<double>> curve(std::string condition, int n,
        double precision=5e-3)
{
    int m = 2 * n - 1;
    int n1, n2;
    double right = (m + 1) * pi / 2.0 / l1;

    std::vector<double> U1;
    std::vector<double> U2;
    double u1, u2;

    while (u1 < right)
    {
        n1 = (int) (2.0 * l1 * u1 / pi);
        n2 = m - n1;
        double bottom = n2 * pi / 2.0 / l2;
        double top = (n2 + 1) * pi / 2.0 / l2;
        u2 = bisection([u1, condition](double u2)
                -> double
                {
                    if (condition == "e") { return e_condition(u1, u2);};
                    if (condition == "m") { return m_condition(u1, u2);};
                    return u1 + u2; // for amazing results
                },
                bottom + precision, top - precision);
        if (u2)
        {
            if (VERBOSE_CURVE) {std::cout << u1 << ", " << u2 << std::endl;}
            U1.push_back(u1);
            U2.push_back(u2);
        }
        u1 += precision;
    }
    std::vector<std::vector<double>>result = {U1, U2};
    return result;
}

void plot(mglGraph *gr, std::string condition, int number)
{
    auto X = curve(condition, number);
    mglData x;
    mglData y;
    x.Set(X[0]);  // convert to internal format
    y.Set(X[1]);  // convert to internal format
    gr->Plot(x,y,"k");   // plot it
}

void plot_transversal_wavenumbers(std::string condition, const int number)
{
    mglGraph gr; // create canvas
    gr.SetRanges(0, number * pi / l1 , 0, number * pi / l2);
    for (int i = 1; i <= number; i++) {
        plot(&gr, condition, i);
    }
    gr.Axis();
    gr.Label('x', "u_1, cm^{-1}", 0);
    gr.Label('y', "u_2, cm^{-1}", 0);
    gr.Grid("k");

    auto name = condition + ".eps";
    gr.WriteFrame(name.c_str()); // save it
}

int main(int argc, const char *argv[])
{
    plot_transversal_wavenumbers("e", 5);
    plot_transversal_wavenumbers("m", 5);
    return 0;
}


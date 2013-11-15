#include <mgl2/mgl.h>
#include <cmath>
#include <iostream>

/******************************************************************************
 *                                 Для отладки
 *****************************************************************************/
bool VERBOSE_BISECTION = false;
bool VERBOSE_curves = false;

/******************************************************************************
 *                          Свойства и размеры сред
 *****************************************************************************/
const double e1 = 1.0;
const double m1 = 1.0;
const double l1 = 3.5;

const double e2 = 5.0;
const double m2 = 1.0;
const double l2 = 1.5;

// поперечный размер, см
const double b = 3.0;

/******************************************************************************
 *                          Универсальные постоянные
 *****************************************************************************/
// скорость света, см/с
const double c = 3e10;
// число π
const double pi = std::acos(-1);


/******************************************************************************
 *                           Дисперсионные уравнения
 *****************************************************************************/
double e_condition(double u1, double u2)
{
    return u1 * std::tan(u1 * l1) / e1 + u2 * std::tan(u2 * l2) / e2;
}

double m_condition(double u1, double u2)
{
    return m1 * std::tan(u1 * l1) / u1 + m2 * std::tan(u2 * l2) / u2;
}

double transversal_wavenumbers_relation(double omega, double u1, double u2)
{
    return pow(omega/c, 2)*(e1 * m1 - e2 * m2) - (pow(u1,2) - pow(u2,2));
}

/*****************************************************************************/

// Метод бисекции поиска корня (трансцендентного) уравнения
template <typename F>
double bisection(F f, double left, double right, double precision=5e-3)
{
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

// нахождение поперечного волнового числа во 2 области для известного в 1
// области
//
// condition ∈ {"e", "m"}
// number — номер кривой, на которой искать точку
double second_from_first(std::string condition, double u1, int number,
        double precision=5e-3)
{
    double n1 = (int) (2.0 * l1 * u1 / pi);
    double n2 = 2 * number - 1 - n1;
    double bottom = n2 * pi / 2.0 / l2;
    double top = (n2 + 1) * pi / 2.0 / l2;
    double u2 = bisection([u1, condition](double u2)
            -> double
            {
                if (condition == "e") { return e_condition(u1, u2);};
                if (condition == "m") { return m_condition(u1, u2);};
                return u1 + u2; // for amazing results
            },
            bottom + precision, top - precision);
    return u2;
}

// определение координат для кривых, определяемых 1 и 2 дисперсионными
// уравнениями
//
// condition ∈ {"e", "m"}
// number — номер кривой
std::vector<std::vector<double>> curves(std::string condition, int number,
        double precision=5e-3)
{
    double right = number * pi / l1;

    std::vector<double> U1;
    std::vector<double> U2;
    double u1, u2;

    while (u1 < right)
    {
        u2 = second_from_first(condition, u1, number);
        if (u2)
        {
            if (VERBOSE_curves) {std::cout << u1 << ", " << u2 << std::endl;}
            U1.push_back(u1);
            U2.push_back(u2);
        }
        u1 += precision;
    }
    std::vector<std::vector<double>>result = {U1, U2};
    return result;
}

// определение поперечных волновых чисел для заданной частоты при помощи
// 3-го дисперсионного уравнения
std::pair<double, double> transversal_wavenumbers(std::string condition,
        int number, double omega, double precision=1e-5)
{
    int m = 2 * number - 1;
    double left = 0;
    double right = (m + 1) * pi / 2.0 / l1;
    double u1, u2;
    u1 = bisection([condition, omega, number](double u1)
            -> double
            {
                double u2 = second_from_first(condition, u1, number);
                return transversal_wavenumbers_relation(omega, u1, u2);
            },
            left + precision, right - precision);

    u2 = second_from_first(condition, u1, number);
    std::pair<double, double>result = {u1, u2};
    return result;
}

// продольное волновое число
std::vector<std::vector<double>> longitudinal_wavenumber(std::string condition,
        int n, int k)
{
    std::vector<double> H, O;
    double omega = 2e10;
    double sqr_h, u1, u2;
    std::pair<double,double> tw;
    while (omega < 6e10)
    {
        tw = transversal_wavenumbers(condition, k, omega, 1e-7);
        u1 = tw.first; u2 = tw.second;
        sqr_h = pow(omega / c, 2) * (e1 * m1) - u1 * u1 -
        pow((pi * n / b), 2);
        if (u1 && u2 && (sqr_h > 0))
        {
            O.push_back(omega);
            H.push_back(sqrt(sqr_h));
        }
        omega +=1e9;
    }
    std::vector<std::vector<double>> result;
    result = {O, H};
    return result;
}

void plot(mglGraph *gr, std::string condition,
        std::vector<double> x, std::vector<double> y)
{
    mglData X;
    mglData Y;
    X.Set(x);  // convert to internal format
    Y.Set(y);  // convert to internal format
    gr->Plot(X, Y, "k");   // plot it
}

void plot_curves(std::string condition, const int number)
{
    mglGraph gr; // create canvas
    gr.SetRanges(0, number * pi / l1 , 0, number * pi / l2);
    for (int i = 1; i <= number; i++) {
        auto data = curves(condition, i);
        plot(&gr, condition, data[0], data[1]);
    }
    gr.Axis();
    gr.Label('x', "u_1, cm^{-1}", 0);
    gr.Label('y', "u_2, cm^{-1}", 0);
    gr.Grid("k");

    auto name = condition + ".eps";
    gr.WriteFrame(name.c_str()); // save it
}

void plot_dispersion_relation(std::string condition, int n, std::vector<int>N)
{
    mglGraph gr; // create canvas
    gr.SetRanges(0, 6e10 , 0, 2);
    for (int i=0; i<N.size(); ++i)
    {
        auto data = longitudinal_wavenumber(condition, n, N[i]);
        plot(&gr, condition, data[0], data[1]);
    }
    gr.Axis();
    gr.Label('x', "\\omega, rad/s}", 0);
    gr.Label('y', "h, cm^{-1}", 0);
    gr.Grid("k");

    auto name = condition + "_h.eps";
    gr.WriteFrame(name.c_str()); // save it
}


int main(int argc, const char *argv[])
{
    plot_curves("e", 5);
    plot_curves("m", 5);
    plot_dispersion_relation("e", 1, {1,2,3});
    plot_dispersion_relation("m", 1, {1,2,3});
    return 0;
}


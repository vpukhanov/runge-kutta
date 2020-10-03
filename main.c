#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX_STEPS 100

int runge_kutta(double, double, double, double, double *);
double f(double, double);
double y_exp(double);

int main(void) {
    double series[MAX_STEPS];
    double x0 = 0;
    double y0 = 1;
    double xn = 2 * M_PI;
    double h = M_PI / 8;
    int steps, i;

    printf("Solving y' = sin(x), y(0) = 0 using Runge-Kutta method\n");
    printf("(with x0 = %.5f, y0 = %.5f, h = %.5f, xn = %.5f)\n\n", x0, y0, h, xn);

    steps = runge_kutta(x0, y0, h, xn, series);

    printf("x\t\t y(x)\t\t\t expected\n");
    for (i = 0; i < steps; i++) {
        double x = x0 + h * i;
        printf("%7.5f\t\t%8.5f\t\t%8.5f\n", x, series[i], y_exp(x));
    }

    return EXIT_SUCCESS;
}

int runge_kutta(double x0, double y0, double h, double xn, double *result) {
    int steps = 0;
    double x = x0;
    double y = y0;

    do {
        double k1 = f(x, y);
        double k2 = f(x + h / 2, y + h * k1 / 2);
        double k3 = f(x + h / 2, y + h * k2 / 2);
        double k4 = f(x + h, y + h * k3);

        result[steps++] = y;

        y += h * (k1 + 2 * k2 + 2 * k3 + k4) / 6;
        x += h;
    } while (x < xn);

    return steps;
}

double f(double x, double y) {
    /* y is unused in this example */
    return -sin(x);
}

double y_exp(double x) {
    /* in this example y = cos(x) */
    return cos(x);
}

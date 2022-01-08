
#include <cmath>
#include "gnuplot-iostream.hpp"
#include <iostream>
#include <fstream>
#include <vector>



double EE (unsigned int tt, double dt, double y_old, const std::function<double(const double &, const double &)> &f);

double predictor_corrector (unsigned int tt, double dt, double y_old, const std::function<double(const double &, const double &)> &f);

double Heun (unsigned int tt, double dt, double y_old, const std::function<double(const double &, const double &)> &f);

double IE (unsigned int tt, double dt, double y_old, const std::function<double(const double &, const double &)> &f);

double compute_J(double dt, unsigned int tt, double y);

double reaction_rate (unsigned int tt, double dt, double y_old, const std::function<double(const double &, const double &)> &f);

void output_results(const std::vector<unsigned int>& time, const std::vector<double>& y_EE, const std::vector<double>& y_IE, const std::vector<double>& y_PC, const std::vector<double>& y_H, const std::vector<double>& y_ex);

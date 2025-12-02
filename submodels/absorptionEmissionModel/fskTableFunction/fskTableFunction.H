/* Header file with data structures and functions required for utilizing
 * FSK Table v3 within OpenFOAM.
 * 
 * Copyright @ Nehal Jajal, The Ohio State University
 */

#ifndef fskTableFunction_H
#define fskTableFunction_H

#include <iostream>
#include <vector>
#include <string>

struct GasMixInfo
{
    double P,T,fv,xCO2,xH2O,xCO,xNH3,xCH4;
};

void loadIntoMem(std::vector<float>& kq, std::string fpath);

void loadIntoMemV2(std::vector<float>& kq, std::string fpath, const double Tref);

void loadIntoMemFuel(std::vector<float>& kq, std::string fpath, const double Tref);

void loadIntoMemProduct(std::vector<float>& kq, std::string fpath, const double Tref);

void get_ka(const GasMixInfo& Mix_Info, const double Tref, const std::vector<float>& kq, const std::vector<double>& g, const std::vector<double>& weight, const bool isV2, std::vector<double>& k, std::vector<double>& a);

void get_ka_mix(const GasMixInfo& Mix_Info, const double Tref, const std::vector<float>& kq, const std::vector<float>& kq2, const std::vector<double>& g, const std::vector<double>& weight, std::vector<double>& k, std::vector<double>& a, const int mix_scheme);

void quadgen2(const bool cheb2, std::vector<double>& g, std::vector<double>& w, const int points, const double alpha);

std::vector<double> get_k(const GasMixInfo& Mix_Info, double Tref, const std::vector<float>& kq);

std::vector<double> get_k_fuel(const GasMixInfo& Mix_Info, double Tref, const std::vector<float>& kq);

std::vector<double> get_k_product(const GasMixInfo& Mix_Info, double Tref, const std::vector<float>& kq);

std::vector<double> get_k_V2(const GasMixInfo& Mix_Info, double Tref, const std::vector<float>& kq);

std::vector<double> afun(const std::vector<double> grr,const std::vector<double> krr,const std::vector<double> grc,const std::vector<double> krc,const std::vector<double> g, const std::vector<double> w);

std::vector<double> linearInterpMono(const std::vector<double> xx,const std::vector<double> yy,const std::vector<double> xi);

void coef_setup(const GasMixInfo& Mix_Info, const double Tref, std::vector<std::vector<int>>& grid, std::vector<double>& coef);

void coef_setupV2(const GasMixInfo& Mix_Info, const double Tref, std::vector<std::vector<int>>& grid, std::vector<double>& coef);

std::vector<int> detmLP(std::vector<double> arr); // convert to vector

int detmID(std::vector<double> arr, double val);

std::vector<double> get_k_mix(const GasMixInfo& Mix_Info, double Tref, const std::vector<double>& g, const std::vector<float>& kq, const std::vector<float>& kq2, const std::vector<double>& weight, const int mix_scheme);

std::vector<double> multiplication(
    const std::vector<double>& k1,
    const std::vector<double>& g1,
    const std::vector<double>& k2,
    const std::vector<double>& g2,
    const std::vector<double>& k);

std::vector<double> superposition(
    const std::vector<double>& k1,
    const std::vector<double>& g1,
    const std::vector<double>& k2,
    const std::vector<double>& g2,
    const std::vector<double>& k);

std::vector<double> modestRiazzi(
    const std::vector<double>& k1,
    const std::vector<double>& g1,
    const std::vector<double>& k2,
    const std::vector<double>& g2,
    const std::vector<double>& k,
    const std::vector<double>& wq);

void bound(const std::vector<double>& k, std::vector<double>& g, double kmin, double kmax);

std::vector<double> kPowerLaw(double kmin, double kmax, int n, double pwr);

#endif

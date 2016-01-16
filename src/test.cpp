#include <iostream>
#include <armadillo>
#include <complex>
#include <cmath>

using namespace std;
using namespace arma;

const double PlankConstant = 6.62606876e-34;
const double hbar = PlankConstant/M_PI;
const double T = 1e3;


mat Hamilton(double s)
{
    mat res = {{1+s,s-1},{s-1,1-s}};
    return res;
}

cx_vec next_timestep(cx_vec state,double time,double dt=1e-2)
{
    return expmat(complex<double>(0,1/hbar*dt)*Hamilton(time/T))*state;
}

int main()
{
    double dt = 1e-2;
    cx_vec state = {1,0};
    for(int i=0;i<1e3;i++)
    {
        state = next_timestep(state,i*dt,dt);
        cout<<conj(state.t())*Hamilton(i*dt/T)*state<<endl;
    }
}
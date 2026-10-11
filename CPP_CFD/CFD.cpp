#include <chrono>
#include <iostream>
#include <limits>
#include <string>
#include <vector>
#include <numbers>
#include <algorithm>
#include "matplotlibcpp.h"
namespace plt = matplotlibcpp;


double pi = std::numbers::pi;
template <typename T>
class Helper{
    public:
        
        
        T max_element_2D(std::vector<std::vector<T>> &mat){
            T abs_max = std::numeric_limits<T>::lowest();
            for(auto& row : mat){
                if(!row.empty()){
                    T row_max = *std::max_element(
                        row.begin(), row.end());
                    abs_max = std::max(abs_max, row_max);
                }
            }
        return abs_max; 
        }

        static std::pair<std::vector<std::vector<T>>, 
        std::vector<std::vector<T>>>
        CylidericalMeshGrid(int Nr, int Nphi, double ri, double ro){
            double dr = (ro - ri)/Nr;
            std::vector<T> r_vec(Nr + 1);
            std::vector<T> theta_vec(Nphi);
            for(int i = 0; i <= Nr; ++i){
                r_vec[i] = ri + i*dr;
            }
            for(int j = 0; j < Nphi; ++j){
                theta_vec[j] = 2.0*pi*j/(Nphi - 1);
            }

            std::vector<std::vector<T>> r(Nr + 1, 
            std::vector<T>(Nphi));
            std::vector<std::vector<T>> theta(Nr + 1,
            std::vector<T>(Nphi));
            for(int i = 0; i <= Nr; ++i){
                for(int j = 0; j < Nphi ; ++j){
                    r[i][j] = r_vec[i];
                    theta[i][j] = theta_vec[j];
                }
            }

            std::vector<std::vector<T>> x(Nr + 1, 
            std::vector<T>(Nphi));
            std::vector<std::vector<T>> y(Nr + 1,
            std::vector<T>(Nphi));
            
            for(int i = 0; i <= Nr; i++){
                for(int j = 0; j < Nphi; j++){
                    x[i][j] = r[i][j]*std::cos(theta[i][j]);
                    y[i][j] = r[i][j]*std::sin(theta[i][j]);
                }
            }
            return {x, y};
    }
};


class CylinderialTemp:public Helper<double>{
    public:
        // --------------- Physics Parametrs -------------
        double a = 1.16e-5;
        double ri = 0.02, ro = 0.05;
        double k = 45.0, Ta = 20.0, Tw = 5.0, 
        ha = 50.0, hw = 250;
        // ---------- Mathematical quantaties -------------

    std::pair<std::vector<std::vector<double>>, double>
    Explicit(double Tin, 
        int Nr, int Nphi, int tstps, double dt){
        double dr = (ro - ri)/Nr;
        double dphi = 2*pi/Nphi;
        double dr_ri = dr/ri, dr_dphi = dr/dphi;
        double d = a * dt/(dr*dr);
        double Time_taken = 0;
        std::vector<std::vector<double>> T(Nr + 1, 
            std::vector<double>(Nphi, Tin));
        std::vector<std::vector<double>> Tprev = T;
        
        for(int stp = 0; stp < tstps; stp++){
            std::swap(T, Tprev);
            # pragma omp parallel for num_threads(8)
            for(int i = 1; i < Nr; i++){
                double r = ri + i*dr;
                // # pragma omp parallel for
                for(int j = 0; j < Nphi; j++){
                    int jf = (j + 1)%(Nphi);
                    int jp = (j - 1 + Nphi)%(Nphi);
                    T[i][j] = Tprev[i][j] +
                            + d*((Tprev[i+1][j] - 2*Tprev[i][j] + Tprev[i - 1][j])
                            + (0.5*dr/r)*(Tprev[i + 1][j] - Tprev[i - 1][j])
                            + (dr_dphi*dr_dphi/(r*r))*(Tprev[i][jf] - 2*Tprev[i][j] + Tprev[i][jp]));
                }
            } 
            T[0] = T[1];
            // outer boundary condition
            # pragma omp parallel for num_threads(8)
            for(int j = 0; j < Nphi; j++){
                int i = Nr;
                double r = ri + i*dr;
                int jf = (j + 1)%Nphi;
                int jp = (j - 1 + Nphi)%Nphi;
                if(j <= Nphi/2){
                    double Td = -(2*dr*ha/k)*(Tprev[i][j] - Ta) + Tprev[i-1][j];
                    T[i][j] = Tprev[i][j] +
                            + d*((Td - 2*Tprev[i][j] + Tprev[i - 1][j])
                            + (0.5*dr/r)*(Td - Tprev[i - 1][j])
                            + (dr_dphi*dr_dphi/(r*r))*(Tprev[i][jf] - 2*Tprev[i][j] + Tprev[i][jp]));
                }
                else{
                    double Td = -(2*dr*hw/k)*(Tprev[i][j] - Tw) + Tprev[i-1][j];
                    T[i][j] = Tprev[i][j] +
                            + d*((Td - 2*Tprev[i][j] + Tprev[i - 1][j])
                            + (0.5*dr/r)*(Td - Tprev[i - 1][j])
                            + (dr_dphi*dr_dphi/(r*r))*(Tprev[i][jf] - 2*Tprev[i][j] + Tprev[i][jp]));
                }

            }
            if(stp%100 == 0){
                std::cout<<"Time step: " << stp << " seconds "<< "Temperature " <<
                max_element_2D(T) <<std::endl;
                
            }
            if(max_element_2D(T) < 100 ){
                Time_taken = stp * dt;
                std::cout << "Temperature reached 100" <<" "<< max_element_2D(T) <<std::endl;
                break;
            }
            
        }
        return {T, Time_taken};
    }
};



int main()
{
    double Tin = 200;
    int Nr = 30;
    int Nphi = 50;
    int tstps = 10000000;
    double dt = 1/100.0;
    double ri = 0.02, ro = 0.05;
    std::cout <<"Starting solver"<<std::endl;
    auto [x, y] = Helper<double>::CylidericalMeshGrid(Nr, Nphi, ri, ro);
    auto start = std::chrono::high_resolution_clock::now();
    CylinderialTemp obj = CylinderialTemp();
    auto [T, Time_taken] = obj.Explicit(Tin , 
        Nr ,Nphi ,tstps, dt);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    int colorinterpolation = 50;
    plt::figure_size(700, 500);
    PyObject* contour = plt::contourf(
        x, y, T, colorinterpolation, "jet"
    );
    // plt::colorbar();
    plt::colorbar(contour);

Py_DECREF(contour);
    plt::title(
        "Temperatur profile with grid (Nr, Nphi) = (" + 
        std::to_string(Nr) + "," + 
        std::to_string(Nphi) + 
        ") at Time = " +
        std::to_string(Time_taken)+ "\n"
        "Time taken by solver =" +
        std::to_string(elapsed.count())
    );
    plt::xlabel("x (m)");
    plt::ylabel("y (m)");
    plt::axis("equal");

    plt::plot(
        std::vector<double> {-ro, ro},
        std::vector<double> {0.0, 0.0},
        "k--"
    );
    plt::text(
        0.0,
        ro/2.0,
        "Air (20 C)"
    );
    plt::text(
        0.0,
        -ro/2.0,
        "Water (5 C)"
    );
    plt::show();
    return 0;

}
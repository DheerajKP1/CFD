#include <chrono>
#include <iostream>
#include <limits>
#include <string>
#include <vector>
#include <numbers>
#include <algorithm>
#include "matplotlibcpp.h"
#include "mpi.h"

namespace plt = matplotlibcpp;
int LEFT = 100;
int RIGHT = 101;
int Npes = 3;
int NR = 30;
int NPHI = 50;

double pi = std::numbers::pi;
template <typename T>
class Helper{
    public:
        T max_element_2D(std::vector<std::vector<T>> &mat){
            T abs_max = std::numeric_limits<T>::lowest();
            for(int i = 1; i < NR/Npes + 1; i ++){
                if(!mat[i].empty()){
                    double row_max = *std::max_element(mat[i].begin(), mat[i].end());
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
        int Nr, int Nphi, int tstps, double dt, int mype, 
        int npes, MPI_Status status){  
        double dr = (ro - ri)/NR;
        double dphi = 2*pi/Nphi;
        double dr_ri = dr/ri, dr_dphi = dr/dphi;
        double d = a * dt/(dr*dr);
        double offset = Nr*mype;
        double Time_taken = 0;
        std::vector<std::vector<double>> T(Nr + 2, 
            std::vector<double>(Nphi, Tin));
        std::vector<std::vector<double>> Tprev = T;
        std::cout<<"Nr"<<Nr<<"Nphi"<<Nphi<<"tstps"<<tstps<<"dt"<<dt<<"mype"<<mype<<std::endl;
        for(int stp = 0; stp < tstps; stp++){
            std::swap(T, Tprev);
            // # pragma omp parallel for num_threads(8)
            for(int i = 1; i <= Nr; i++){
                double r = ri + (offset + i - 1)*dr;
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
            if(mype < npes - 1){
                MPI_Send(&T[Nr][0], Nphi, MPI_DOUBLE, mype + 1, RIGHT, MPI_COMM_WORLD);
            }
            if(mype != 0){
                MPI_Send(&T[1][0], Nphi, MPI_DOUBLE, mype - 1, LEFT, MPI_COMM_WORLD);
            }
            if(mype != 0){
                MPI_Recv(&T[0][0], Nphi, MPI_DOUBLE, MPI_ANY_SOURCE, RIGHT, MPI_COMM_WORLD,
                &status);
            }
            if(mype != npes - 1){
                MPI_Recv(&T[Nr + 1][0], Nphi, MPI_DOUBLE, MPI_ANY_SOURCE, LEFT, MPI_COMM_WORLD,
                &status);
            }
            if(mype == 0){
                T[1] = T[2];
                // T[0] = T[2];
            }
            // outer boundary condition
            // # pragma omp parallel for num_threads(8)
            if(mype == npes - 1){
                for(int j = 0; j < Nphi; j++){
                    int i = Nr;
                    double r = ri + (offset + i - 1)*dr;
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
                //  T[Nr + 1] = T[Nr]; 
            }

            if(stp%100 == 0){
                std::cout<<"Time step: " << stp*dt << " seconds "<< "Temperature " <<
                max_element_2D(T) <<std::endl;
                
            }
            // if(max_element_2D(T) < 100 ){
            //     Time_taken = stp * dt;
            //     std::cout << "Temperature reached 100" <<std::endl;
                
            //     break;
            // }
            double local_max = max_element_2D(T);
            double global_max = 0;
            
            MPI_Allreduce(
                &local_max,
                &global_max,
                1, 
                MPI_DOUBLE,
                MPI_MAX,
                MPI_COMM_WORLD
            );
            if(global_max < 100.0){
                Time_taken = stp*dt;
                if(mype == 0){
                    std::cout << "Temperature reached 100"
                              << " at time " <<Time_taken
                              << " seconds. Global max = "
                              << global_max <<std::endl;
                }
                break;
            }
        // MPI_Barrier(MPI_COMM_WORLD); 
        }
        return {T, Time_taken};
    }    
};



int main(int argc, char* argv[])
{
    int npes;
    int mype;
    int stat;
    int niter;
    MPI_Status status;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &npes);
    MPI_Comm_rank(MPI_COMM_WORLD, &mype);


    double Tin = 200;
    int Nr = NR/Npes;
    int Nphi = NPHI;
    int tstps = 10000000;
    double dt = 1/100.0;
    double ri = 0.02, ro = 0.05;
    std::cout <<"Starting solver"<<std::endl;
    
    auto start = std::chrono::high_resolution_clock::now();
    CylinderialTemp obj = CylinderialTemp();
    auto [T, Time_taken] = obj.Explicit(Tin , 
        Nr ,Nphi ,tstps, dt, mype, npes, status);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    
    std::vector<double> local_T(Nr*Nphi);
    for(int i = 1; i <= Nr; i++){
        for(int j = 0; j < Nphi; j ++){
            local_T[(i - 1)*Nphi + j] = T[i][j];
        }
    }
    std::vector<double> global_T;
    if(mype == 0){
        global_T.resize(NR*Nphi);
    }
    MPI_Gather(
        local_T.data(),
        Nr*Nphi,
        MPI_DOUBLE,
        mype == 0?global_T.data():nullptr,
        Nr*Nphi, MPI_DOUBLE,
        0, MPI_COMM_WORLD
    );
    if(mype == 0){
        std::vector<std::vector<double>> T(NR,
        std::vector<double>(Nphi));
        for(int i = 0; i < NR; i++){
            for(int j = 0; j <Nphi; j++){
                T[i][j] = global_T[i*Nphi + j];
            }
        }

    auto [x, y] = Helper<double>::CylidericalMeshGrid(NR - 1, Nphi, ri, ro);
    std::cout<<" T sixe "<<T.size()<< " "<< T[0].size()
            <<" x " << x.size()<< " " << x[0].size()<<" y "<< y.size()<<std::endl;
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
            std::to_string(NR) + "," + 
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
    }
    MPI_Finalize();
    return 0;

}
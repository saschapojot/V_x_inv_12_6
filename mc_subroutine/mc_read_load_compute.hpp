//
// Created by polya on 7/19/24.
//

#ifndef V_INV_12_6_MC_READ_LOAD_COMPUTE_HPP
#define V_INV_12_6_MC_READ_LOAD_COMPUTE_HPP
#include "../potentialFunction/potentialFunctionPrototype.hpp"
#include <boost/filesystem.hpp>

#include <boost/math/quadrature/trapezoidal.hpp>
#include <chrono>
#include <cstdlib>
#include <cxxabi.h>
#include <fstream>
#include <initializer_list>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <random>
#include <regex>
#include <sstream>
#include <string>
#include <typeinfo>
#include <vector>

namespace fs = boost::filesystem;
const auto PI=M_PI;


class mc_computation {
public:
    mc_computation(const std::string &cppInParamsFileName) {
        std::ifstream file(cppInParamsFileName);
        if (!file.is_open()) {
            std::cerr << "Failed to open the file." << std::endl;
            std::exit(20);
        }
        std::string line;

        int paramCounter = 0;

        while (std::getline(file, line)) {
            // Check if the line is empty
            if (line.empty()) {
                continue; // Skip empty lines
            }
            std::istringstream iss(line);

            //read T
            if (paramCounter == 0) {
                iss >> T;
                if (T <= 0) {
                    std::cerr << "T must be >0" << std::endl;
                    std::exit(1);
                }//end if
                std::cout << "T=" << T << std::endl;
                this->beta = 1 / T;
                double stepForT1 = 0.1;
                double h_threshhold=0.01;
                this->h=h_threshhold;
//                this->h = stepForT1 * T > h_threshhold ? h_threshhold : stepForT1 * T;//stepSize;
                std::cout << "h=" << h << std::endl;
                this->M = std::pow(2.0 * PI, 0.5) * h * 1.001;
                std::cout<<"M="<<M<<std::endl;
                paramCounter++;
                continue;
            }//end reading T
            //read coefficients
            if(paramCounter==1){
                iss>>coefsToPotFunc;
                paramCounter++;
                continue;

            }// end reading coefficients


            //read potential function name
            if(paramCounter==2){
                iss>>potFuncName;
                paramCounter++;
                continue;
            }//end reading potential function name

            //read initial values
            if(paramCounter==3){
                std::string temp;
                if (std::getline(iss, temp, ',')){
                    LInit=std::stod(temp);
                }
                if (std::getline(iss, temp, ',')){
                    x0AInit=std::stod(temp);
                }
                if (std::getline(iss, temp, ',')){
                    x0BInit=std::stod(temp);
                }
                if (std::getline(iss, temp, ',')){
                    x1AInit=std::stod(temp);
                }
                if (std::getline(iss, temp, ',')){
                    x1BInit=std::stod(temp);
                }
                paramCounter++;
                continue;




            }//end reading initial values


            //read loopToWrite
            if(paramCounter==4){
                iss>>loopToWrite;
                paramCounter++;
                continue;
            }//end reading loopToWrite

            //read newFlushNum
            if(paramCounter==5){
                iss>>newFlushNum;
                paramCounter++;
                continue;
            }//end reading newFlushNum

            //read loopLastFile
            if(paramCounter==6){
                //if loopLastFileStr is "-1", loopLastFile uses the overflowed value
                //and loopLastFile+1 will be 0
                iss>>loopLastFile;
                paramCounter++;
                continue;
            }//end reading loopLastFile

            //read TDirRoot
            if (paramCounter==7){
                iss>>TDirRoot;
                paramCounter++;
                continue;
            }//end reading TDirRoot

            //read U_dist_dataDir
            if(paramCounter==8){
                iss>>U_dist_dataDir;
                paramCounter++;
                continue;
            }//end reading U_dist_dataDir



        }//end while
        this->potFuncPtr = createPotentialFunction(potFuncName, coefsToPotFunc);
        potFuncPtr->init();
        this->varNum = 6;
        try {
            this->U_dist_ptr= std::shared_ptr<double[]>(new double[loopToWrite * varNum],
                                                        std::default_delete<double[]>());
        }
        catch (const std::bad_alloc &e) {
            std::cerr << "Memory allocation error: " << e.what() << std::endl;
            std::exit(2);
        } catch (const std::exception &e) {
            std::cerr << "Exception: " << e.what() << std::endl;
            std::exit(2);
        }
        std::cout<<"LInit="<<LInit<<", x0AInit="<<x0AInit
                 <<", x0BInit="<<x0BInit<<", x1AInit="<<x1AInit<<", x1BInit="<<x1BInit<<std::endl;

        std::cout<<"loopToWrite="<<loopToWrite<<std::endl;
        std::cout<<"newFlushNum="<<newFlushNum<<std::endl;
        std::cout<<"loopLastFile+1="<<loopLastFile+1<<std::endl;
        std::cout<<"TDirRoot="<<TDirRoot<<std::endl;
        std::cout<<"U_dist_dataDir="<<U_dist_dataDir<<std::endl;

    }//end constructor




public:



    ///
    /// @param a
    /// @param b
    /// @return a % b
    double positive_fmod(const double& a,const  double& b) {
        double result = fmod(a, b);
        if (result < 0) {
            result += b;
        }
        return result;
    }

    ///
    /// @param LCurr
    /// @param x0ACurr
    /// @param x0BCurr
    /// @param x1ACurr
    /// @param x1BCurr
    /// @param LNext
    /// @param x0ANext
    /// @param x0BNext
    /// @param x1ANext
    /// @param x1BNext
    void mod_proposal(const double &LCurr, const double &x0ACurr, const double &x0BCurr, const double& x1ACurr,const double &x1BCurr,
                  double &LNext, double &x0ANext, double &x0BNext, double &x1ANext, double &x1BNext);

    ///
/// @param y
/// @param x center
/// @param a left end
///@param b right end
/// @return known proposal function, which is normal distribution
    double Q(const double &y, const double &x, const double &a, const double &b);


    ///
    /// @param y
    /// @param x center
    /// @param a left end
    /// @param b right end
    /// @return truncated Gaussian
    double f(const double &y, const double &x, const double &a, const double &b);

    double generate_nearby_normal(const double &x);

    ///
    /// @param x center
    /// @param a left end
    /// @param b right end
    /// @return integral
    double zVal(const double& x,const double &a, const double &b);


    ///
    /// @param x center
    /// @param a left end
    /// @param b right end
    /// @return random number from truncated Gaussian
    double reject_sampling_one_data(const double &x,const double &a, const double &b);
    ///
    /// @param y
    /// @param x center
    /// @param a left end
    /// @param b right end
    /// @return
    double integrand(const double &y, const double& x,const double &a, const double &b);
    void execute_mc(const double& L,const double &x0A, const double &x0B, const double& x1A, const double& x1B, const size_t & loopInit, const size_t & flushNum);


    static void saveArrayToCSV(const std::shared_ptr<double[]>& array, const  size_t& arraySize, const std::string& filename, const size_t& numbersPerRow) ;

    ///
    /// @param LCurr
    /// @param x0ACurr
    /// @param x0BCurr
    /// @param x1ACurr
    /// @param x1BCurr
    /// @param UCurr
    /// @param LNext
    /// @param x0ANext
    /// @param x0BNext
    /// @param x1ANext
    /// @param x1BNext
    /// @param UNext
    /// @return
    double mod_acceptanceRatio(const double &LCurr,const double &x0ACurr, const double &x0BCurr,
                           const double &x1ACurr, const double &x1BCurr,const double& UCurr,
                           const double & LNext, const double &x0ANext, const double &x0BNext,
                           const double &x1ANext, const double &x1BNext,double &UNext);

    void init_and_run();

public:
    double T;// temperature
    double beta;
    double h;// step size
    size_t loopToWrite;
    size_t newFlushNum;
    size_t loopLastFile;
    std::shared_ptr<potentialFunction> potFuncPtr;
    std::string TDirRoot;
    std::string U_dist_dataDir;
    std::shared_ptr<double[]> U_dist_ptr;
    size_t varNum;
    double LInit;
    double x0AInit;
    double x0BInit;
    double x1AInit;
    double x1BInit;
    std::string coefsToPotFunc;
    std::string potFuncName;
    double M;

};


#endif //V_INV_12_6_MC_READ_LOAD_COMPUTE_HPP

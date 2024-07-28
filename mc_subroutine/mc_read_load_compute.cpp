//
// Created by polya on 7/19/24.
//

#include "mc_read_load_compute.hpp"



void mc_computation::execute_mc(const double& L,const double &x0A, const double &x0B, const double& x1A, const double& x1B, const size_t & loopInit, const size_t & flushNum){

    double LCurr = L;
    double x0ACurr = x0A;
    double x0BCurr = x0B;
    double x1ACurr = x1A;
    double x1BCurr=x1B;
    std::cout<<"Before mc: "<<"LCurr="<<LCurr<<", x0ACurr="<<x0ACurr<<", x0BCurr="<<x0BCurr<<", x1ACurr="<<x1ACurr<<", x1BCurr="<<x1BCurr<<std::endl;
    double UCurr;// = (*potFuncPtr)(LCurr, y0Curr, z0Curr, y1Curr);
    std::random_device rd;
    std::ranlux24_base e2(rd());
    std::uniform_real_distribution<> distUnif01(0, 1);//[0,1)
    size_t loopStart = loopInit;
    for (size_t fls = 0; fls < flushNum; fls++) {
        const auto tMCStart{std::chrono::steady_clock::now()};
        for (size_t j = 0; j < loopToWrite; j++) {
            //propose a move
            double LNext;
            double x0ANext;
            double x0BNext;
            double x1ANext;
            double x1BNext;

            this->mod_proposal(LCurr,x0ACurr,x0BCurr,x1ACurr,x1BCurr,
                               LNext,x0ANext,x0BNext,x1ANext,x1BNext);
            double UNext;
            UCurr=((*potFuncPtr))(LCurr,x0ACurr,x0BCurr,x1ACurr,x1BCurr);
            double r= mod_acceptanceRatio(LCurr,x0ACurr,x0BCurr,x1ACurr,x1BCurr,UCurr,
                                          LNext,x0ANext, x0BNext,x1ANext,x1BNext,UNext);
            double u = distUnif01(e2);
            if (u <= r) {
                LCurr = LNext;
                x0ACurr=x0ANext;
                x0BCurr=x0BNext;
                x1ACurr=x1ANext;
                x1BCurr=x1BNext;
                UCurr = UNext;

            }//end of accept-reject
            U_dist_ptr[varNum*j+0]=UCurr;
            U_dist_ptr[varNum*j+1]=LCurr;
            U_dist_ptr[varNum*j+2]=x0ACurr;
            U_dist_ptr[varNum*j+3]=x0BCurr;
            U_dist_ptr[varNum*j+4]=x1ACurr;
            U_dist_ptr[varNum*j+5]=x1BCurr;
        }//end for loop
        size_t loopEnd = loopStart + loopToWrite - 1;
        std::string fileNameMiddle = "loopStart" + std::to_string(loopStart) + "loopEnd" + std::to_string(loopEnd);
        std::string out_U_distPickleFileName = this->U_dist_dataDir + "/" + fileNameMiddle + ".U_dist.csv";

        //save U_dist_ptr
        saveArrayToCSV(U_dist_ptr,varNum * loopToWrite,out_U_distPickleFileName,varNum);
        const auto tMCEnd{std::chrono::steady_clock::now()};
        const std::chrono::duration<double> elapsed_secondsAll{tMCEnd - tMCStart};
        std::cout << "loop " + std::to_string(loopStart) + " to loop " + std::to_string(loopEnd) + ": "
                  << elapsed_secondsAll.count() / 3600.0 << " h" << std::endl;

        loopStart = loopEnd + 1;
    }//end flush for loop

    std::cout<<"mc executed for "<<flushNum<<" flushes."<<std::endl;


}







void mc_computation::saveArrayToCSV(const std::shared_ptr<double[]>& array, const  size_t& arraySize, const std::string& filename, const size_t& numbersPerRow) {

    std::ofstream outFile(filename);

    if (!outFile.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }
    outFile << std::setprecision(std::numeric_limits<double>::digits10 + 1) << std::fixed;

    outFile<<"U,"<<"L,"<<"y0,"<<"z0,"<<"y1"<<"\n";
    for (size_t i = 0; i < arraySize; ++i) {
        outFile << array[i];
        if ((i + 1) % numbersPerRow == 0) {
            outFile << '\n';
        } else {
            outFile << ',';
        }
    }

    // If the last row isn't complete, end the line
    if (arraySize % numbersPerRow != 0) {
        outFile << '\n';
    }

    outFile.close();


}

void mc_computation::init_and_run(){
    this->execute_mc(LInit,x0AInit,x0BInit,x1AInit,x1BInit,loopLastFile+1,newFlushNum);


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
void mc_computation::mod_proposal(const double &LCurr, const double &x0ACurr, const double &x0BCurr, const double& x1ACurr,const double &x1BCurr,
              double &LNext, double &x0ANext, double &x0BNext, double &x1ANext, double &x1BNext){




    double lm = potFuncPtr->getLm();
    //for coordinates, the proposed value are Gaussian
    x0ANext=this->generate_nearby_normal(x0ACurr);
    x0BNext=this->generate_nearby_normal(x0BCurr);
    x1ANext=this->generate_nearby_normal(x1ACurr);
    x1BNext=this->generate_nearby_normal(x1BCurr);

    x0ANext=positive_fmod(x0ANext,lm);
    x0BNext=positive_fmod(x0BNext,lm);
    x1ANext=positive_fmod(x1ANext,lm);
    x1BNext=positive_fmod(x1BNext,lm);


    LNext= reject_sampling_one_data(LCurr,0,lm);

}


double mc_computation::generate_nearby_normal(const double &x){
    std::random_device rd;  // Create a random device object
    std::ranlux24_base engine(rd());  // Seed the engine with the random device

    std::normal_distribution<> normal_dist(x,h);

    double xNext=normal_dist(engine);
    return xNext;

}


///
/// @param y
/// @param x center
/// @return known proposal function, which is normal distribution
double mc_computation::Q(const double &y, const double &x, const double &a, const double &b){

    double val=1/(std::pow(2.0*PI,0.5)*h)
               *std::exp(-1/(2*std::pow(h,2))*std::pow(y-x,2.0));

    return val;

}


///
/// @param y
/// @param x center
/// @param a left end
/// @param b right end
/// @return truncated Gaussian
double mc_computation::f(const double &y, const double &x, const double &a, const double &b){


    if(y<=a or y>=b){
        return 0;
    }else{

        double val=std::exp(-1.0/(2.0*std::pow(h,2))*std::pow(y-x,2));
        return val;
    }

}

///
/// @param x center
/// @param a left end
/// @param b right end
/// @return random number from truncated Gaussian
double mc_computation::reject_sampling_one_data(const double &x,const double &a, const double &b){

    std::random_device rd;  // Create a random device object
    std::ranlux24_base engine(rd());  // Seed the engine with the random device

    std::normal_distribution<> normal_dist(x,h);
    std::uniform_real_distribution<> distUnif01(0, 1);//[0,1)
    double y=normal_dist(engine);
    double u=distUnif01(engine);

    while(u>=f(y,x,a,b)/(M* Q(y,x,a,b))){
        y=normal_dist(engine);
        u=distUnif01(engine);

    }

    return y;

}

///
/// @param x center
/// @param a left end
/// @param b right end
/// @return integral
double mc_computation::zVal(const double& x,const double &a, const double &b){

    auto integrandWithParam=[x,a,b, this](const double &y){
        return this->integrand(y,x,a,b);
    };
    double result = boost::math::quadrature::trapezoidal(integrandWithParam,a,b);

    return result;

}

///
/// @param y
/// @param x center
/// @param a left end
/// @param b right end
/// @return
double mc_computation::integrand(const double &y, const double& x,const double &a, const double &b){

    return f(y,x,a,b);


}


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
double mc_computation::mod_acceptanceRatio(const double &LCurr,const double &x0ACurr, const double &x0BCurr,
                       const double &x1ACurr, const double &x1BCurr,const double& UCurr,
                       const double & LNext, const double &x0ANext, const double &x0BNext,
                       const double &x1ANext, const double &x1BNext,double &UNext) {
    double lm = potFuncPtr->getLm();
    UNext = (*potFuncPtr)(LNext, x0ANext, x0BNext, x1ANext, x1BNext);
    double numerator = -this->beta * UNext;
    double denominator = -this->beta * UCurr;
    double R = std::exp(numerator - denominator);

    double zLCurr = zVal(LCurr, 0, lm);
    double zLNext = zVal(LNext, 0, lm);

    double ratio_L = zLCurr / zLNext;
    R *= ratio_L;
    return std::min(1.0, R);

}
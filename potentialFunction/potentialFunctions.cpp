//
// Created by polya on 7/19/24.
//

#include "potentialFunctionPrototype.hpp"

class V_inv_12_6:public potentialFunction {

public:
    V_inv_12_6(const std::string &coefsStr):potentialFunction(){
        this->coefsInStr=coefsStr;
    }

public:
    void json2Coefs(const std::string &coefsStr)override{
        std::stringstream iss;
        iss<<coefsStr;
        std::string temp;
        //read a1
        if (std::getline(iss, temp, ',')){
            this->a1=std::stod(temp);
        }

        //read b1

        if (std::getline(iss, temp, ',')){
            this->b1=std::stod(temp);
        }

        //read a2
        if (std::getline(iss, temp, ',')){
            this->a2=std::stod(temp);
        }

        //read b2

        if (std::getline(iss, temp, ',')){
            this->b2=std::stod(temp);
        }
    }

    void init() override{
        this->json2Coefs(coefsInStr);
        this->r1=std::pow(2.0*a1/b1,1.0/6.0);
        this->r2=std::pow(2.0*a2/b2,1.0/6.0);
        this->lm=(2*(r1+r2))*1.5;
        this->eps=((r1+r2)/2.0)/8;
        std::cout << "a1=" << a1 << ", b1=" << b1 << ", a2=" << a2 << ", b2=" << b2 << std::endl;
        std::cout<<"r1="<<r1<<", r2="<<r2<<std::endl;
        std::cout<<"lm="<<lm<<std::endl;
//        std::cout<<"eps="<<eps<<std::endl;


    }


    double operator()(const double &L, const double &x0A, const double &x0B, const double &x1A, const double& x1B) const override {
        double y0 = x0B - x0A;
        double y1 = x1B - x1A;
        double z0 = x1A - x0B;
//        std::cout<<"L="<<L<<", y0="<<y0<<", z0="<<z0<<", y1="<<y1<<std::endl;
        double val = a1 / std::pow(y0, 12.0) - b1 / std::pow(y0, 6.0)
                     + a2 / std::pow(z0, 12.0) - b2 / std::pow(z0, 6.0)
                     + a1 / std::pow(y1, 12.0) - b1 / std::pow(y1, 6.0)
                     + a2 / std::pow(-y0 - z0 - y1 + L, 12.0) - b2 / std::pow(-y0 - z0 - y1 + L, 6.0);

//        std::cout<<"a1/std::pow(y0,12.0)-b1/std::pow(y0,6.0)="<<a1/std::pow(y0,12.0)-b1/std::pow(y0,6.0)<<std::endl;
//        std::cout<<"a2/std::pow(z0,12.0)-b2/std::pow(z0,6.0)="<<a2/std::pow(z0,12.0)-b2/std::pow(z0,6.0)<<std::endl;
//        std::cout<<"a1/std::pow(y1,12.0)-b1/std::pow(y1,6.0)="<<a1/std::pow(y1,12.0)-b1/std::pow(y1,6.0)<<std::endl;
//        std::cout<<"a2/std::pow(-y0-z0-y1+L,12.0)-b2/std::pow(-y0-z0-y1+L,6.0)="<<a2/std::pow(-y0-z0-y1+L,12.0)-b2/std::pow(-y0-z0-y1+L,6.0)<<std::endl;
//        std::cout<<"val="<<val<<std::endl;
        return val;
    }
    double getLm() const override {
        return lm;
    }
    double get_eps() const override {
        return eps;
    }
public:
    double a1;
    double a2;
    double b1;
    double b2;
    std::string coefsInStr;
    double r1;//min position of V1
    double r2;//min position of V2
    double lm;//range of distances
    double eps;//half interval length of uniform distribution
};

std::shared_ptr<potentialFunction>  createPotentialFunction(const std::string& funcName, const std::string &coefsJsonStr) {
    if (funcName == "V_inv_12_6") {

        return std::make_shared<V_inv_12_6>(coefsJsonStr);
    }

    else {
        throw std::invalid_argument("Unknown potential function type");
    }
}
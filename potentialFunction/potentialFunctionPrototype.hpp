//
// Created by polya on 7/19/24.
//

#ifndef V_INV_12_6_POTENTIALFUNCTIONPROTOTYPE_HPP
#define V_INV_12_6_POTENTIALFUNCTIONPROTOTYPE_HPP
#include <fstream>
#include <iostream>
#include <math.h>
#include <memory>
#include <regex>
#include <stdexcept>
#include <string>


class potentialFunction {
//base class for potential function
public:
    virtual double operator()(const double&L,const double &x0A, const double &x0B, const double& x1A, const double &x1B) const = 0;
    virtual void json2Coefs(const std::string &coefsStr)=0;
    virtual  void init()=0;
    virtual double getLm() const = 0; //  method to get lm
    virtual double get_eps() const = 0; //  method to get eps
    virtual ~ potentialFunction() {};
};


std::shared_ptr<potentialFunction>  createPotentialFunction(const std::string& funcName, const std::string &row) ;

#endif //V_INV_12_6_POTENTIALFUNCTIONPROTOTYPE_HPP

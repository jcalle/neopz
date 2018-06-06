//
//  TestRandom.cpp
//  PZ
//
//  Created by Thiago Quinelato on 5/dec/17.
//  Copyright 2011 UNICAMP. All rights reserved.
//

#include "pzvec.h"
#include "TPZUniformRandom.h"
#include "TPZNormalRandom.h"
#include "TPZConstrainedNormalRandom.h"

#ifdef USING_BOOST

#ifndef WIN32
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MAIN pz random tests

#include <boost/test/unit_test.hpp>


void PrintDistribution(REAL& begin, unsigned int& steps, TPZVec<int>& slots, REAL& increment){
    REAL max = 0;
    for (unsigned int i = 0; i < steps; ++i) {
        max = std::max(max, (REAL)slots[i]);
    }
    for (unsigned int i = 0; i < steps; ++i) {
        int nstars = slots[i]*50/max;
        printf("%3.1f :", begin + i*increment);
        for (unsigned int j = 0; j < nstars; ++j) {
            std::cout << "*";
        }
        std::cout << std::endl;
    }    
    std::cout << std::endl;
}

void TestUniform(REAL begin, REAL end){
    unsigned int numbers = 1000;
    unsigned int steps = 10;
    TPZVec<int> slots(steps,0);
    REAL increment = (end-begin)/steps;
    TPZUniformRandom<REAL> r(begin, end);
    for (unsigned int i = 0; i < numbers; ++i) {
        REAL value = r.next();
        BOOST_ASSERT(value > begin);
        BOOST_ASSERT(value < end);
        slots[std::floor((value-begin)/increment)]++;
    }

//    for (unsigned int i = 0; i < steps; ++i) {
//      BOOST_ASSERT(slots[i] > .8*numbers/steps);
//  	BOOST_ASSERT(slots[i] < 1.2*numbers/steps);
//    }
}

void TestNormal(REAL mean, REAL stdev){
    unsigned int numbers = 1000;
    unsigned int steps = 20;
    TPZVec<int> slots(steps,0);
    REAL begin = mean - 5;
    REAL end = mean + 5;
    REAL increment = (end-begin)/steps;
    TPZNormalRandom<REAL> r(mean, stdev);
    for (unsigned int i = 0; i < numbers; ++i) {
        REAL value;
        do {
            value = r.next();
        } while (value <= begin || value >= end);
        slots[std::floor((value-begin)/increment)]++;
    }
    //PrintDistribution(begin, steps, slots, increment);
}

void TestConstrainedNormal(REAL begin, REAL end, REAL mean, REAL stdev){
    unsigned int numbers = 1000;
    unsigned int steps = 10;
    TPZVec<int> slots(steps,0);
    REAL increment = (end-begin)/steps;
    TPZConstrainedNormalRandom<REAL> r(begin, end, mean, stdev);
    for (unsigned int i = 0; i < numbers; ++i) {
        REAL value = r.next();
        BOOST_ASSERT(value > begin);
        BOOST_ASSERT(value < end);
        slots[std::floor((value-begin)/increment)]++;
    }
    //PrintDistribution(begin, steps, slots, increment);
}


BOOST_AUTO_TEST_SUITE(random_generator_tests)

BOOST_AUTO_TEST_CASE(uniform_random_generator_tests)
{
    TestUniform(0,1);
    TestUniform(1,2);
	TestUniform(-10,2);
}

BOOST_AUTO_TEST_CASE(normal_random_generator_tests)
{
    TestNormal(0,1);
	TestNormal(1,2);
	TestNormal(10,3);
}

BOOST_AUTO_TEST_CASE(constrained_normal_random_generator_tests)
{
    TestConstrainedNormal(-1,1,0.2,0.1);
	TestConstrainedNormal(-1,1,1,2);
    TestConstrainedNormal(0,10,12,3);
}

BOOST_AUTO_TEST_SUITE_END()

#endif

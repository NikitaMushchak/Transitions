#ifndef MATH_CONSTANTS_H_INCLUDED
#define MATH_CONSTANTS_H_INCLUDED
#include <math.h>
#include <fstream>
#include <istream>
#include <iostream>
#include <fstream>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <boost/qvm/mat.hpp>
#include <boost/qvm/mat_access.hpp>
#include <boost/qvm/vec.hpp>
#include <boost/qvm/map_vec_mat.hpp>
#include <boost/qvm/mat_operations.hpp>
#include <boost/qvm/mat_traits_array.hpp>
#include "DataTypes.h"
// Type of val --- unsigned long long

const double MC_s6 = sqrt(6.0);
const double MC_1d3 = 1.0/3.0;
const double MC_1d6 = 1.0/6.0;
const double MC_1d9 = 1.0/9.0;
const double MC_1d18 = 1.0/18.0;
const double MC_1d27 = 1.0/27.0;
const double MC_2d3 = 2.0/3.0;
const double MC_2d6 = 2.0/6.0;
const double MC_1ds3 = 1.0/sqrt(3.0);
const double MC_2ds3 = 2.0/sqrt(3.0);
const double MC_s3 = sqrt(3.0);
const double MC_s3d2 = sqrt(3.0)/2.0;
const double MC_s3d3 = sqrt(3.0)/3.0;
const double MC_s3d4 = sqrt(3.0)/4.0;
const double MC_s2 = sqrt(2.0);
const double MC_s2d2 = sqrt(2.0)/2.0;
const double MC_1ds2 = MC_s2d2;
const double MC_2s2 = 2.0*sqrt(2.0);
const double MC_2s3 = 2.0*sqrt(3.0);
const double MC_2s2d3 = 2.0*sqrt(2.0)/3.0;
const double MC_4s2d3 = 4.0*sqrt(2.0)/3.0;
const double cos30 = MC_s3d2;
const double cos45 = MC_s2d2;
const double cos60 = 0.5;
const double sin30 = 0.5;
const double sin45 = MC_s2d2;
const double sin60 = MC_s3d2;
const double tan30 = sin30/cos30;
const double tan45 = sin45/cos45;
const double tan60 = sin60/cos60;
const double pi = 3.1415926535897932384626433832795;
const double e = 2.7182818284590452353602874713527;
const double min_double = -1.7e-308;
const double max_double = 1.7e+308;
const double zero_double = 1e-15;
const double MC_1d72 = 1.0/72.0;
const double RAND_MAX_double = (double)RAND_MAX;
const double _1d_RAND_MAX_double = 1.0/RAND_MAX_double;
//const uint_fast32_t UINT64_MAX = 0xffffffffffffffff;
const double near_zero = 1e-14;

struct Tensor2D
{
    double xx,xy,yy;
};

struct Tensor2D_full
{
    double xx,xy,yx,yy;
};

struct Vector2Di
{
    int_fast32_t x,y;
};

struct Vector2D
{
    double x,y;
    bool operator ==(const Vector2D& obj) const
    {
        return (x == obj.x && y == obj.y)? true:false;
    }
};
struct Vector2Dui
{
    unsigned int x,y;
    bool operator ==(const Vector2Dui& obj) const
    {
        return (x == obj.x && y == obj.y)? true:false;
    }
};
struct Tensor_3r
{
    double xx,yy,zz,xy,xz,yz;
};

const Vector2D ZeroVector2D= {0,0};
const Vector2D zero_vector = {0,0};
const Vector2D vector_11 = {1,1};
const Vector2D vector_10 = {1,0};
const Vector2D vector_01 = {0,1};

const Vector2Dui ZeroVector2Dui = {0,0};

const Tensor_3r zero_tensor = {0,0,0,0,0,0};

//const boost::qvm::mat<double,3,3> IndentityT33= {1.0,0,0, 0,1.0,0, 0,0,1.0};
const boost::qvm::mat<double,2,2> MC_1T22= {1.0,0, 0,1.0};
const boost::qvm::mat<double,2,2> MC_0T22= {0,0,0,0};
const boost::qvm::vec<int_fast32_t,3> MC_i1XYZV3= {1,1,1};
const boost::qvm::vec<int_fast32_t,3> MC_i110XYZV3= {1,1,0};

const boost::qvm::mat<double,3,3> MC_1T33= {1.0,0,0, 0,1.0,0, 0,0,1.0};
const boost::qvm::mat<double,3,3> MC_0T33= {0,0,0, 0,0,0, 0,0,0};


boost::qvm::mat<double,2,2> tens(const boost::qvm::vec<double,2> &a, const boost::qvm::vec<double,2> &b);
boost::qvm::mat<double,4,4> tens(const boost::qvm::mat<double,2,2> &a, const boost::qvm::mat<double,2,2> &b);
boost::qvm::mat<double,2,2> dotdot(const boost::qvm::mat<double,4,4> &a, const boost::qvm::mat<double,2,2> &b);

boost::qvm::mat<double,3,3> tens(const boost::qvm::vec<double,3> &a, const boost::qvm::vec<double,3> &b);
boost::qvm::mat<double,9,9> tens(const boost::qvm::mat<double,3,3> &a, const boost::qvm::mat<double,3,3> &b);
boost::qvm::mat<double,3,3> dotdot(const boost::qvm::mat<double,9,9> &a, const boost::qvm::mat<double,3,3> &b);


void toVoigtNotationC2D(const boost::qvm::mat<double,4,4> &a, boost::qvm::mat<double,3,3> &b);
void fromVoigtNotationC2D(const boost::qvm::mat<double,3,3> &a, boost::qvm::mat<double,4,4> &b);
void toVoigtNotationS2D(const boost::qvm::mat<double,4,4> &a, boost::qvm::mat<double,3,3> &b);
void fromVoigtNotationS2D(const boost::qvm::mat<double,3,3> &a, boost::qvm::mat<double,4,4> &b);
void toElastic2D(const boost::qvm::mat<double,3,3> &a, ElasticModules2D &b);


void toVoigtNotationC3D(const boost::qvm::mat<double,9,9> &a, boost::qvm::mat<double,6,6> &b);
void fromVoigtNotationC3D(const boost::qvm::mat<double,6,6> &a, boost::qvm::mat<double,9,9> &b);
void toVoigtNotationS3D(const boost::qvm::mat<double,9,9> &a, boost::qvm::mat<double,6,6> &b);
void fromVoigtNotationS3D(const boost::qvm::mat<double,6,6> &a, boost::qvm::mat<double,9,9> &b);
void toElastic3D(const boost::qvm::mat<double,6,6> &a, ElasticModules3D &b);
#endif // MATH_CONSTANTS_H_INCLUDED

//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// Planet - An atmospheric code for planetary bodies, adapted to Titan
//
// Copyright (C) 2013 The PECOS Development Team
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the Version 2.1 GNU Lesser General
// Public License as published by the Free Software Foundation.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor,
// Boston, MA  02110-1301  USA
//
//-----------------------------------------------------------------------el-

#ifndef _ABSORPTION_GRID_
#define _ABSORPTION_GRID_

//Antioch
#include "antioch/metaprogramming.h"

//Planet

//C++
#include <vector>
#include <iostream>

namespace Planet{

template <typename CoeffType = double, typename VectorCoeffType = std::vector<double> >
class AbsorptionGrid
{
     public:
        AbsorptionGrid(const VectorCoeffType &x, const VectorCoeffType &y);
        AbsorptionGrid();
        ~AbsorptionGrid();

        //!
        void set_x_grid(const VectorCoeffType &x);
        //!
        void set_y_grid(const VectorCoeffType &y);

        //!
        void y_on_x_grid(const VectorCoeffType &x);

        //!
        const VectorCoeffType x() const {return _x;}
        //!
        const VectorCoeffType y() const {return _y;}
        //!
        const VectorCoeffType y_on_custom() const {return _y_on_custom;}

     private:
        VectorCoeffType _x;
        VectorCoeffType _y;
        VectorCoeffType _y_on_custom;
};

template <typename CoeffType, typename VectorCoeffType>
inline
AbsorptionGrid<CoeffType,VectorCoeffType>::AbsorptionGrid()
{
  return;
}

template <typename CoeffType, typename VectorCoeffType>
inline
AbsorptionGrid<CoeffType,VectorCoeffType>::AbsorptionGrid(const VectorCoeffType &x, const VectorCoeffType &y)
{
  set_x_grid(x);
  set_y_grid(y);
  return;
}

template <typename CoeffType, typename VectorCoeffType>
inline
AbsorptionGrid<CoeffType,VectorCoeffType>::~AbsorptionGrid()
{
  return;
}

template <typename CoeffType, typename VectorCoeffType>
inline
void AbsorptionGrid<CoeffType,VectorCoeffType>::set_x_grid(const VectorCoeffType &x)
{
 _x = x;
}

template <typename CoeffType, typename VectorCoeffType>
inline
void AbsorptionGrid<CoeffType,VectorCoeffType>::set_y_grid(const VectorCoeffType &y)
{
  _y = y;
  _y_on_custom = _y;
}

template <typename CoeffType, typename VectorCoeffType>
inline
void AbsorptionGrid<CoeffType,VectorCoeffType>::y_on_x_grid(const VectorCoeffType &x)
{
  
  _y_on_custom.clear();
  _y_on_custom.resize(x.size(),0.L);

  unsigned int ilow(0);
  while(x[ilow] <= _x[0])ilow++; //skipping too low bins

  unsigned int j(1);
  for(unsigned int i = ilow; i < x.size() - 1; i++)//bin per bin, right stairs: _y[i] from _x[i] to _x[i+1]
  {
     while(_x[j] < x[i]) //find lowest j / _x[j] > x[i]
     {
       j++;
      if(j >= _x.size() - 1)return;
     }
     CoeffType bin;
     Antioch::set_zero(bin);

//here we are: _x[j-1] =< x[i] < _x[j] with j-1 >= 0

     //targeted bin within stored bin: _x[j-1] =< x[i] < x[i+1] =< _x[j], take all of them
     while(_x[j] >= x[i+1])
     {
       _y_on_custom[i] = _y[j-1]; //rectangle from i to i+1, same height
       i = i + 1;
       if(i >= x.size() - 2)break;
     }

     //_x[j-1] < x[i] < _x[j] < x[i+1], calculating rectangle from x[i] to _x[j], height is _y[j-1]
     bin = _y[j-1] * (_x[j] - x[i]); 

     j++;
     while(_x[j] < x[i+1])// finding lowest j / x[i+1] < _x[j], adding all the k cases _x[j-1] < x[i] < _x[j] < _x[j+1] < ... < _x[j+k] < x[i+1]
     {
        bin += _y[j-1] * (_x[j] - _x[j-1]); // adding contained bins
        j++;
        if(j >= _x.size()-1)break;
     }

//now we have found n, we calculate rectangle from _x[j+n-1] to x[i+1]
     if(j < _x.size() - 1)bin += _y[j-2] * (x[i+1] - _x[j-1]); //if exist, above rectangle
     _y_on_custom[i] = bin/(x[i+1] - x[i]); //rectangle from i to i+1

  }

  return;
}

}

#endif

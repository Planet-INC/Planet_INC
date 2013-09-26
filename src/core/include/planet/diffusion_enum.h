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

#ifndef _DIFFUSION_ENUM_
#define _DIFFUSION_ENUM_

//Antioch
//Planet
//C++

namespace Planet
{
  enum DiffusionType
        {
          Massman = 0, /* Dij(D01,beta) = D01[0 K,1 atm] * P_normal/P * (T/T_standard)^beta */
          Wilson,      /* Dij(A,s) = A * (T/Tw[1K])^s */
          Wakeham,     /* Dij(A,s) = A * T^s/n */
          NoData       /* Dij(Mi,Mj) = Dii * sqrt((Mj/Mi+1)/2)  if Mj < Mi
                                     = Dii * sqrt(Mj/Mi)        if Mj >= Mi */
        };
}

#endif

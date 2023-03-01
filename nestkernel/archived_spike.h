/*
 *  ArchivedSpikeTrace.h
 *
 *  This file is part of NEST.
 *
 *  Copyright (C) 2004 The NEST Initiative
 *
 *  NEST is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  NEST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with NEST.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

/**
 * \file ArchivedSpikeTrace.h
 * Part of definition of ArchivingNode which is capable of recording and managing a spike history.
 * \author Moritz Helias, Abigail Morrison
 * \note moved to separate file to avoid circular inclusion in node.h
 * \date april 2006
 */

#ifndef HISTENTRY_H
#define HISTENTRY_H

// Includes from nestkernel:
#include "nest_types.h"

namespace nest
{

class ArchivedSpikeBase
{
  public:
    ArchivedSpikeBase( double t )
      : t( t ){}

    double t;               //!< point in time when spike occurred (in ms)
};

// entry in the spiking history
class ArchivedSpikeTrace : public ArchivedSpikeBase
{
public:
  ArchivedSpikeTrace( double t, double Kminus, double Kminus_triplet )
    : ArchivedSpikeBase( t )
    , Kminus( Kminus )
    , Kminus_triplet( Kminus_triplet ){}

  double Kminus;          //!< value of Kminus at that time
  double Kminus_triplet;  //!< value of triplet STDP Kminus at that time
};

// entry in the history of plasticity rules which consider additional factors
class ArchivedSpikeGeneric : public ArchivedSpikeBase
{
public:
  ArchivedSpikeGeneric( double t, double dw )
    : ArchivedSpikeBase( t )
    , value( dw ){}

  double value;              //!< value dependent on the additional factor
};
}

#endif

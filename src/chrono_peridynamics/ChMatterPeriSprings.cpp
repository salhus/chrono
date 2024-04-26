// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2014 projectchrono.org
// All rights reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
// Authors: Alessandro Tasora, Radu Serban 
// =============================================================================

#include <cstdlib>
#include <algorithm>

#include "chrono_peridynamics/ChMatterPeriSprings.h"
#include "chrono_peridynamics/ChProximityContainerPeridynamics.h"

namespace chrono {

using namespace fea;
using namespace peridynamics;
using namespace geometry;

// Register into the object factory, to enable run-time dynamic creation and persistence
CH_FACTORY_REGISTER(ChMatterPeriSprings)

// Register into the object factory, to enable run-time dynamic creation and persistence
CH_FACTORY_REGISTER(ChMatterPeriSpringsBreakable)



} // end namespace chrono

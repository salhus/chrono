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
// Authors: Radu Serban
// =============================================================================

#ifndef CH_COLLISION_SHAPE
#define CH_COLLISION_SHAPE

#include "chrono/core/ChApiCE.h"
#include "chrono/physics/ChMaterialSurface.h"

namespace chrono {
namespace collision {

/// @addtogroup chrono_collision
/// @{

/// Class defining a collision shape.
class ChApi ChCollisionShape {
  public:
    /// Supported collision shape types.
    enum Type {
        SPHERE,
        ELLIPSOID,
        BOX,
        CYLINDER,
        CYLSHELL,
        CONVEXHULL,
        TRIANGLEMESH,
        BARREL,  // Not supported in Chrono collision system
        POINT,
        TRIANGLE,
        CAPSULE,
        CONE,         // Not implemented in Bullet collision system
        ROUNDEDBOX,   // Not implemented in Bullet collision system
        ROUNDEDCYL,   // Not implemented in Bullet collision system
        TETRAHEDRON,  // Not implemented in Bullet collision system
        PATH2D,
        UNKNOWN_SHAPE
    };

    ChCollisionShape(Type type = Type::UNKNOWN_SHAPE);
    ChCollisionShape(Type type, std::shared_ptr<ChMaterialSurface> material);
    virtual ~ChCollisionShape() {}

    Type GetType() const { return m_type; }

    std::shared_ptr<ChMaterialSurface> GetMaterial() const { return m_material; }
    ChContactMethod GetContactMethod() const { return m_material->GetContactMethod(); }

    /// Method to allow serialization of transient data to archives.
    virtual void ArchiveOut(ChArchiveOut& marchive);

    /// Method to allow de-serialization of transient data from archives.
    virtual void ArchiveIn(ChArchiveIn& marchive);

  protected:
    Type m_type;                                    ///< type of collision shape
    std::shared_ptr<ChMaterialSurface> m_material;  ///< surface contact material

    friend class ChCollisionModel;
};

/// @} chrono_collision

}  // end namespace collision
}  // end namespace chrono

#endif

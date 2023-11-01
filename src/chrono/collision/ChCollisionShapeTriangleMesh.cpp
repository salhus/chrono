// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2023 projectchrono.org
// All rights reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
// Authors: Radu Serban
// =============================================================================

#include "chrono/collision/ChCollisionShapeTriangleMesh.h"

namespace chrono {
namespace collision {

ChCollisionShapeTriangleMesh::ChCollisionShapeTriangleMesh() : ChCollisionShape(Type::TRIANGLEMESH), trimesh(nullptr){};

ChCollisionShapeTriangleMesh::ChCollisionShapeTriangleMesh(std::shared_ptr<ChMaterialSurface> material,
                                                           std::shared_ptr<geometry::ChTriangleMesh> mesh,
                                                           bool is_static,
                                                           bool is_convex,
                                                           double thickness)
    : ChCollisionShape(Type::TRIANGLEMESH, material), trimesh(mesh) {
    this->is_static = is_static;
    this->is_convex = is_convex;
    this->thickness = thickness;
}

void ChCollisionShapeTriangleMesh::ArchiveOut(ChArchiveOut& marchive) {
    // version number
    marchive.VersionWrite<ChCollisionShapeTriangleMesh>();
    // serialize parent class
    ChCollisionShape::ArchiveOut(marchive);
    // serialize all member data:
    marchive << CHNVP(trimesh);
    marchive << CHNVP(is_static);
    marchive << CHNVP(is_convex);
    marchive << CHNVP(thickness);
}

void ChCollisionShapeTriangleMesh::ArchiveIn(ChArchiveIn& marchive) {
    // version number
    /*int version =*/marchive.VersionRead<ChCollisionShapeTriangleMesh>();
    // deserialize parent class
    ChCollisionShape::ArchiveIn(marchive);
    // stream in all member data:
    marchive >> CHNVP(trimesh);
    marchive >> CHNVP(is_static);
    marchive >> CHNVP(is_convex);
    marchive >> CHNVP(thickness);
}

}  // end namespace collision
}  // end namespace chrono

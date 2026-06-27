/*******************************************************************************
 * This file is part of Tissue Forge.
 * Copyright (c) 2022-2024 T.J. Sego and Tien Comlekoglu
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

/**
 * @file tfMeshReconnection.h
 *
 * @brief Native 3D T1 / reversible network reconnection (RNR): the Okuda I<->H face<->edge swap
 *        (Okuda et al. 2013, Biomech Model Mechanobiol 12:627-644). This is the 3D topology
 *        operation TissueForge otherwise lacks (its other 3D quality ops are only degenerate
 *        collapses). The RNR neighborhood walks, Condition-4 vetoes, Okuda Appendix-1 vertex
 *        placement, and surface surgery live in tfMeshReconnection.cpp; the MeshQuality scheduler
 *        (tfMeshQuality.cpp) drives them via the single entry point below.
 */

#ifndef _MODELS_VERTEX_SOLVER_TFMESHRECONNECTION_H_
#define _MODELS_VERTEX_SOLVER_TFMESHRECONNECTION_H_

#include "tfMeshQuality.h"

#include <tf_port.h>

#include <vector>


namespace TissueForge::models::vertex {


    class Mesh;


    /**
     * @brief Build the native RNR (Okuda I<->H) reconnection operations for the current mesh.
     *
     * Scans the mesh for reconnection candidates at the given trigger length (Okuda Condition 2,
     * Delta_l_th) -- short interior edges (I->H) and small triangular faces (H->I) -- applies the
     * Condition-4 vetoes, and creates one ReconnectionOperation per surviving site, de-duplicated
     * and indexed by trigger-surface id. The caller (MeshQuality::doQuality's reconnection pass)
     * then builds the dependency chains and executes them.
     *
     * @param mesh the mesh to scan
     * @param passMask per-surface mask; masked (excluded) surfaces are skipped
     * @param reconnectLength trigger length (Okuda Delta_l_th); <= 0 leaves ops empty (no-op)
     * @param reconnectHysteresis placement hysteresis for created features
     * @param reconnectEnergyGate optional greedy energy gate (a DEPARTURE from Okuda; default off)
     * @param ops [out] sized to mesh->sizeSurfaces(); filled with new ReconnectionOperation* (as
     *            MeshQualityOperation*) at each trigger-surface slot. Ownership passes to the caller.
     */
    HRESULT MeshReconnection_buildOperations(
        Mesh *mesh,
        const std::vector<bool> &passMask,
        const FloatP_t &reconnectLength,
        const FloatP_t &reconnectHysteresis,
        const bool &reconnectEnergyGate,
        std::vector<MeshQualityOperation*> &ops
    );


}

#endif // _MODELS_VERTEX_SOLVER_TFMESHRECONNECTION_H_

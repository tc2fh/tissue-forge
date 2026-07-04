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
 * C++ gate for the 3D vertex-model reversible network reconnection (RNR) MeshQuality op.
 *
 * Exercises the native reconnection C++ API directly so an in-tree build covers the public C++
 * surface (the exhaustive behavioral coverage -- that RNR actually reconnects, sorts, is
 * reversible, etc. -- lives in the Python gates, rnr/tests). Builds a cubic lattice of cells
 * (createPLPDMesh -- an interior edge is shared by four cells, the Okuda [I] configuration), then:
 *   1. round-trips the RNR + per-pass + serialization API surface (reconnect length/hysteresis/
 *      energy-gate/interval, the per-pass stock_{vertex,surface,body}_operations flags, the
 *      volume-repair flag, and the active-motility drive parameters);
 *   2. runs the live pipeline -- the native RNR reconnection pass coexisting with the safe stock
 *      quality passes (body-demote excluded, per H3) plus the active-motility drive -- for several
 *      steps, asserting the whole vertex sub-engine advances inside Universe::step() with no error
 *      and leaves the mesh valid (all cell volumes positive).
 *
 * NB the physical force actors (VolumeConstraint, etc.) are not CAPI_EXPORT, so they cannot be
 * constructed from an external translation unit; the Python gates cover the force-driven behavior.
 */

#include "tfTest.h"

#include <tfUniverse.h>

#include <models/vertex/solver/tfMeshSolver.h>
#include <models/vertex/solver/tfMesh.h>
#include <models/vertex/solver/tfMeshQuality.h>
#include <models/vertex/solver/tfBody.h>
#include <models/vertex/solver/tfSurface.h>
#include <models/vertex/solver/tf_mesh_create.h>
#include <models/vertex/solver/actors/tfVolumeConstraint.h>

#include <cmath>


using namespace TissueForge;
namespace vm = TissueForge::models::vertex;


#define TF_TEST_EXPECT(cond) { if(!(cond)) { TF_TEST_REPORTERR(); return E_FAIL; } }


HRESULT test_rnr_api_and_live_pipeline() {
    Simulator::Config conf;
    conf.setWindowless(true);
    conf.universeConfig.dim = FVector3(20, 20, 20);
    conf.universeConfig.dt = 0.01;
    TF_TEST_CHECK(tfTest_init(conf));
    TF_TEST_CHECK(vm::MeshSolver::init());

    vm::Mesh *mesh = vm::Mesh::get();
    TF_TEST_EXPECT(mesh != NULL);

    // --- types: a cell and its interface surface ------------------------------------------
    // Construct with noReg=true so the type is NOT auto-registered by its constructor; we
    // register it explicitly below to exercise the MeshSolver::registerType API path.
    vm::BodyType *btype = new vm::BodyType(true);
    btype->name = "Cell";
    vm::SurfaceType *stype = new vm::SurfaceType(true);
    stype->name = "Interface";
    TF_TEST_CHECK(vm::MeshSolver::registerType(btype));
    TF_TEST_CHECK(vm::MeshSolver::registerType(stype));

    // --- 4x4x4 cubic lattice of cells -----------------------------------------------------
    auto grid = vm::createPLPDMesh(btype, stype, FVector3(3, 3, 3), 4, 4, 4, 1.0, 1.0, 1.0);
    std::vector<vm::BodyHandle> bodies;
    for(auto &plane : grid)
        for(auto &row : plane)
            for(auto &b : row)
                bodies.push_back(b);
    TF_TEST_EXPECT(bodies.size() == 64);

    // --- RNR + per-pass + serialization API round-trip ------------------------------------
    vm::MeshQuality q;
    TF_TEST_CHECK(q.setReconnectLength(0.2));
    TF_TEST_EXPECT(std::abs(q.getReconnectLength() - 0.2) < 1e-6);
    TF_TEST_CHECK(q.setReconnectHysteresis(0.2));
    TF_TEST_EXPECT(std::abs(q.getReconnectHysteresis() - 0.2) < 1e-6);
    TF_TEST_CHECK(q.setReconnectEnergyGate(false));
    TF_TEST_EXPECT(q.getReconnectEnergyGate() == false);
    TF_TEST_CHECK(q.setReconnectInterval(2));
    TF_TEST_EXPECT(q.getReconnectInterval() == 2);
    TF_TEST_CHECK(q.setReconnectInterval(1));

    // H3 per-pass stock-op flags: keep the safe passes, exclude only the body-demote pass
    TF_TEST_CHECK(q.setStockQualityOps(true));
    TF_TEST_EXPECT(q.getStockQualityOps() == true);
    TF_TEST_CHECK(q.setStockVertexOps(true));
    TF_TEST_EXPECT(q.getStockVertexOps() == true);
    TF_TEST_CHECK(q.setStockSurfaceOps(true));
    TF_TEST_EXPECT(q.getStockSurfaceOps() == true);
    TF_TEST_CHECK(q.setStockBodyOps(false));
    TF_TEST_EXPECT(q.getStockBodyOps() == false);

    // H5 volume-orientation repair flag round-trip
    TF_TEST_CHECK(vm::MeshSolver::setVolumeRepair(false));
    TF_TEST_EXPECT(vm::MeshSolver::getVolumeRepair() == false);
    TF_TEST_CHECK(vm::MeshSolver::setVolumeRepair(true));
    TF_TEST_EXPECT(vm::MeshSolver::getVolumeRepair() == true);

    // H5 active-motility drive parameters
    TF_TEST_CHECK(vm::MeshSolver::setMotility(0.2, 1.0, 5));
    TF_TEST_EXPECT(std::abs(vm::MeshSolver::getMotilityV0() - 0.2) < 1e-5);
    TF_TEST_EXPECT(std::abs(vm::MeshSolver::getMotilityDr() - 1.0) < 1e-5);
    TF_TEST_EXPECT(vm::MeshSolver::getMotilitySeeded() == true);
    TF_TEST_EXPECT(vm::MeshSolver::getMotilityRngState().size() > 0);

    // --- live pipeline: the RNR reconnection pass + safe stock passes + motility -----------
    vm::MeshQuality *quality = new vm::MeshQuality();
    quality->setReconnectLength(0.2);
    quality->setReconnectHysteresis(0.2);
    quality->setStockQualityOps(true);
    quality->setStockBodyOps(false);   // H3: exclude only the finite-block-unsafe body-demote pass
    TF_TEST_CHECK(mesh->setQuality(quality));

    for(int i = 0; i < 30; i++)
        TF_TEST_CHECK(Universe::step());

    // the vertex sub-engine advanced with no error and the mesh stayed valid
    unsigned int nliving = 0;
    for(auto &b : bodies) {
        vm::Body *body = b.body();
        if(body && body->objectId() >= 0) {
            nliving++;
            TF_TEST_EXPECT(body->getVolume() > 0);
        }
    }
    TF_TEST_EXPECT(nliving == bodies.size());

    TF_Log(LOG_INFORMATION) << "vertex RNR C++ gate passed: living cells = " << nliving
                            << ", vertices = " << mesh->sizeVertices();
    return S_OK;
}


int main(int argc, char const *argv[]) {
    if(test_rnr_api_and_live_pipeline() != S_OK) {
        TF_TEST_REPORTERR();
        return 1;
    }
    return 0;
}

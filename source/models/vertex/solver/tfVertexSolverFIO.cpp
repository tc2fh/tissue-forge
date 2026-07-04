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

#include "tfVertexSolverFIO.h"

#include "tfMeshSolver.h"
#include "tf_mesh_io.h"

#include <tfError.h>


using namespace TissueForge;
namespace VMod = TissueForge::models::vertex;


std::string VMod::io::VertexSolverFIOModule::moduleName() { return "VertexSolver"; }

HRESULT VMod::io::VertexSolverFIOModule::toFile(const TissueForge::io::MetaData &metaData, TissueForge::io::IOElement &fileElement) {

    VMod::MeshSolver *solver = VMod::MeshSolver::get();

    if(!solver) 
        return tf_error(E_FAIL, "toFile requested without a solver");

    // Store types

    std::vector<VMod::BodyType> bodyTypes;
    std::vector<int> bodyTypeIds;
    for(size_t i = 0; i < solver->numBodyTypes(); i++) {
        auto o = solver->getBodyType(i);
        bodyTypes.push_back(*o);
        bodyTypeIds.push_back(o->id);
    }
    TF_IOTOEASY(fileElement, metaData, "bodyTypes", bodyTypes);
    TF_IOTOEASY(fileElement, metaData, "bodyTypeIds", bodyTypeIds);

    std::vector<VMod::SurfaceType> surfaceTypes;
    std::vector<int> surfaceTypeIds;
    for(size_t i = 0; i < solver->numSurfaceTypes(); i++) {
        auto o = solver->getSurfaceType(i);
        surfaceTypes.push_back(*o);
        surfaceTypeIds.push_back(o->id);
    }
    TF_IOTOEASY(fileElement, metaData, "surfaceTypes", surfaceTypes);
    TF_IOTOEASY(fileElement, metaData, "surfaceTypeIds", surfaceTypeIds);
    TF_IOTOEASY(fileElement, metaData, "mesh", solver->mesh);

    // Native active-motility drive (PORTING_NOTES §6o): the per-cell directors are saved
    // with each Body; here we save the solver-level drive parameters and the RNG stream so
    // that a restored run continues the exact director trajectory rather than re-seeding.
    TF_IOTOEASY(fileElement, metaData, "motilityV0", VMod::MeshSolver::getMotilityV0());
    TF_IOTOEASY(fileElement, metaData, "motilityDr", VMod::MeshSolver::getMotilityDr());
    TF_IOTOEASY(fileElement, metaData, "motilitySeeded", VMod::MeshSolver::getMotilitySeeded());
    TF_IOTOEASY(fileElement, metaData, "motilityRngState", VMod::MeshSolver::getMotilityRngState());

    return S_OK;
}

HRESULT VMod::io::VertexSolverFIOModule::fromFile(const TissueForge::io::MetaData &metaData, const TissueForge::io::IOElement &fileElement) {

    VMod::MeshSolver *solver = VMod::MeshSolver::get();
    if(!solver && (VMod::MeshSolver::init() != S_OK || (solver = VMod::MeshSolver::get()) == NULL)) 
        return tf_error(E_FAIL, "Could not initialize solver");

    if(clearImport() != S_OK) 
        return tf_error(E_FAIL, "Failed clearing previous import");
    VMod::io::VertexSolverFIOModule::importSummary = new VMod::io::VertexSolverFIOImportSummary();

    // Load and register types and store their id maps

    std::vector<int> bodyTypeIds;
    TF_IOFROMEASY(fileElement, metaData, "bodyTypeIds", &bodyTypeIds);
    for(size_t i = 0; i < bodyTypeIds.size(); i++) 
        VMod::io::VertexSolverFIOModule::importSummary->bodyTypeIdMap.insert({bodyTypeIds[i], solver->numBodyTypes() + i});
    std::vector<VMod::BodyType*> bodyTypes;
    TF_IOFROMEASY(fileElement, metaData, "bodyTypes", &bodyTypes);
    for(size_t i = 0; i < bodyTypes.size(); i++) 
        solver->registerType(bodyTypes[i]);

    std::vector<int> surfaceTypeIds;
    TF_IOFROMEASY(fileElement, metaData, "surfaceTypeIds", &surfaceTypeIds);
    for(size_t i = 0; i < surfaceTypeIds.size(); i++) 
        VMod::io::VertexSolverFIOModule::importSummary->surfaceTypeIdMap.insert({surfaceTypeIds[i], solver->numSurfaceTypes() + i});
    std::vector<VMod::SurfaceType*> surfaceTypes;
    TF_IOFROMEASY(fileElement, metaData, "surfaceTypes", &surfaceTypes);
    for(size_t i = 0; i < surfaceTypes.size(); i++) 
        solver->registerType(surfaceTypes[i]);

    // Load meshes

    TF_IOFROMEASY(fileElement, metaData, "mesh", solver->mesh);

    // Restore the native active-motility drive (guarded for backward compatibility with
    // files written before motility was serialized). Directors themselves are restored
    // per-cell in Body::fromFile; here we restore the drive parameters + the RNG stream.
    TissueForge::io::IOChildMap fec = TissueForge::io::IOElement::children(fileElement);
    if(fec.find("motilityV0") != fec.end()) {
        FloatP_t motilityV0 = 0, motilityDr = 1;
        bool motilitySeeded = false;
        std::string motilityRngState;
        TF_IOFROMEASY(fileElement, metaData, "motilityV0", &motilityV0);
        TF_IOFROMEASY(fileElement, metaData, "motilityDr", &motilityDr);
        TF_IOFROMEASY(fileElement, metaData, "motilitySeeded", &motilitySeeded);
        if(fec.find("motilityRngState") != fec.end())
            TF_IOFROMEASY(fileElement, metaData, "motilityRngState", &motilityRngState);
        VMod::MeshSolver::restoreMotility(motilityV0, motilityDr, motilitySeeded, motilityRngState);
    }

    // Refresh derived geometry (centroids/areas/volumes/orientSign) so the restored mesh is
    // immediately consistent. Without this the bodies carry stale values computed while the
    // mesh was still being assembled (before periodic_geometry was set), so the first read
    // of a restored volume/centroid would be wrong until the next step's positionChanged.
    VMod::MeshSolver::positionChanged();

    return S_OK;
}

bool VMod::io::VertexSolverFIOModule::hasImport() { return VMod::io::VertexSolverFIOModule::importSummary != NULL; }

HRESULT VMod::io::VertexSolverFIOModule::clearImport() {
    if(hasImport()) {
        delete VMod::io::VertexSolverFIOModule::importSummary;
        VMod::io::VertexSolverFIOModule::importSummary = NULL;
    }
    return S_OK;
}

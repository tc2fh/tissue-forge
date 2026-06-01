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

#include "tfMeshQuality.h"

#include "tfMesh.h"
#include "tfMeshSolver.h"
#include "tf_mesh_metrics.h"
#include "tf_mesh_io.h"
#include "tfVertexSolverFIO.h"
#include "tf_mesh_ops.h"

#include <tfError.h>
#include <tfUniverse.h>
#include <tf_metrics.h>
#include <tfTaskScheduler.h>
#include <io/tfIO.h>
#include <io/tfFIO.h>

#include <Magnum/Math/Math.h>

#include <algorithm>
#include <atomic>
#include <array>
#include <map>
#include <sstream>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>


using namespace TissueForge;
using namespace TissueForge::models::vertex;


//////////////////////////
// MeshQualityOperation //
//////////////////////////


HRESULT MeshQualityOperation_checkChain( 
    MeshQualityOperation *op, 
    std::vector<MeshQualityOperation*> &ops
) {
    for(auto &t : op->targets) {
        MeshQualityOperation *t_op = ops[t];
        if(t_op && op != t_op) 
            op->appendNext(t_op);
        
    }

    return S_OK;
}

MeshQualityOperation::MeshQualityOperation(Mesh *_mesh) : 
    flags{Flag::None}, 
    mesh{_mesh}
{}

HRESULT MeshQualityOperation::appendNext(MeshQualityOperation *_next) {
    // Prevent loops
    std::set<MeshQualityOperation*> us = this->upstreams();
    if(std::find(us.begin(), us.end(), _next) != us.end()) 
        return S_OK;

    _next->lock.lock();
    _next->prev.insert(this);
    _next->lock.unlock();

    next.insert(_next);
    
    return S_OK;
}

HRESULT MeshQualityOperation::removeNext(MeshQualityOperation *_next) {
    auto itr = std::find(next.begin(), next.end(), _next);
    if(itr == next.end()) 
        return E_FAIL;
    next.erase(itr);

    itr = std::find(_next->prev.begin(), _next->prev.end(), this);
    if(itr != _next->prev.end()) {
        _next->lock.lock();
        _next->prev.erase(itr);
        _next->lock.unlock();
    }

    return S_OK;
}

static void MeshQuality_upstreams(const MeshQualityOperation *op, std::set<MeshQualityOperation*> &ops) {
    for(auto *op_u : op->prev) {
        if(std::find(ops.begin(), ops.end(), op_u) == ops.end()) {
            ops.insert(op_u);
            MeshQuality_upstreams(op_u, ops);
        }
    }
}

std::set<MeshQualityOperation*> MeshQualityOperation::upstreams() const {
    std::set<MeshQualityOperation*> result;
    MeshQuality_upstreams(this, result);
    return result;
}

static void MeshQuality_downstreams(const MeshQualityOperation *op, std::set<MeshQualityOperation*> &ops) {
    for(auto *op_d : op->next) {
        if(std::find(ops.begin(), ops.end(), op_d) == ops.end()) {
            ops.insert(op_d);
            MeshQuality_downstreams(op_d, ops);
        }
    }
}

std::set<MeshQualityOperation*> MeshQualityOperation::downstreams() const {
    std::set<MeshQualityOperation*> result;
    MeshQuality_downstreams(this, result);
    return result;
}

template <typename T> 
std::set<MeshQualityOperation*> MeshQualityOperation_headOperations(const T &ops) {
    std::set<MeshQualityOperation*> result;
    for(auto &op : ops) 
        if(op->prev.empty()) 
            result.insert(op);
    return result;
}

std::set<MeshQualityOperation*> MeshQualityOperation::headOperations() const { 
    std::set<MeshQualityOperation*> us = upstreams();
    us.insert(const_cast<MeshQualityOperation*>(this));
    return MeshQualityOperation_headOperations(std::vector<MeshQualityOperation*>(us.begin(), us.end()));
}


////////////////
// Operations //
////////////////


/** Merges two vertices */
struct VertexMergeOperation : MeshQualityOperation {

    int vId, v2Id;
    Vertex *v, *v2;
    std::unordered_set<int> affectedChildren;

    VertexMergeOperation(Mesh *_mesh, Vertex *_source, Vertex *_target) : MeshQualityOperation(_mesh) {
        flags = MeshQualityOperation::Flag::Active;
        vId = _source->objectId();
        v2Id = _target->objectId();
        for(auto &s : _source->sharedSurfaces(_target)) 
            targets.push_back(s->objectId());
    };

    void prep() override {
        v =  mesh->getVertex(vId);
        v2 = mesh->getVertex(v2Id);
        for(auto &t : targets) 
            for(auto &b : mesh->getSurface(t)->getBodies()) 
                affectedChildren.insert(b->objectId());
    }

    std::vector<int> implement() override { 
        MeshSolver::engineLock();

        HRESULT res = v->merge(v2);
        
        MeshSolver::engineUnlock();

        if(res == S_OK) {
            next.clear();
            return std::vector<int>(affectedChildren.begin(), affectedChildren.end());
        }
        return {};
    }
};


/** Creates and inserts a vertex between two vertices */
struct VertexCreateInsertOperation : MeshQualityOperation {

    int v1Id, v2Id;
    Vertex *v1, *v2;
    std::unordered_set<int> affectedChildren;

    VertexCreateInsertOperation(Mesh *_mesh, Vertex *_source, Vertex *_target) : MeshQualityOperation(_mesh) {
        flags = MeshQualityOperation::Flag::Active;
        v1Id = _source->objectId();
        v2Id = _target->objectId();
        for(auto &s : _source->sharedSurfaces(_target)) 
            targets.push_back(s->objectId());
    };

    size_t numNewVertices() const override { return 1; }

    void prep() override {
        v1 = mesh->getVertex(v1Id);
        v2 = mesh->getVertex(v2Id);
        for(auto &t : targets) 
            for(auto &b : mesh->getSurface(t)->getBodies()) 
                affectedChildren.insert(b->objectId());
    }

    std::vector<int> implement() override { 

        MeshSolver::engineLock();
        
        HRESULT res = Vertex::insert((v1->getPosition() + v2->getPosition()) * 0.5, v1, v2) != NULL ? S_OK : E_FAIL;
        
        MeshSolver::engineUnlock();

        if(res == S_OK) {
            next.clear();
            return std::vector<int>(affectedChildren.begin(), affectedChildren.end());
        }
        return {};
    }
};

/** Inserts a vertex between two vertices */
struct VertexInsertOperation : MeshQualityOperation {

    int vId, vaId, vbId;
    Vertex *v, *va, *vb;
    std::unordered_set<int> affectedChildren;

    VertexInsertOperation(Mesh *_mesh, Vertex *_source, Surface *_target, Vertex *_va, Vertex *_vb) : MeshQualityOperation(_mesh) {
        flags = MeshQualityOperation::Flag::Active;
        vId = _source->objectId();
        vaId = _va->objectId();
        vbId = _vb->objectId();

        std::unordered_set<int> target_surfs = {_target->objectId()};
        for(auto &s : _source->getSurfaces()) 
            target_surfs.insert(s->objectId());
        for(auto &s : _target->connectedSurfaces({_va, _vb})) 
            target_surfs.insert(s->objectId());
        targets = std::vector<int>(target_surfs.begin(), target_surfs.end());
    };

    void prep() override {
        v = mesh->getVertex(vId);
        va = mesh->getVertex(vaId);
        vb = mesh->getVertex(vbId);
        for(auto &t : targets) 
            for(auto &b : mesh->getSurface(t)->getBodies()) 
                affectedChildren.insert(b->objectId());
    }

    std::vector<int> implement() override {
        MeshSolver::engineLock();

        HRESULT res = v->insert(va, vb);

        MeshSolver::engineUnlock();

        if(res == S_OK) {
            next.clear();
            return std::vector<int>(affectedChildren.begin(), affectedChildren.end());
        }
        return {};
    }
};

/** Converts a body to a vertex */
struct BodyDemoteOperation : MeshQualityOperation {

    int toReplaceId;
    Body *toReplace;

    BodyDemoteOperation(Mesh *_mesh, Body *_source) : MeshQualityOperation(_mesh) {
        flags = MeshQualityOperation::Flag::Active;
        toReplaceId = _source->objectId();
        for(auto &nb : adjacentTo(removedBodiesByB2V(_source))) 
            targets.push_back(nb->objectId());
    }

    size_t numNewVertices() const override { return 1; };

    void prep() override {
        toReplace = mesh->getBody(toReplaceId);
    }

    std::vector<int> implement() override {
        const FVector3 toReplaceCentroid = toReplace->getCentroid();

        MeshSolver::engineLock();

        HRESULT res = Vertex::replace(toReplaceCentroid, toReplace) != NULL ? S_OK : E_FAIL;

        MeshSolver::engineUnlock();

        if(res == S_OK) 
            next.clear();
        return {};
    };
};

/** Converts a surface to a vertex */
struct SurfaceDemoteOperation : MeshQualityOperation {

    int toReplaceId;
    Surface *toReplace;
    std::unordered_set<int> affectedChildren;

    SurfaceDemoteOperation(Mesh *_mesh, Surface *_source) : MeshQualityOperation(_mesh) {
        flags = MeshQualityOperation::Flag::Active;
        toReplaceId = _source->objectId();
        std::unordered_set<Surface*> removed = removedSurfacesByS2V(_source);
        for(auto &nb : removed) 
            targets.push_back(nb->objectId());
        for(auto &nb : connectedSurfacesToS2V(removed)) 
            if(nb->objectId() != _source->objectId()) 
                targets.push_back(nb->objectId());
    }

    size_t numNewVertices() const override { return 1; };

    void prep() override {
        toReplace = mesh->getSurface(toReplaceId);
        for(auto &t : targets) 
            for(auto &b : mesh->getSurface(t)->getBodies()) 
                affectedChildren.insert(b->objectId());
    }

    std::vector<int> implement() override {
        const FVector3 toReplaceCentroid = toReplace->getCentroid();

        MeshSolver::engineLock();
        
        HRESULT res = Vertex::replace(toReplaceCentroid, toReplace) != NULL ? S_OK : E_FAIL;

        MeshSolver::engineUnlock();

        if(res == S_OK) {
            next.clear();
            return std::vector<int>(affectedChildren.begin(), affectedChildren.end());
        }
        return {};
    }
};

/** Converts a surface edge to a vertex */
struct EdgeDemoteOperation : MeshQualityOperation {

    int vId, v2Id;
    Vertex *v, *v2;
    std::unordered_set<int> affectedChildren;

    EdgeDemoteOperation(Mesh *_mesh, Vertex *_v1, Vertex *_v2) : MeshQualityOperation(_mesh) {
        flags = MeshQualityOperation::Flag::Active;
        vId = _v1->objectId();
        v2Id = _v2->objectId();

        std::unordered_set<int> targets_set{vId, v2Id};
        for(auto &c : _v1->connectedVertices()) 
            targets_set.insert(c->objectId());
        for(auto &c : _v2->connectedVertices()) 
            targets_set.insert(c->objectId());
        targets = std::vector<int>(targets_set.begin(), targets_set.end());
    }

    void prep() override {
        v = mesh->getVertex(vId);
        v2 = mesh->getVertex(v2Id);
        for(auto &t : targets) 
            for(auto &s : mesh->getVertex(t)->getSurfaces()) 
                affectedChildren.insert(s->objectId());
    }

    std::vector<int> implement() override {
        MeshSolver::engineLock();
        HRESULT res = v->merge(v2);
        MeshSolver::engineUnlock();
        if(res == S_OK) {
            next.clear();
            return std::vector<int>(affectedChildren.begin(), affectedChildren.end());
        }
        return {};
    }
};

/** Converts a vertex to an edge */
struct EdgeInsertOperation : MeshQualityOperation {

    int v0Id, v1Id, v2Id;
    Vertex *v0, *v1, *v2;
    std::unordered_set<int> affectedChildren;

    EdgeInsertOperation(Mesh *_mesh, Vertex *_source, Vertex *_target1, Vertex *_target2) : MeshQualityOperation(_mesh) {
        flags = MeshQualityOperation::Flag::Active;
        v0Id = _source->objectId();
        v1Id = _target1->objectId();
        v2Id = _target2->objectId();
        std::unordered_set<int> _targets;
        for(auto &s : _source->sharedSurfaces(_target1)) 
            _targets.insert(s->objectId());
        for(auto &s : _source->sharedSurfaces(_target2)) 
            _targets.insert(s->objectId());
        targets = std::vector<int>(_targets.begin(), _targets.end());
    }

    void prep() override {
        v0 = mesh->getVertex(v0Id);
        v1 = mesh->getVertex(v1Id);
        v2 = mesh->getVertex(v2Id);
        for(auto &t : targets) 
            for(auto &s : mesh->getVertex(t)->getSurfaces()) 
                affectedChildren.insert(s->objectId());
    }

    size_t numNewVertices() const override { return 2; };

    std::vector<int> implement() override {
        const FVector3 pos0 = v0->getPosition();
        const FVector3 pos1 = (pos0 + v1->getPosition()) * 0.5;
        const FVector3 pos2 = (pos0 + v2->getPosition()) * 0.5;

        MeshSolver::engineLock();

        Vertex::insert(pos1, v0, v1);
        Vertex::insert(pos2, v0, v2);
        v0->destroy();

        MeshSolver::engineUnlock();

        next.clear();
        return std::vector<int>(affectedChildren.begin(), affectedChildren.end());
    }

};


static bool MeshQuality_vertexSplitTest(
    Vertex *v, 
    Mesh *m, 
    const FloatP_t &edgeSplitDist, 
    FVector3 &sep, 
    const std::vector<Vertex*> &v_nbs, 
    std::vector<Vertex*> &vert_nbs, 
    std::vector<Vertex*> &new_vert_nbs
) {
    Particle *p = v->particle()->part();

    // Calculate current relative force
    FVector3 force_rel(0);
    for(auto &vn : v_nbs) {
        force_rel += vn->particle()->getForce();
    }
    force_rel -= p->force * v_nbs.size();
    FPTYPE mask[] = {
        (p->flags & PARTICLE_FROZEN_X) ? 0.0f : 1.0f,
        (p->flags & PARTICLE_FROZEN_Y) ? 0.0f : 1.0f,
        (p->flags & PARTICLE_FROZEN_Z) ? 0.0f : 1.0f
    };
    for(int k = 0; k < 3; k++) force_rel[k] *= mask[k];
    if(force_rel.isZero()) 
        return false;

    sep = force_rel.normalized() * edgeSplitDist;
    if(sep.isZero()) 
        return false;

    // Get split plan along direction of force
    v->splitPlan(sep, vert_nbs, new_vert_nbs);

    // Enforce that a split must create a new shared edge
    if(vert_nbs.size() < 2 || new_vert_nbs.size() < 2) 
        return false;

    // Calculate relative force on each vertex of a new edge
    FVector3 vert_force_rel, new_vert_force_rel;
    for(auto &vn : vert_nbs) 
        vert_force_rel += vn->particle()->getForce();
    for(auto &vn : new_vert_nbs) 
        new_vert_force_rel += vn->particle()->getForce();
    vert_force_rel -= p->force * 0.5 * vert_nbs.size();
    new_vert_force_rel -= p->force * 0.5 * new_vert_nbs.size();

    // Test whether the new edge would be in tension and return true if so
    return sep.dot(vert_force_rel) < 0 && sep.dot(new_vert_force_rel) > 0;
}

/** Splits a vertex into an edge */
struct VertexSplitOperation : MeshQualityOperation { 

    int vId;
    FVector3 sep;
    Vertex *v;
    std::vector<int> vert_nbsIds, new_vert_nbsIds;
    std::vector<Vertex*> vert_nbs, new_vert_nbs;
    std::unordered_set<int> affectedChildren;

    VertexSplitOperation(Mesh *_mesh, Vertex *_source, const FVector3 &_sep, std::vector<Vertex*> _vert_nbs, std::vector<Vertex*> _new_vert_nbs) : 
        MeshQualityOperation(_mesh), 
        sep{_sep}
    {
        flags = MeshQualityOperation::Flag::Active;
        vId = _source->objectId();
        for(auto &_v : _vert_nbs) {
            vert_nbsIds.push_back(_v->objectId());
            targets.push_back(_v->objectId());
        }
        for(auto &_v : _new_vert_nbs) {
            new_vert_nbsIds.push_back(_v->objectId());
            targets.push_back(_v->objectId());
        }
    }

    size_t numNewVertices() const override { return 1; };

    void prep() override {
        v = mesh->getVertex(vId);
        vert_nbs.reserve(vert_nbsIds.size());
        new_vert_nbs.reserve(new_vert_nbsIds.size());
        for(auto &vn : vert_nbsIds) 
            vert_nbs.push_back(mesh->getVertex(vn));
        for(auto &vn : new_vert_nbsIds) 
            new_vert_nbs.push_back(mesh->getVertex(vn));
        for(auto &t : targets) 
            for(auto &s : mesh->getVertex(t)->getSurfaces()) 
                affectedChildren.insert(s->objectId());
    }

    std::vector<int> implement() override { 

        // Create a candidate vertex
        MeshSolver::engineLock();
        Vertex *new_v = v->splitExecute(sep, vert_nbs, new_vert_nbs);
        MeshSolver::engineUnlock();

        // Only invalidate if a vertex was created, since some requested configurations are invalid and subsequently ignored
        if(new_v) {
            next.clear();
            return std::vector<int>(affectedChildren.begin(), affectedChildren.end());
        }
        return {};
    }
};


////////////////////////////////////////////////////////////////////////////////////////
// RNR (Okuda I<->H) topology walk + Condition-4 guards -- Phase B                       //
////////////////////////////////////////////////////////////////////////////////////////
//
// The native 3D T1 / reversible network reconnection (RNR), i.e. the Okuda I<->H
// face<->edge swap (Okuda et al. 2013, Biomech Model Mechanobiol 12:627-644). This is the
// topology operation TissueForge otherwise lacks; the other 3D quality ops are only
// degenerate collapses (BodyDemote/SurfaceDemote/EdgeDemote).
//
// This block is the READ-ONLY half: it gathers + validates the [I]/[H] neighborhood of a
// candidate reconnection site (Okuda Fig. 3) and vetoes the topologically-irreversible
// patterns (Okuda Condition 4, Figs. 6/9). It is a faithful translation of our validated
// Python prototype rnr/topology.py (the neighborhood walk) and rnr/conditions.py (the
// Condition-4 guards) and the Okuda 2013 equations -- NOT copied from the GPL tvm/3DVertVor
// reference (license boundary, CLAUDE.md). The MUTATE half (Okuda Appendix-1 vertex
// placement + manual surface surgery) lives in ReconnectionOperation::implement() (Phase C,
// done); Phase D wires the scan below into doQuality's live quality loop as the active pass.
//
// TissueForge has no explicit Edge object: an "edge" is an ordered consecutive vertex pair
// in a surface's vertex ring. The whole walk uses only TF adjacency helpers
// (Vertex::getBodies/connectedVertices/sharedSurfaces, Surface::getBodies/getVertices/
// neighborVertices/numSharedContiguousVertexSets, Body::findInterface).

namespace {

// --- small adjacency helpers (read-only); mirror topology.py -----------------

/** Stable set of body ids touching a vertex (objectId is the only stable identity). */
static std::set<int> rnr_bodyIds(Vertex *v) {
    std::set<int> out;
    for(auto *b : v->getBodies()) out.insert(b->objectId());
    return out;
}

/**
 * Edge length under the mesh-level periodic geometry convention. The Python oracle uses a
 * plain coordinate difference in finite-cluster tests; when periodic geometry is enabled this
 * becomes the minimum-image length needed for bulk periodic meshes.
 */
static FloatP_t rnr_edgeLength(Vertex *a, Vertex *b) {
    return meshRelativePosition(a->getPosition(), b->getPosition()).length();
}

/** True iff va, vb are consecutive (cyclically) in `s`'s vertex ring (topology.is_consecutive). */
static bool rnr_isConsecutive(Surface *s, Vertex *va, Vertex *vb) {
    std::tuple<Vertex*, Vertex*> nb = s->neighborVertices(va);
    Vertex *vp = std::get<0>(nb), *vn = std::get<1>(nb);
    const int bId = vb->objectId();
    return (vp && vp->objectId() == bId) || (vn && vn->objectId() == bId);
}

/**
 * The ring-neighbor of v in `s` other than excludeId (null unless exactly one remains);
 * mirrors topology.other_neighbor. A polygon vertex has exactly two ring-neighbors; one is
 * the short-edge partner (excluded), the other is the surface's outer vertex on v's side.
 */
static Vertex *rnr_otherNeighbor(Surface *s, Vertex *v, int excludeId) {
    std::tuple<Vertex*, Vertex*> nb = s->neighborVertices(v);
    Vertex *vp = std::get<0>(nb), *vn = std::get<1>(nb);
    Vertex *found = nullptr;
    int n = 0;
    if(vp && vp->objectId() != excludeId) { found = vp; n++; }
    if(vn && vn->objectId() != excludeId) { found = vn; n++; }
    return n == 1 ? found : nullptr;
}

// --- Condition-4 primitives (conditions.py) ----------------------------------

/** 4(iii), the tvm "extra rule": two cells share >= 2 faces. */
static bool rnr_cellsShareMultipleFaces(Body *a, Body *b) { return a->findInterface(b).size() >= 2; }

/** 4(ii): two faces share >= 2 separate edges (>= 2 maximal runs of contiguous shared verts). */
static bool rnr_facesShareMultipleEdges(Surface *a, Surface *b) { return a->numSharedContiguousVertexSets(b) >= 2; }

// --- the [I] neighborhood (edge state); mirrors topology.i_neighbourhood ------

/** One of the 3 arms of the [I] neighborhood, keyed by a SIDE surface bearing the short edge. */
struct RNR_Arm {
    Surface *sideSurface = nullptr;
    int sideSurfaceId = -1;
    Vertex *outerTop = nullptr;   // outer vertex on the v10 side (within this side surface)
    Vertex *outerBot = nullptr;   // outer vertex on the v11 side
};

/** Full [I] neighborhood of a short edge (v10, v11). Build-time handles; the operation that
 *  outlives a mutation persists only ids and re-walks by id (handles invalidate on mutate). */
struct RNR_IConfig {
    bool valid = false;
    Vertex *v10 = nullptr, *v11 = nullptr;
    int v10Id = -1, v11Id = -1;
    std::vector<int> sideCellIds;            // the 3 side cells (sorted)
    Body *capTop = nullptr, *capBot = nullptr;
    int capTopId = -1, capBotId = -1;
    std::vector<RNR_Arm> arms;               // 3 arms (one per side surface)
    std::map<int, Surface*> topFaces;        // sideCellId -> interface(side, capTop)
    std::map<int, Surface*> bottomFaces;     // sideCellId -> interface(side, capBot)
    FloatP_t length = 0;
};

/**
 * Gather + validate the [I] neighborhood of the ordered short edge (v10, v11). Returns false
 * unless (v10, v11) is a clean interior 3-cell edge in the canonical Okuda [I] config (Fig. 3):
 * both endpoints interior (4 cells), edge shared by exactly the 3 side cells with one cap each,
 * consecutive in exactly the 3 (2-body) side surfaces, and the 3 top + 3 bottom side<->cap
 * interfaces unique. Purely structural; the Condition-2 threshold + Condition-4 vetoes are
 * applied separately.
 */
static bool rnr_iNeighbourhood(Vertex *v10, Vertex *v11, RNR_IConfig &cfg) {
    // both endpoints interior: shared by exactly 4 cells.
    std::set<int> b10 = rnr_bodyIds(v10), b11 = rnr_bodyIds(v11);
    if(b10.size() != 4 || b11.size() != 4) return false;

    // edge shared by exactly the 3 SIDE cells (intersection); each endpoint one unique CAP cell.
    std::set<int> side, capTopIds, capBotIds;
    for(int id : b10) { if(b11.count(id)) side.insert(id); else capTopIds.insert(id); }
    for(int id : b11) { if(!b10.count(id)) capBotIds.insert(id); }
    if(side.size() != 3 || capTopIds.size() != 1 || capBotIds.size() != 1) return false;

    // short edge consecutive in exactly 3 surfaces, all interior (2-body).
    std::vector<Surface*> shared = v10->sharedSurfaces(v11);
    std::vector<Surface*> sideSurfs;
    for(auto *s : shared) if(rnr_isConsecutive(s, v10, v11)) sideSurfs.push_back(s);
    if(sideSurfs.size() != 3) return false;
    for(auto *s : sideSurfs) if(s->getBodies().size() != 2) return false;

    const int capTopId = *capTopIds.begin(), capBotId = *capBotIds.begin();
    Body *capTop = nullptr, *capBot = nullptr;
    for(auto *b : v10->getBodies()) if(b->objectId() == capTopId) { capTop = b; break; }
    for(auto *b : v11->getBodies()) if(b->objectId() == capBotId) { capBot = b; break; }
    if(!capTop || !capBot) return false;

    // body-id -> handle for the side cells (the side surfaces between them cover all 3).
    std::map<int, Body*> bid;
    for(auto *s : sideSurfs) for(auto *b : s->getBodies()) bid[b->objectId()] = b;

    // build 3 arms: each side surface -> its two outer vertices (toward each cap).
    cfg.arms.clear();
    for(auto *s : sideSurfs) {
        Vertex *ot = rnr_otherNeighbor(s, v10, v11->objectId());
        Vertex *ob = rnr_otherNeighbor(s, v11, v10->objectId());
        if(!ot || !ob) return false;
        RNR_Arm a;
        a.sideSurface = s; a.sideSurfaceId = s->objectId();
        a.outerTop = ot; a.outerBot = ob;
        cfg.arms.push_back(a);
    }

    // 3 top + 3 bottom faces (side<->cap), each a unique surface (the 4(iii) extra rule).
    std::map<int, Surface*> topFaces, bottomFaces;
    for(int scId : side) {
        Body *sc = bid[scId];
        std::vector<Surface*> it = sc->findInterface(capTop);
        std::vector<Surface*> ib = sc->findInterface(capBot);
        if(it.size() != 1 || ib.size() != 1) return false;
        topFaces[scId] = it[0];
        bottomFaces[scId] = ib[0];
    }

    cfg.valid = true;
    cfg.v10 = v10; cfg.v11 = v11; cfg.v10Id = v10->objectId(); cfg.v11Id = v11->objectId();
    cfg.sideCellIds.assign(side.begin(), side.end());   // std::set iterates ascending
    cfg.capTop = capTop; cfg.capBot = capBot; cfg.capTopId = capTopId; cfg.capBotId = capBotId;
    cfg.topFaces = topFaces; cfg.bottomFaces = bottomFaces;
    cfg.length = rnr_edgeLength(v10, v11);
    return true;
}

/** I->H veto (conditions.i_to_h_veto): "" if legal, else the Okuda-irreversible reason. */
static std::string rnr_iToHVeto(const RNR_IConfig &cfg) {
    // [beta]: the two caps must not already be in contact (else the new triangle is a 2nd face).
    if(cfg.capTop->findInterface(cfg.capBot).size() != 0)
        return "caps c123,c456 already share a face (would create double trigonal face, 4(iii)/[beta])";

    // [alpha]: reconnecting an edge of a trigonal (3-vertex) side face creates a double edge.
    for(auto &a : cfg.arms)
        if(a.sideSurface->getVertices().size() == 3)
            return "side surface " + std::to_string(a.sideSurfaceId) +
                   " is a triangle (reconnecting its edge -> double edge, 4(i)/[alpha])";

    // 4(iii) among side-cell pairs (each pair should share exactly the one side face).
    std::map<int, Body*> bid;
    for(auto &a : cfg.arms) for(auto *b : a.sideSurface->getBodies()) bid[b->objectId()] = b;
    const std::vector<int> &side = cfg.sideCellIds;
    for(size_t i = 0; i < side.size(); i++)
        for(size_t j = i + 1; j < side.size(); j++)
            if(rnr_cellsShareMultipleFaces(bid[side[i]], bid[side[j]]))
                return "side cells " + std::to_string(side[i]) + "," + std::to_string(side[j]) +
                       " already share >=2 faces (4(iii))";

    // 4(iii) side cell vs cap (each side<->cap interface must be a single face).
    for(int scId : side) {
        Body *sc = bid[scId];
        if(rnr_cellsShareMultipleFaces(sc, cfg.capTop) || rnr_cellsShareMultipleFaces(sc, cfg.capBot))
            return "side cell " + std::to_string(scId) + " already shares >=2 faces with a cap (4(iii))";
    }

    return "";
}

// --- the [H] neighborhood (triangle state); mirrors topology.h_neighbourhood --

/** One arm of the [H] neighborhood, keyed by a triangle VERTEX (v7/v8/v9). */
struct RNR_HArm {
    Vertex *triVertex = nullptr;
    int triVertexId = -1;
    Surface *sideSurface = nullptr;          // the side face holding this triangle vertex
    int sideSurfaceId = -1;
    Vertex *outerTop = nullptr;              // outer vertex toward cap_top
    Vertex *outerBot = nullptr;              // outer vertex toward cap_bot
};

struct RNR_HConfig {
    bool valid = false;
    Surface *triangle = nullptr;
    int triangleId = -1;
    std::vector<int> triVertexIds;
    Body *capTop = nullptr, *capBot = nullptr;
    int capTopId = -1, capBotId = -1;
    std::vector<int> sideCellIds;            // the 3 side cells (sorted)
    std::vector<RNR_HArm> arms;              // 3 arms (one per triangle vertex / side face)
    std::map<int, Surface*> topFaces;        // sideCellId -> interface(side, capTop)
    std::map<int, Surface*> bottomFaces;     // sideCellId -> interface(side, capBot)
    FloatP_t maxEdge = 0;
};

/**
 * Gather + validate the [H] neighborhood of a triangular face (the [H] state). Mirrors
 * i_neighbourhood so the reverse H->I collapse gets the same handles the forward I->H build
 * used. A reconnection-site triangle has exactly 3 vertices, is shared by exactly 2 cells (the
 * caps), each triangle edge is additionally shared by exactly one side cell, each triangle
 * vertex sits in exactly one side face, and the surrounding 9 faces resolve uniquely.
 */
static bool rnr_hNeighbourhood(Surface *tri, RNR_HConfig &cfg) {
    std::vector<Vertex*> verts = tri->getVertices();
    if(verts.size() != 3) return false;
    std::vector<Body*> caps = tri->getBodies();
    if(caps.size() != 2) return false;
    Body *capTop = caps[0], *capBot = caps[1];
    for(auto *v : verts) if(rnr_bodyIds(v).size() != 4) return false;

    const std::vector<int> vids = {verts[0]->objectId(), verts[1]->objectId(), verts[2]->objectId()};
    std::map<int, Vertex*> vbyid;
    for(auto *v : verts) vbyid[v->objectId()] = v;
    const int capTopId = capTop->objectId(), capBotId = capBot->objectId();

    // the three triangle edges, as ordered pairs of vertex ids.
    const int E[3][2] = {{vids[0], vids[1]}, {vids[1], vids[2]}, {vids[2], vids[0]}};

    // side cell per triangle edge = common cell of the edge's endpoints minus the caps.
    int edgeSideCell[3];
    std::set<int> sideIds;
    for(int e = 0; e < 3; e++) {
        std::set<int> ba = rnr_bodyIds(vbyid[E[e][0]]);
        std::set<int> bb = rnr_bodyIds(vbyid[E[e][1]]);
        std::vector<int> rest;
        for(int id : ba) if(bb.count(id) && id != capTopId && id != capBotId) rest.push_back(id);
        if(rest.size() != 1) return false;
        edgeSideCell[e] = rest[0];
        sideIds.insert(rest[0]);
    }
    if(sideIds.size() != 3) return false;

    // body-id -> handle for caps + side cells.
    std::map<int, Body*> bid;
    bid[capTopId] = capTop; bid[capBotId] = capBot;
    for(auto *v : verts) for(auto *b : v->getBodies()) bid[b->objectId()] = b;

    FloatP_t maxEdge = 0;
    for(int e = 0; e < 3; e++) maxEdge = std::max(maxEdge, rnr_edgeLength(vbyid[E[e][0]], vbyid[E[e][1]]));

    // each triangle vertex sits in exactly one SIDE face: the interface between the two side
    // cells flanking it (the side cells of its two incident triangle edges).
    cfg.arms.clear();
    for(auto *v : verts) {
        const int vId = v->objectId();
        std::vector<int> flank;
        for(int e = 0; e < 3; e++) if(E[e][0] == vId || E[e][1] == vId) flank.push_back(edgeSideCell[e]);
        if(flank.size() != 2) return false;
        Body *scA = bid[flank[0]], *scB = bid[flank[1]];
        std::vector<Surface*> iface = scA->findInterface(scB);
        Surface *sideSurf = nullptr;
        for(auto *s : iface) {
            for(auto *w : s->getVertices()) if(w->objectId() == vId) { sideSurf = s; break; }
            if(sideSurf) break;
        }
        if(!sideSurf) return false;
        // the two ring-neighbors of v in that side face are its outer vertices (one per cap).
        std::tuple<Vertex*, Vertex*> nb = sideSurf->neighborVertices(v);
        Vertex *vp = std::get<0>(nb), *vn = std::get<1>(nb);
        if(!vp || !vn) return false;
        Vertex *outerTop = nullptr, *outerBot = nullptr;
        for(Vertex *w : {vp, vn}) {
            std::set<int> wb = rnr_bodyIds(w);
            if(wb.count(capTopId)) outerTop = w;
            else if(wb.count(capBotId)) outerBot = w;
        }
        if(!outerTop || !outerBot) return false;
        RNR_HArm a;
        a.triVertex = v; a.triVertexId = vId;
        a.sideSurface = sideSurf; a.sideSurfaceId = sideSurf->objectId();
        a.outerTop = outerTop; a.outerBot = outerBot;
        cfg.arms.push_back(a);
    }

    // 3 top + 3 bottom faces (side<->cap), each a unique surface.
    std::map<int, Surface*> topFaces, bottomFaces;
    for(int scId : sideIds) {
        Body *sc = bid[scId];
        std::vector<Surface*> it = sc->findInterface(capTop);
        std::vector<Surface*> ib = sc->findInterface(capBot);
        if(it.size() != 1 || ib.size() != 1) return false;
        topFaces[scId] = it[0];
        bottomFaces[scId] = ib[0];
    }

    cfg.valid = true;
    cfg.triangle = tri; cfg.triangleId = tri->objectId();
    cfg.triVertexIds = vids;
    cfg.capTop = capTop; cfg.capBot = capBot; cfg.capTopId = capTopId; cfg.capBotId = capBotId;
    cfg.sideCellIds.assign(sideIds.begin(), sideIds.end());
    cfg.topFaces = topFaces; cfg.bottomFaces = bottomFaces;
    cfg.maxEdge = maxEdge;
    return true;
}

/** H->I veto (conditions.h_to_i_veto): "" if legal, else the Okuda-irreversible reason. */
static std::string rnr_hToIVeto(const RNR_HConfig &cfg) {
    if(rnr_cellsShareMultipleFaces(cfg.capTop, cfg.capBot))
        return "caps share >=2 faces (would leave a second contact, 4(iii)/[beta])";

    std::map<int, Body*> bid;
    bid[cfg.capTopId] = cfg.capTop; bid[cfg.capBotId] = cfg.capBot;
    for(auto *v : cfg.triangle->getVertices()) for(auto *b : v->getBodies()) bid[b->objectId()] = b;

    const std::vector<int> &side = cfg.sideCellIds;
    for(size_t i = 0; i < side.size(); i++)
        for(size_t j = i + 1; j < side.size(); j++)
            if(rnr_cellsShareMultipleFaces(bid[side[i]], bid[side[j]]))
                return "side cells " + std::to_string(side[i]) + "," + std::to_string(side[j]) +
                       " share >=2 faces (4(iii))";

    for(size_t i = 0; i < cfg.arms.size(); i++)
        for(size_t j = i + 1; j < cfg.arms.size(); j++)
            if(rnr_facesShareMultipleEdges(cfg.arms[i].sideSurface, cfg.arms[j].sideSurface))
                return "side faces " + std::to_string(cfg.arms[i].sideSurfaceId) + "," +
                       std::to_string(cfg.arms[j].sideSurfaceId) +
                       " already share >=2 edges (collapse would double an edge, 4(ii)/[alpha])";

    return "";
}

// --- scanners (Condition-2 triggers); mirror topology.find_short_edges/find_small_triangles --

/** All interior short edges (length < threshold) that form a valid [I] config. Each undirected
 *  edge is visited once (from its smaller-id endpoint). */
static void rnr_findShortEdges(Mesh *mesh, FloatP_t threshold, std::vector<RNR_IConfig> &out) {
    const unsigned int n = mesh->sizeVertices();
    for(unsigned int i = 0; i < n; i++) {
        Vertex *v = mesh->getVertex(i);
        if(!v || v->getBodies().size() != 4) continue;
        for(auto *w : v->connectedVertices()) {
            if(w->objectId() <= v->objectId()) continue;
            if(rnr_edgeLength(v, w) >= threshold) continue;
            RNR_IConfig cfg;
            if(rnr_iNeighbourhood(v, w, cfg)) out.push_back(cfg);
        }
    }
}

/** All triangular faces whose MAX edge < threshold (Okuda Condition 2 for the reverse direction
 *  -- max, NOT min) that form a valid [H] config. */
static void rnr_findSmallTriangles(Mesh *mesh, FloatP_t threshold, std::vector<RNR_HConfig> &out) {
    const unsigned int n = mesh->sizeSurfaces();
    for(unsigned int i = 0; i < n; i++) {
        Surface *s = mesh->getSurface(i);
        if(!s || s->getVertices().size() != 3) continue;
        RNR_HConfig cfg;
        if(rnr_hNeighbourhood(s, cfg) && cfg.maxEdge < threshold) out.push_back(cfg);
    }
}

/** The deterministic surface slot an I->H op is keyed at: the smallest of its 3 side-surface ids
 *  (a surface the op touches, so the dependency graph serializes conflicts correctly). */
static int rnr_iTriggerSlot(const RNR_IConfig &cfg) {
    int slot = -1;
    for(auto &a : cfg.arms) if(slot < 0 || a.sideSurfaceId < slot) slot = a.sideSurfaceId;
    return slot;
}

// --- JSON emit (diagnostic; debug entry points cross-checked vs the Python oracle) -----------

static std::string rnr_jsonEscape(const std::string &s) {
    std::string out;
    for(char c : s) {
        if(c == '"' || c == '\\') { out.push_back('\\'); out.push_back(c); }
        else out.push_back(c);
    }
    return out;
}

static std::string rnr_jsonIntArray(const std::vector<int> &v) {
    std::stringstream ss;
    ss << "[";
    for(size_t i = 0; i < v.size(); i++) { if(i) ss << ", "; ss << v[i]; }
    ss << "]";
    return ss.str();
}

static std::string rnr_iConfigJson(const RNR_IConfig &cfg, const std::string &veto) {
    std::stringstream ss;
    ss << "{\"valid\": true, \"kind\": \"I\""
       << ", \"v10_id\": " << cfg.v10Id
       << ", \"v11_id\": " << cfg.v11Id
       << ", \"triangle_id\": -1"
       << ", \"trigger_surface_id\": " << rnr_iTriggerSlot(cfg)
       << ", \"cap_top_id\": " << cfg.capTopId
       << ", \"cap_bot_id\": " << cfg.capBotId
       << ", \"side_cell_ids\": " << rnr_jsonIntArray(cfg.sideCellIds)
       << ", \"length\": " << cfg.length
       << ", \"legal\": " << (veto.empty() ? "true" : "false")
       << ", \"veto_reason\": \"" << rnr_jsonEscape(veto) << "\"}";
    return ss.str();
}

static std::string rnr_hConfigJson(const RNR_HConfig &cfg, const std::string &veto) {
    std::stringstream ss;
    ss << "{\"valid\": true, \"kind\": \"H\""
       << ", \"v10_id\": -1, \"v11_id\": -1"
       << ", \"triangle_id\": " << cfg.triangleId
       << ", \"trigger_surface_id\": " << cfg.triangleId
       << ", \"cap_top_id\": " << cfg.capTopId
       << ", \"cap_bot_id\": " << cfg.capBotId
       << ", \"side_cell_ids\": " << rnr_jsonIntArray(cfg.sideCellIds)
       << ", \"length\": " << cfg.maxEdge
       << ", \"legal\": " << (veto.empty() ? "true" : "false")
       << ", \"veto_reason\": \"" << rnr_jsonEscape(veto) << "\"}";
    return ss.str();
}

static std::string rnr_reconnectResultJson(
    bool ok,
    const std::string &reason,
    int newSurfaceId,
    const std::vector<int> &newVertexIds
) {
    std::stringstream ss;
    ss << "{\"ok\": " << (ok ? "true" : "false")
       << ", \"reason\": \"" << rnr_jsonEscape(reason) << "\""
       << ", \"new_surface_id\": " << newSurfaceId
       << ", \"new_vertex_ids\": " << rnr_jsonIntArray(newVertexIds)
       << "}";
    return ss.str();
}

// --- Phase-C placement math (Okuda 2013 Appendix 1); mirrors reconnect.py ----

static FVector3 rnr_unit(const FVector3 &v) {
    return v.isZero() ? v : v.normalized();
}

/** I->H vertex placement, Okuda Eqs. 46-56. */
static bool rnr_placeIToH(const RNR_IConfig &cfg, const FloatP_t &dlTh, std::array<FVector3, 3> &out) {
    if(cfg.arms.size() != 3) return false;

    const FVector3 p10 = cfg.v10->getPosition();
    const FVector3 p11 = meshPositionNear(cfg.v11->getPosition(), p10);
    const FVector3 r0 = p10 + (p11 - p10) * 0.5;           // Eq. 50: edge midpoint
    const FVector3 uT = rnr_unit(meshRelativePosition(p10, p11)); // Eq. 49: edge axis
    if(uT.isZero()) return false;

    std::array<FVector3, 3> vproj;
    for(size_t i = 0; i < cfg.arms.size(); i++) {
        const RNR_Arm &a = cfg.arms[i];
        const FVector3 dTop = rnr_unit(meshRelativePosition(a.outerTop->getPosition(), r0));
        const FVector3 dBot = rnr_unit(meshRelativePosition(a.outerBot->getPosition(), r0));
        const FVector3 w = (dTop + dBot) * 0.5;            // Eqs. 54-56
        vproj[i] = w - uT * w.dot(uT);                    // Eqs. 51-53: project off edge
    }

    FloatP_t lMax = 0;
    for(size_t i = 0; i < vproj.size(); i++)
        for(size_t j = i + 1; j < vproj.size(); j++)
            lMax = std::max(lMax, (vproj[i] - vproj[j]).length());
    if(lMax == 0) lMax = 1;

    for(size_t i = 0; i < vproj.size(); i++)
        out[i] = meshWrapPosition(r0 + vproj[i] * (dlTh / lMax)); // Eqs. 46-48
    return true;
}

/** H->I vertex placement, Okuda Eqs. 42-45. */
static bool rnr_placeHToI(const RNR_HConfig &cfg, const FloatP_t &dlTh, FVector3 &p10, FVector3 &p11) {
    if(cfg.arms.size() != 3) return false;

    std::array<FVector3, 3> p = {
        cfg.arms[0].triVertex->getPosition(),
        meshPositionNear(cfg.arms[1].triVertex->getPosition(), cfg.arms[0].triVertex->getPosition()),
        meshPositionNear(cfg.arms[2].triVertex->getPosition(), cfg.arms[0].triVertex->getPosition())
    };
    const FVector3 r0 = (p[0] + p[1] + p[2]) / 3.0;        // Eq. 45: triangle centroid
    FVector3 n = rnr_unit(Magnum::Math::cross(p[1] - p[0], p[2] - p[0])); // Eq. 44
    if(n.isZero()) return false;

    FVector3 topMean(0);
    for(auto &a : cfg.arms)
        topMean += meshPositionNear(a.outerTop->getPosition(), r0);
    topMean /= (FloatP_t)cfg.arms.size();
    if((topMean - r0).dot(n) < 0)
        n = n * -1;

    const FloatP_t half = 0.5 * dlTh;
    p10 = meshWrapPosition(r0 + n * half);                 // Eq. 42
    p11 = meshWrapPosition(r0 - n * half);                 // Eq. 43
    return true;
}

// --- low-level manual surgery (maintain both sides of every adjacency) -------

static bool rnr_vertexHasSurface(Vertex *v, Surface *s) {
    if(!v || !s) return false;
    for(auto *vs : v->getSurfaces())
        if(vs && vs->objectId() == s->objectId())
            return true;
    return false;
}

static bool rnr_surfaceHasBody(Surface *s, Body *b) {
    if(!s || !b) return false;
    for(auto *sb : s->getBodies())
        if(sb && sb->objectId() == b->objectId())
            return true;
    return false;
}

static bool rnr_bodyHasSurface(Body *b, Surface *s) {
    if(!b || !s) return false;
    for(auto *bs : b->getSurfaces())
        if(bs && bs->objectId() == s->objectId())
            return true;
    return false;
}

static HRESULT rnr_replaceV(Surface *s, Vertex *oldV, Vertex *newV) {
    if(!s || !oldV || !newV) return E_FAIL;
    s->replace(newV, oldV);                 // TF edits only the surface ring here.
    if(!rnr_vertexHasSurface(newV, s)) newV->add(s);
    if(rnr_vertexHasSurface(oldV, s)) oldV->remove(s);
    return S_OK;
}

static HRESULT rnr_insertBetween(Surface *s, Vertex *newV, Vertex *v1, Vertex *v2) {
    if(!s || !newV || !v1 || !v2) return E_FAIL;
    if(s->insert(newV, v1, v2) != S_OK) return E_FAIL;
    if(!rnr_vertexHasSurface(newV, s)) newV->add(s);
    return S_OK;
}

static HRESULT rnr_dropV(Surface *s, Vertex *v) {
    if(!s || !v) return E_FAIL;
    if(s->remove(v) != S_OK) return E_FAIL;
    if(rnr_vertexHasSurface(v, s)) v->remove(s);
    return S_OK;
}

static HRESULT rnr_attachBody(Surface *s, Body *b) {
    if(!s || !b) return E_FAIL;
    if(!rnr_surfaceHasBody(s, b) && s->add(b) != S_OK) return E_FAIL;
    if(!rnr_bodyHasSurface(b, s) && b->add(s) != S_OK) return E_FAIL;
    return S_OK;
}

static HRESULT rnr_detachBody(Surface *s, Body *b) {
    if(!s || !b) return E_FAIL;
    if(rnr_surfaceHasBody(s, b) && s->remove(b) != S_OK) return E_FAIL;
    if(rnr_bodyHasSurface(b, s) && b->remove(s) != S_OK) return E_FAIL;
    return S_OK;
}

static void rnr_refreshAfterSurgery(const std::vector<Surface*> &surfaces, const std::vector<Body*> &extraBodies) {
    std::unordered_set<int> seenSurfaces;
    std::vector<Surface*> liveSurfaces;
    liveSurfaces.reserve(surfaces.size());
    for(auto *s : surfaces) {
        if(!s || s->objectId() < 0 || seenSurfaces.count(s->objectId())) continue;
        seenSurfaces.insert(s->objectId());
        liveSurfaces.push_back(s);
        s->positionChanged();
    }
    for(auto *s : liveSurfaces)
        s->refreshBodies();

    std::unordered_set<int> seenBodies;
    std::vector<Body*> liveBodies;
    for(auto *b : extraBodies) {
        if(!b || b->objectId() < 0 || seenBodies.count(b->objectId())) continue;
        seenBodies.insert(b->objectId());
        liveBodies.push_back(b);
    }
    for(auto *s : liveSurfaces)
        for(auto *b : s->getBodies()) {
            if(!b || b->objectId() < 0 || seenBodies.count(b->objectId())) continue;
            seenBodies.insert(b->objectId());
            liveBodies.push_back(b);
        }
    for(auto *b : liveBodies)
        b->positionChanged();

    std::unordered_set<int> seenVertices;
    std::vector<Vertex*> liveVertices;
    for(auto *s : liveSurfaces)
        for(auto *v : s->getVertices()) {
            if(!v || v->objectId() < 0 || seenVertices.count(v->objectId())) continue;
            seenVertices.insert(v->objectId());
            liveVertices.push_back(v);
        }
    for(auto *v : liveVertices)
        v->updateConnectedVertices();
}

} // anonymous namespace


/**
 * The native 3D T1 / Okuda I<->H reconnection operation. Phase B: the candidate (a short
 * interior edge for I->H, or a small triangular face for H->I) is identified by the scan and
 * stored by id; check() re-validates by re-walking the neighborhood and re-applying the
 * Condition-4 veto (handles may have been invalidated by an earlier op in the same pass), and
 * prep() re-fetches + gathers the affected bodies. implement() is still a no-op stub -- the
 * Okuda Appendix-1 vertex placement + manual surface surgery land in Phase C, mirroring the
 * structure of the ops above. Translated from our validated Python prototype rnr/reconnect.py,
 * not copied from the GPL tvm/3DVertVor reference.
 */
struct ReconnectionOperation : MeshQualityOperation {

    enum Kind { I2H = 0, H2I = 1 };

    Kind kind;
    int v10Id, v11Id;          // I2H: the short-edge endpoints (H2I: -1)
    int triId;                 // H2I: the triangle surface (I2H: -1)
    FloatP_t dlTh;             // reconnect length (Okuda Delta_l_th), for Phase-C placement
    FloatP_t hysteresis;       // placement hysteresis (Phase C)
    bool energyGate;           // optional greedy gate (Phase C+); OFF by default
    bool enforceTrigger;       // live doQuality ops re-check Condition 2; force ops bypass it
    std::unordered_set<int> affectedChildren;
    bool lastOk = false;
    std::string lastReason;
    int lastNewSurfaceId = -1;
    std::vector<int> lastNewVertexIds;

    /** I->H: short edge (cfg.v10, cfg.v11) -> triangular face. */
    ReconnectionOperation(Mesh *_mesh, const RNR_IConfig &cfg,
                          FloatP_t _dlTh, FloatP_t _hys, bool _eg,
                          bool _enforceTrigger=true) : MeshQualityOperation(_mesh) {
        flags = MeshQualityOperation::Flag::Active;
        kind = I2H;
        v10Id = cfg.v10Id; v11Id = cfg.v11Id; triId = -1;
        dlTh = _dlTh; hysteresis = _hys; energyGate = _eg;
        enforceTrigger = _enforceTrigger;
        // targets = every surface this op modifies (3 side + 3 top + 3 bottom faces), so the
        // dependency graph serializes any other op touching one of them.
        std::unordered_set<int> t;
        for(auto &a : cfg.arms) t.insert(a.sideSurfaceId);
        for(auto &kv : cfg.topFaces) t.insert(kv.second->objectId());
        for(auto &kv : cfg.bottomFaces) t.insert(kv.second->objectId());
        targets.assign(t.begin(), t.end());
    }

    /** H->I: triangular face -> short edge. */
    ReconnectionOperation(Mesh *_mesh, const RNR_HConfig &cfg,
                          FloatP_t _dlTh, FloatP_t _hys, bool _eg,
                          bool _enforceTrigger=true) : MeshQualityOperation(_mesh) {
        flags = MeshQualityOperation::Flag::Active;
        kind = H2I;
        v10Id = -1; v11Id = -1; triId = cfg.triangleId;
        dlTh = _dlTh; hysteresis = _hys; energyGate = _eg;
        enforceTrigger = _enforceTrigger;
        std::unordered_set<int> t;
        t.insert(cfg.triangleId);
        for(auto &a : cfg.arms) t.insert(a.sideSurfaceId);
        for(auto &kv : cfg.topFaces) t.insert(kv.second->objectId());
        for(auto &kv : cfg.bottomFaces) t.insert(kv.second->objectId());
        targets.assign(t.begin(), t.end());
    }

    bool check() override {
        if(kind == I2H) {
            Vertex *v10 = mesh->getVertex(v10Id);
            Vertex *v11 = mesh->getVertex(v11Id);
            if(!v10 || !v11) return false;
            RNR_IConfig cfg;
            if(!rnr_iNeighbourhood(v10, v11, cfg)) return false;
            if(enforceTrigger && cfg.length >= dlTh) return false;
            return rnr_iToHVeto(cfg).empty();
        } else {
            Surface *tri = mesh->getSurface(triId);
            if(!tri) return false;
            RNR_HConfig cfg;
            if(!rnr_hNeighbourhood(tri, cfg)) return false;
            if(enforceTrigger && cfg.maxEdge >= dlTh) return false;
            return rnr_hToIVeto(cfg).empty();
        }
    }

    void prep() override {
        affectedChildren.clear();
        if(kind == I2H) {
            Vertex *v10 = mesh->getVertex(v10Id);
            Vertex *v11 = mesh->getVertex(v11Id);
            if(v10) for(auto *b : v10->getBodies()) affectedChildren.insert(b->objectId());
            if(v11) for(auto *b : v11->getBodies()) affectedChildren.insert(b->objectId());
        } else {
            Surface *tri = mesh->getSurface(triId);
            RNR_HConfig cfg;
            if(tri && rnr_hNeighbourhood(tri, cfg)) {
                affectedChildren.insert(cfg.capTopId);
                affectedChildren.insert(cfg.capBotId);
                for(int id : cfg.sideCellIds) affectedChildren.insert(id);
            }
            else if(tri) for(auto *b : tri->getBodies()) affectedChildren.insert(b->objectId());
        }
    }

    size_t numNewVertices() const override { return kind == I2H ? 3 : 2; }
    size_t numNewSurfaces() const override { return kind == I2H ? 1 : 0; }

    FloatP_t placementLength() const { return dlTh * (1 + hysteresis); }

    std::vector<int> fail(const std::string &reason) {
        lastOk = false;
        lastReason = reason;
        lastNewSurfaceId = -1;
        lastNewVertexIds.clear();
        return {};
    }

    std::vector<int> implementIToH() {
        Vertex *v10 = mesh->getVertex(v10Id);
        Vertex *v11 = mesh->getVertex(v11Id);
        if(!v10 || !v11) return fail("missing vertex");

        RNR_IConfig cfg;
        if(!rnr_iNeighbourhood(v10, v11, cfg)) return fail("no valid I-configuration");
        std::string veto = rnr_iToHVeto(cfg);
        if(!veto.empty()) return fail(veto);

        const FloatP_t placeDl = placementLength();
        if(placeDl <= 0) return fail("reconnectLength must be > 0");

        std::array<FVector3, 3> positions;
        if(!rnr_placeIToH(cfg, placeDl, positions)) return fail("I->H placement failed");

        SurfaceType *stype = cfg.arms[0].sideSurface->type();
        if(!stype) return fail("trigger surface has no SurfaceType");

        std::unordered_map<int, size_t> armByOuterTop, armByOuterBot;
        for(size_t i = 0; i < cfg.arms.size(); i++) {
            armByOuterTop[cfg.arms[i].outerTop->objectId()] = i;
            armByOuterBot[cfg.arms[i].outerBot->objectId()] = i;
        }

        // Preflight the ordered face edits before mutating; these mirror reconnect.py.
        for(auto &kv : cfg.topFaces) {
            Vertex *nextV, *prevV;
            std::tie(nextV, prevV) = kv.second->neighborVertices(cfg.v10);
            if(!prevV || !nextV || !armByOuterTop.count(prevV->objectId()) || !armByOuterTop.count(nextV->objectId()))
                return fail("top face ring-neighbours of v10 are not arm outer_top verts");
        }
        for(auto &kv : cfg.bottomFaces) {
            Vertex *nextV, *prevV;
            std::tie(nextV, prevV) = kv.second->neighborVertices(cfg.v11);
            if(!prevV || !nextV || !armByOuterBot.count(prevV->objectId()) || !armByOuterBot.count(nextV->objectId()))
                return fail("bottom face ring-neighbours of v11 are not arm outer_bot verts");
        }

        MeshSolver::engineLock();

        std::array<Vertex*, 3> tri = {nullptr, nullptr, nullptr};
        std::vector<Vertex*> created;
        for(size_t i = 0; i < positions.size(); i++) {
            VertexHandle vh = Vertex::create(positions[i]);
            tri[i] = vh ? vh.vertex() : nullptr;
            if(!tri[i]) {
                for(auto *v : created) v->destroy();
                MeshSolver::engineUnlock();
                return fail("failed to create triangle vertex");
            }
            created.push_back(tri[i]);
        }

        SurfaceHandle tHandle = (*stype)({
            VertexHandle(tri[0]->objectId()),
            VertexHandle(tri[1]->objectId()),
            VertexHandle(tri[2]->objectId())
        });
        Surface *T = tHandle ? tHandle.surface() : nullptr;
        if(!T) {
            for(auto *v : created) v->destroy();
            MeshSolver::engineUnlock();
            return fail("failed to create triangle surface");
        }

        std::vector<Surface*> touchedSurfs;
        touchedSurfs.reserve(10);

        // (1) SIDE faces: [outer_top, v10, v11, outer_bot] -> [outer_top, vt_k, outer_bot].
        for(size_t i = 0; i < cfg.arms.size(); i++) {
            Surface *s = cfg.arms[i].sideSurface;
            rnr_replaceV(s, cfg.v10, tri[i]);
            rnr_dropV(s, cfg.v11);
            touchedSurfs.push_back(s);
        }

        // (2) TOP faces: v10 -> triangle edge (vt_prev, vt_next).
        for(auto &kv : cfg.topFaces) {
            Surface *face = kv.second;
            Vertex *nextV, *prevV;
            std::tie(nextV, prevV) = face->neighborVertices(cfg.v10);
            Vertex *vtPrev = tri[armByOuterTop[prevV->objectId()]];
            Vertex *vtNext = tri[armByOuterTop[nextV->objectId()]];
            rnr_replaceV(face, cfg.v10, vtPrev);
            rnr_insertBetween(face, vtNext, vtPrev, nextV);
            touchedSurfs.push_back(face);
        }

        // (3) BOTTOM faces: v11 -> triangle edge, mirror of top.
        for(auto &kv : cfg.bottomFaces) {
            Surface *face = kv.second;
            Vertex *nextV, *prevV;
            std::tie(nextV, prevV) = face->neighborVertices(cfg.v11);
            Vertex *vtPrev = tri[armByOuterBot[prevV->objectId()]];
            Vertex *vtNext = tri[armByOuterBot[nextV->objectId()]];
            rnr_replaceV(face, cfg.v11, vtPrev);
            rnr_insertBetween(face, vtNext, vtPrev, nextV);
            touchedSurfs.push_back(face);
        }

        // (4) New triangular cap-cap contact.
        rnr_attachBody(T, cfg.capTop);
        rnr_attachBody(T, cfg.capBot);
        touchedSurfs.push_back(T);

        lastNewSurfaceId = T->objectId();
        lastNewVertexIds = {tri[0]->objectId(), tri[1]->objectId(), tri[2]->objectId()};

        // (5) Destroy the now-orphaned short-edge vertices, after all rewiring is done.
        cfg.v10->destroy();
        cfg.v11->destroy();

        rnr_refreshAfterSurgery(touchedSurfs, {cfg.capTop, cfg.capBot});

        MeshSolver::engineUnlock();

        lastOk = true;
        lastReason = "";
        next.clear();
        return std::vector<int>(affectedChildren.begin(), affectedChildren.end());
    }

    std::vector<int> implementHToI() {
        Surface *triSurf = mesh->getSurface(triId);
        if(!triSurf) return fail("missing triangle surface");

        RNR_HConfig cfg;
        if(!rnr_hNeighbourhood(triSurf, cfg)) return fail("no valid H-configuration");
        std::string veto = rnr_hToIVeto(cfg);
        if(!veto.empty()) return fail(veto);

        const FloatP_t placeDl = placementLength();
        if(placeDl <= 0) return fail("reconnectLength must be > 0");

        FVector3 p10, p11;
        if(!rnr_placeHToI(cfg, placeDl, p10, p11)) return fail("H->I placement failed");

        const std::set<int> triIds(cfg.triVertexIds.begin(), cfg.triVertexIds.end());
        for(auto &kv : cfg.topFaces) {
            int n = 0;
            for(auto *v : kv.second->getVertices()) if(triIds.count(v->objectId())) n++;
            if(n != 2) return fail("top face expected 2 triangle vertices");
        }
        for(auto &kv : cfg.bottomFaces) {
            int n = 0;
            for(auto *v : kv.second->getVertices()) if(triIds.count(v->objectId())) n++;
            if(n != 2) return fail("bottom face expected 2 triangle vertices");
        }

        MeshSolver::engineLock();

        VertexHandle nv10Handle = Vertex::create(p10);
        VertexHandle nv11Handle = Vertex::create(p11);
        Vertex *nv10 = nv10Handle ? nv10Handle.vertex() : nullptr;
        Vertex *nv11 = nv11Handle ? nv11Handle.vertex() : nullptr;
        if(!nv10 || !nv11) {
            if(nv10) nv10->destroy();
            if(nv11) nv11->destroy();
            MeshSolver::engineUnlock();
            return fail("failed to create recovered edge vertices");
        }

        std::vector<Surface*> touchedSurfs;
        touchedSurfs.reserve(9);

        // (1) SIDE faces: [outer_top, vt_k, outer_bot] -> [outer_top, nv10, nv11, outer_bot].
        for(auto &a : cfg.arms) {
            Surface *s = a.sideSurface;
            rnr_replaceV(s, a.triVertex, nv10);
            rnr_insertBetween(s, nv11, nv10, a.outerBot);
            touchedSurfs.push_back(s);
        }

        // (2) TOP faces: triangle edge -> nv10.
        for(auto &kv : cfg.topFaces) {
            Surface *face = kv.second;
            std::vector<Vertex*> present;
            for(auto *v : face->getVertices()) if(triIds.count(v->objectId())) present.push_back(v);
            rnr_replaceV(face, present[0], nv10);
            rnr_dropV(face, present[1]);
            touchedSurfs.push_back(face);
        }

        // (3) BOTTOM faces: triangle edge -> nv11.
        for(auto &kv : cfg.bottomFaces) {
            Surface *face = kv.second;
            std::vector<Vertex*> present;
            for(auto *v : face->getVertices()) if(triIds.count(v->objectId())) present.push_back(v);
            rnr_replaceV(face, present[0], nv11);
            rnr_dropV(face, present[1]);
            touchedSurfs.push_back(face);
        }

        std::vector<int> oldTriVertexIds = cfg.triVertexIds;
        Surface *T = cfg.triangle;
        rnr_detachBody(T, cfg.capTop);
        rnr_detachBody(T, cfg.capBot);
        T->destroy();

        // SurfaceHandle::destroy/member Surface::destroy does not cascade-delete orphan
        // vertices; explicitly remove the three now-orphaned triangle vertices.
        for(int id : oldTriVertexIds) {
            Vertex *v = mesh->getVertex(id);
            if(v && v->getSurfaces().empty())
                v->destroy();
        }

        lastNewSurfaceId = -1;
        lastNewVertexIds = {nv10->objectId(), nv11->objectId()};

        rnr_refreshAfterSurgery(touchedSurfs, {cfg.capTop, cfg.capBot});

        MeshSolver::engineUnlock();

        lastOk = true;
        lastReason = "";
        next.clear();
        return std::vector<int>(affectedChildren.begin(), affectedChildren.end());
    }

    std::vector<int> implement() override {
        lastOk = false;
        lastReason = "";
        lastNewSurfaceId = -1;
        lastNewVertexIds.clear();
        return kind == I2H ? implementIToH() : implementHToI();
    }
};


/////////////////
// MeshQuality //
/////////////////


HRESULT MeshQuality_constructChains(std::vector<MeshQualityOperation*> &ops) {
    auto func = [&ops](int i) -> void {
        MeshQualityOperation *op = ops[i];
        if(!op) return;
        MeshQualityOperation_checkChain(op, ops);
    };
    parallel_for(ops.size(), func);
    return S_OK;
}

static HRESULT MeshQuality_constructOperationsVertex(
    Mesh *mesh, 
    const std::vector<bool> &passMask, 
    const FloatP_t &edgeSplitDist, 
    const FloatP_t &vertexMergeDist, 
    std::vector<MeshQualityOperation*> &ops_active, 
    std::vector<MeshQualityOperation*> &op_heads
) {
    std::vector<MeshQualityOperation*> ops(mesh->sizeVertices(), 0);
    const FloatP_t vertexMergeDist2 = vertexMergeDist * vertexMergeDist;

    auto check_verts = [&mesh, &passMask, &ops, edgeSplitDist, vertexMergeDist2](int i) -> void {
        if(passMask[i]) return;

        Vertex *v = mesh->getVertex(i);
        if(!v) return;

        MeshQualityOperation *op = NULL;

        // Check for vertex split if the vertex defines four or more surfaces

        std::vector<Vertex*> v_nbs = v->connectedVertices();
        if(v_nbs.size() > 3) {
            FVector3 sep;
            std::vector<Vertex*> vert_nbs, new_vert_nbs;
            if(MeshQuality_vertexSplitTest(v, mesh, edgeSplitDist, sep, v_nbs, vert_nbs, new_vert_nbs)) {
                ops[i] = new VertexSplitOperation(mesh, v, sep, vert_nbs, new_vert_nbs);
                return;
            }

        }

        // Check vertex merge distance for neighbor vertices with a greater id

        FVector3 vpos = v->getPosition();
        for(auto &nv : v_nbs) {

            if(v->objectId() < nv->objectId()) {

                FVector3 nvpos = nv->getPosition();
                FVector3 nvrelPos = meshRelativePosition(vpos, nvpos);
                FloatP_t nvdist2 = nvrelPos.dot();
                if(nvdist2 < vertexMergeDist2) {
                    TF_Log(LOG_TRACE) << v->objectId() << ", " << nv->objectId() << ", " << nvrelPos;
                    
                    ops[i] = new EdgeDemoteOperation(mesh, v, nv);
                    return;
                }

            }
        }

    };
    parallel_for(mesh->sizeVertices(), check_verts);

    if(MeshQuality_constructChains(ops) != S_OK) { 
        for(size_t i = 0; i < ops.size(); i++) delete ops[i];
        ops.clear();
        return E_FAIL;
    }

    ops_active.clear();
    ops_active.reserve(ops.size());
    for(auto &op : ops) 
        if(op) 
            ops_active.push_back(op);
    std::set<MeshQualityOperation*> op_heads_set = MeshQualityOperation_headOperations(ops_active);
    op_heads = std::vector<MeshQualityOperation*>(op_heads_set.begin(), op_heads_set.end());

    return S_OK;
}

static HRESULT MeshQuality_constructOperationsSurface(
    Mesh *mesh, 
    const std::vector<bool> &passMask, 
    const FloatP_t &surfaceDemoteArea, 
    const bool &collision2D, 
    std::vector<MeshQualityOperation*> &ops_active, 
    std::vector<MeshQualityOperation*> &op_heads
) {
    std::vector<MeshQualityOperation*> ops(mesh->sizeSurfaces(), 0);

    auto check_surfs = [&mesh, &passMask, &ops, surfaceDemoteArea, collision2D](int i) -> void {
        if(passMask[i]) return;

        Surface *s = mesh->getSurface(i);
        if(!s) return;

        // Check for surface demote

        FloatP_t sArea = s->getArea();

        if(sArea < surfaceDemoteArea) {
            TF_Log(LOG_TRACE) << sArea;

            ops[i] = new SurfaceDemoteOperation(mesh, s);
            return;
        }

        if(!collision2D) 
            return;

        // Check for edge penetration

        //  Determine neighborhood search distance
        FloatP_t nbsSearchDist2 = 0;
        FVector3 centroid = s->getCentroid();
        for(auto &v : s->getVertices()) {
            FloatP_t thisVertDist2 = meshRelativePosition(centroid, v->getPosition()).dot();
            if(thisVertDist2 > nbsSearchDist2) 
                nbsSearchDist2 = thisVertDist2;
        }

        //  Get neighbors
        ParticleList nbs = metrics::neighborhoodParticles(centroid, FPTYPE_SQRT(nbsSearchDist2));

        //  Test each neighbor
        for(size_t j = 0; j < nbs.nr_parts; j++) {
            ParticleHandle *nb = nbs.item(j);
            Vertex *v_nb = mesh->getVertexByPID(nb->id);
            if(!v_nb) 
                continue;

            //  No self-intersecting
            if(v_nb->defines(s)) 
                continue;
            
            Vertex *va, *vb;
            if(s->contains(nb->getPosition(), &va, &vb)) {
                ops[i] = new VertexInsertOperation(mesh, v_nb, s, va, vb);
                return;
            }
        }

    };
    parallel_for(mesh->sizeSurfaces(), check_surfs);
    
    if(MeshQuality_constructChains(ops) != S_OK) { 
        for(size_t i = 0; i < ops.size(); i++) delete ops[i];
        ops.clear();
        return E_FAIL;
    }

    ops_active.clear();
    ops_active.reserve(ops.size());
    for(auto &op : ops) 
        if(op) 
            ops_active.push_back(op);
    std::set<MeshQualityOperation*> op_heads_set = MeshQualityOperation_headOperations(ops_active);
    op_heads = std::vector<MeshQualityOperation*>(op_heads_set.begin(), op_heads_set.end());

    return S_OK;
}

static HRESULT MeshQuality_constructOperationsBody(
    Mesh *mesh, 
    const std::vector<bool> &passMask, 
    const FloatP_t &bodyDemoteVolume, 
    std::vector<MeshQualityOperation*> &ops_active, 
    std::vector<MeshQualityOperation*> &op_heads
) {
    std::vector<MeshQualityOperation*> ops(mesh->sizeBodies(), 0);

    auto check_bodys = [&mesh, &passMask, &ops, bodyDemoteVolume](int i) -> void {
        if(passMask[i]) return;

        Body *b = mesh->getBody(i);
        if(!b) return;

        FloatP_t bvol = b->getVolume();
        
        if(bvol < bodyDemoteVolume) {
            TF_Log(LOG_TRACE) << bvol;

            ops[i] = new BodyDemoteOperation(mesh, b);
            return;
        }
    };
    parallel_for(mesh->sizeBodies(), check_bodys);
    
    if(MeshQuality_constructChains(ops) != S_OK) { 
        for(size_t i = 0; i < ops.size(); i++) delete ops[i];
        ops.clear();
        return E_FAIL;
    }

    ops_active.clear();
    ops_active.reserve(ops.size());
    for(auto &op : ops)
        if(op)
            ops_active.push_back(op);
    std::set<MeshQualityOperation*> op_heads_set = MeshQualityOperation_headOperations(ops_active);
    op_heads = std::vector<MeshQualityOperation*>(op_heads_set.begin(), op_heads_set.end());

    return S_OK;
}

static HRESULT MeshQuality_constructOperationsReconnection(
    Mesh *mesh,
    const std::vector<bool> &passMask,
    const FloatP_t &reconnectLength,
    const FloatP_t &reconnectHysteresis,
    const bool &reconnectEnergyGate,
    std::vector<MeshQualityOperation*> &ops_active,
    std::vector<MeshQualityOperation*> &op_heads
) {
    ops_active.clear();
    op_heads.clear();

    // Disabled unless a positive trigger length (Okuda Condition 2, Delta_l_th) is configured.
    if(reconnectLength <= 0)
        return S_OK;

    // Like the other scans, ops is indexed by mesh-object id (here surface id, since both
    // triggers live on a Surface: an I->H short edge is a consecutive vertex pair on a surface,
    // an H->I feature IS a triangular surface). targets reference this same id space, so
    // MeshQuality_constructChains serializes ops with overlapping touched surfaces.
    //
    // The scan is serial (the [I]/[H] walks read shared adjacency and build STL containers); a
    // parallel_for like the sibling scans is a possible later optimization. As of Phase C/D the
    // ops returned here mutate the mesh: each runs the real ReconnectionOperation::implement()
    // surgery when MeshQuality_doOperations walks the chain inside doQuality.
    std::vector<MeshQualityOperation*> ops(mesh->sizeSurfaces(), 0);

    // I->H candidates: short interior edges in a legal [I] neighborhood.
    std::vector<RNR_IConfig> iCands;
    rnr_findShortEdges(mesh, reconnectLength, iCands);
    for(auto &cfg : iCands) {
        if(!rnr_iToHVeto(cfg).empty())                 // Condition-4 veto -> skip
            continue;
        int slot = rnr_iTriggerSlot(cfg);
        if(slot < 0 || (size_t)slot >= ops.size() || passMask[slot] || ops[slot])
            continue;                                   // out of range / excluded / slot taken
        ops[slot] = new ReconnectionOperation(mesh, cfg, reconnectLength, reconnectHysteresis, reconnectEnergyGate);
    }

    // H->I candidates: small triangular faces in a legal [H] neighborhood.
    std::vector<RNR_HConfig> hCands;
    rnr_findSmallTriangles(mesh, reconnectLength, hCands);
    for(auto &cfg : hCands) {
        if(!rnr_hToIVeto(cfg).empty())
            continue;
        int slot = cfg.triangleId;
        if(slot < 0 || (size_t)slot >= ops.size() || passMask[slot] || ops[slot])
            continue;
        ops[slot] = new ReconnectionOperation(mesh, cfg, reconnectLength, reconnectHysteresis, reconnectEnergyGate);
    }

    // Build the dependency chains SERIALLY, not via the parallel MeshQuality_constructChains.
    // Reconnection ops have large (9-surface) target sets with heavy overlap between adjacent
    // edges, so the candidate graph is far denser than the other passes ever produce. The
    // parallel chain builder mutates the shared prev/next graph (and appendNext's loop-check
    // reads it) without fully locking `this->next`/upstreams(), which races on a dense graph and
    // can fabricate a cycle -> unbounded upstreams() recursion -> stack-overflow segfault
    // (observed as a nondeterministic crash on a repeated doQuality). A serial build cannot race
    // and yields a valid DAG, which the parallel MeshQuality_doOperations then walks safely.
    for(size_t i = 0; i < ops.size(); i++)
        if(ops[i])
            MeshQualityOperation_checkChain(ops[i], ops);

    ops_active.clear();
    ops_active.reserve(ops.size());
    for(auto &op : ops)
        if(op)
            ops_active.push_back(op);
    std::set<MeshQualityOperation*> op_heads_set = MeshQualityOperation_headOperations(ops_active);
    op_heads = std::vector<MeshQualityOperation*>(op_heads_set.begin(), op_heads_set.end());

    return S_OK;
}

static HRESULT MeshQuality_doOperations(MeshQualityOperation *op, std::unordered_set<int> &affectedChildren) {
    if(!op->check()) { 
        for(auto &n : op->next) {
            auto itr = std::find(n->prev.begin(), n->prev.end(), op);
            if(itr != n->prev.end()) {
                n->lock.lock();
                n->prev.erase(itr);
                n->lock.unlock();
            }
        }
        return S_OK;
    }

    for(auto &i : op->implement()) 
        affectedChildren.insert(i);

    for(auto &n : op->next) {
        auto itr = std::find(n->prev.begin(), n->prev.end(), op);
        HRESULT res = S_OK;

        n->lock.lock();
        n->prev.erase(itr);
        if(n->prev.empty()) 
            res = MeshQuality_doOperations(n, affectedChildren);
        n->lock.unlock();

        if(res != S_OK) 
            return res;
    }

    return S_OK;
}

static HRESULT MeshQuality_doOperations(
    Mesh *mesh, 
    std::vector<MeshQualityOperation*> &op_active, 
    std::vector<MeshQualityOperation*> &op_heads, 
    std::vector<int> &affectedChildren) 
{
    std::atomic<size_t> atomic_numNewVertices = 0;
    std::atomic<size_t> atomic_numNewSurfaces = 0;
    std::atomic<size_t> atomic_numNewBodies = 0;
    
    auto func_count = [&op_active, &atomic_numNewVertices, &atomic_numNewSurfaces, &atomic_numNewBodies](int tid) -> void {
        size_t _numNewVertices = 0;
        size_t _numNewSurfaces = 0;
        size_t _numNewBodies = 0;

        for(int i = tid; i < op_active.size();) {
            auto op = op_active[i];
            _numNewVertices += op->numNewVertices();
            _numNewSurfaces += op->numNewSurfaces();
            _numNewBodies += op->numNewBodies();
            i += ThreadPool::size();
        }

        atomic_numNewVertices.fetch_add(_numNewVertices);
        atomic_numNewSurfaces.fetch_add(_numNewSurfaces);
        atomic_numNewBodies.fetch_add(_numNewBodies);
    };
    parallel_for(ThreadPool::size(), func_count);

    mesh->ensureAvailableVertices(atomic_numNewVertices);
    mesh->ensureAvailableSurfaces(atomic_numNewSurfaces);
    mesh->ensureAvailableBodies(atomic_numNewBodies);
    
    parallel_for(op_active.size(), [&op_active](int i) -> void { op_active[i]->prep(); });

    static std::mutex affectedChildrenLock;
    std::unordered_set<int> affectedChildrenSet;
    auto func_do = [&op_heads, &affectedChildrenSet](int tid) -> void {
        std::unordered_set<int> affectedChildren_tid;
        for(int i = tid; i < op_heads.size();) {
            MeshQuality_doOperations(op_heads[i], affectedChildren_tid);
            i += ThreadPool::size();
        }
        affectedChildrenLock.lock();
        for(auto &i : affectedChildren_tid) 
            affectedChildrenSet.insert(i);
        affectedChildrenLock.unlock();
    };
    parallel_for(ThreadPool::size(), func_do);
    affectedChildren = std::vector<int>(affectedChildrenSet.begin(), affectedChildrenSet.end());

    return S_OK;
}

// Serial executor for the reconnection pass.
//
// MeshQuality_doOperations (above) runs the head operations' implement() walks IN PARALLEL
// (parallel_for over op_heads). The dependency graph only serializes ops whose `targets` overlap,
// and a ReconnectionOperation's `targets` are just its 9 incident SURFACES -- NOT the 6 outer
// vertices or 5 bodies its I<->H surgery also mutates. So two reconnections with disjoint surface
// targets but a shared outer vertex/body can implement() concurrently and race on that shared
// (non-target) object, corrupting its surface/body lists -> intermittent heap corruption / segfault
// (seen ~once per dozen full gate runs; volumes were perfectly stable right up to the crash, so it
// is a scheduler race, not the dynamics). Widening `targets` to vertices/bodies would change the
// validated walk; instead we run THIS pass serially. It is throttled (reconnectInterval) and sparse
// (a handful of ops per call), so serial cost is negligible. This extends the same serialization
// rationale MeshQuality_constructOperationsReconnection already uses for chain construction.
static HRESULT MeshQuality_doOperationsReconnectionSerial(
    Mesh *mesh,
    std::vector<MeshQualityOperation*> &op_active,
    std::vector<MeshQualityOperation*> &op_heads,
    std::vector<int> &affectedChildren)
{
    size_t numNewVertices = 0, numNewSurfaces = 0, numNewBodies = 0;
    for(auto op : op_active) {
        numNewVertices += op->numNewVertices();
        numNewSurfaces += op->numNewSurfaces();
        numNewBodies += op->numNewBodies();
    }
    mesh->ensureAvailableVertices(numNewVertices);
    mesh->ensureAvailableSurfaces(numNewSurfaces);
    mesh->ensureAvailableBodies(numNewBodies);

    for(auto op : op_active)
        op->prep();

    std::unordered_set<int> affectedChildrenSet;
    for(auto head : op_heads)
        MeshQuality_doOperations(head, affectedChildrenSet);   // single-op recursive walk, serial
    affectedChildren = std::vector<int>(affectedChildrenSet.begin(), affectedChildrenSet.end());

    return S_OK;
}

static HRESULT MeshQuality_clearOperations(std::vector<MeshQualityOperation*> &ops) {
    auto func = [&ops](int i) -> void {
        delete ops[i];
        ops[i] = NULL;
    };
    parallel_for(ops.size(), func);

    return S_OK;
}

MeshQuality::MeshQuality(
    const FloatP_t &vertexMergeDistCf,
    const FloatP_t &surfaceDemoteAreaCf,
    const FloatP_t &bodyDemoteVolumeCf,
    const FloatP_t &_edgeSplitDistCf,
    const FloatP_t &_reconnectLength
) :
    reconnectLength{_reconnectLength},
    reconnectHysteresis{0.0},
    reconnectEnergyGate{false},
    stockQualityOps{true},
    reconnectInterval{1},
    collision2D{true},
    _working{false},
    reconnectCounter{0}
{
    FloatP_t uvolu = Universe::dim().product();
    FloatP_t uleng = std::cbrt(uvolu);
    FloatP_t uarea = uleng * uleng;

    vertexMergeDist = uleng * vertexMergeDistCf;
    surfaceDemoteArea = uarea * surfaceDemoteAreaCf;
    bodyDemoteVolume = uvolu * bodyDemoteVolumeCf;
    edgeSplitDist = _edgeSplitDistCf * vertexMergeDist;
    // reconnectLength is an ABSOLUTE length (Okuda Delta_l_th), not box-scaled like the others;
    // default 0 keeps the reconnection pass disabled (a no-op) until configured.
}

std::string MeshQuality::str() const {
    std::stringstream ss;

    ss << "MeshQuality(";
    ss << "vertexMergeDist="    << this->vertexMergeDist                << ", ";
    ss << "surfaceDemoteArea="  << this->surfaceDemoteArea              << ", ";
    ss << "bodyDemoteVolume="   << this->bodyDemoteVolume               << ", ";
    ss << "edgeSplitDist="      << this->edgeSplitDist                  << ", ";
    ss << "reconnectLength="    << this->reconnectLength                << ", ";
    ss << "reconnectHysteresis="<< this->reconnectHysteresis            << ", ";
    ss << "reconnectEnergyGate="<< (this->reconnectEnergyGate ? "yes" : "no") << ", ";
    ss << "stockQualityOps="    << (this->stockQualityOps ? "yes" : "no") << ", ";
    ss << "reconnectInterval="  << this->reconnectInterval              << ", ";
    ss << "collision2D="        << (this->collision2D ? "yes" : "no")   << ", ";
    ss << "working="            << (this->_working    ? "yes" : "no");
    ss << ")";

    return ss.str();
}

HRESULT MeshQuality::doQuality() { 

    _working = true;

    Mesh *mesh = Mesh::get();

    std::vector<MeshQualityOperation*> op_active, op_heads;
    std::vector<int> affectedChildren;
    std::vector<bool> passMask;

    // Stock TissueForge quality checks. These legacy passes are preserved by default, but
    // Phase-D native RNR harnesses disable them so doQuality() exercises only the Okuda
    // reconnection pass; the stock degenerate-collapse passes are known to crash on finite
    // Kelvin blocks independently of RNR (see rnr/PORTING_NOTES.md section 6b).
    if(stockQualityOps) {

    // Vertex checks
    
    passMask = std::vector<bool>(mesh->sizeVertices(), false);
    for(auto &i : excludedVertices) 
        if(i < passMask.size()) 
            passMask[i] = true;
    if(MeshQuality_constructOperationsVertex(mesh, passMask, edgeSplitDist, vertexMergeDist, op_active, op_heads) != S_OK || 
        MeshQuality_doOperations(mesh, op_active, op_heads, affectedChildren) != S_OK || 
        MeshQuality_clearOperations(op_active) != S_OK) {
        _working = false;
        return E_FAIL;
    }

    // Surface checks

    passMask = std::vector<bool>(mesh->sizeSurfaces(), false);
    std::unordered_set<int> affectedBodiesImpl;
    for(auto &i : affectedChildren) {
        passMask[i] = true;
        Surface *s = mesh->getSurface(i);
        if(s) 
            for(auto &b : s->getBodies()) 
                affectedBodiesImpl.insert(b->objectId());
    }
    for(auto &i : excludedSurfaces) 
        if(i < passMask.size()) 
            passMask[i] = true;
    if(MeshQuality_constructOperationsSurface(mesh, passMask, surfaceDemoteArea, collision2D, op_active, op_heads) != S_OK || 
        MeshQuality_doOperations(mesh, op_active, op_heads, affectedChildren) != S_OK || 
        MeshQuality_clearOperations(op_active) != S_OK) {
        _working = false;
        return E_FAIL;
    }

    // Body checks

    passMask = std::vector<bool>(mesh->sizeBodies(), false);
    for(auto &i : affectedChildren) 
        passMask[i] = true;
    for(auto &i : affectedBodiesImpl) 
        passMask[i] = true;
    for(auto &i : excludedBodies) 
        if(i < passMask.size()) 
            passMask[i] = true;
    if(MeshQuality_constructOperationsBody(mesh, passMask, bodyDemoteVolume, op_active, op_heads) != S_OK ||
        MeshQuality_doOperations(mesh, op_active, op_heads, affectedChildren) != S_OK ||
        MeshQuality_clearOperations(op_active) != S_OK) {
        _working = false;
        return E_FAIL;
    }

    }

    // Reconnection checks (native 3D T1 / Okuda I<->H RNR). The trigger lives on surfaces
    // (a short interior edge is a consecutive vertex pair on a surface; an H-state feature is a
    // triangular surface), so the pass scans surfaces. It is a no-op unless reconnectLength > 0.
    //
    // Throttle (the 3DVertVor oracle's dtr): the reconnection pass runs only on every
    // reconnectInterval-th doQuality() call. Between passes the integrator relaxes the mesh, which
    // breaks the post-reconnection overshoot storm that reconnecting every step produces (oracle
    // reconnects every dtr = 10*dt; see rnr/PORTING_NOTES.md 6c/6d). reconnectInterval = 1 (the
    // default) preserves the original every-step behavior. The counter advances every doQuality()
    // call regardless so the cadence is measured in integration steps.
    const unsigned int interval = reconnectInterval < 1 ? 1 : reconnectInterval;
    const bool reconnectDue = (this->reconnectCounter % interval) == 0;
    this->reconnectCounter++;
    if(reconnectDue) {
        passMask = std::vector<bool>(mesh->sizeSurfaces(), false);
        for(auto &i : excludedSurfaces)
            if(i < passMask.size())
                passMask[i] = true;
        if(MeshQuality_constructOperationsReconnection(mesh, passMask, reconnectLength, reconnectHysteresis, reconnectEnergyGate, op_active, op_heads) != S_OK ||
            MeshQuality_doOperationsReconnectionSerial(mesh, op_active, op_heads, affectedChildren) != S_OK ||
            MeshQuality_clearOperations(op_active) != S_OK) {
            _working = false;
            return E_FAIL;
        }
    }

    _working = false;

    return S_OK;
}

HRESULT MeshQuality::setVertexMergeDistance(const FloatP_t &_val) {
    if(_val < 0) 
        return E_FAIL;
    vertexMergeDist = _val;
    return S_OK;
}

HRESULT MeshQuality::setSurfaceDemoteArea(const FloatP_t &_val) {
    if(_val < 0) 
        return E_FAIL;
    surfaceDemoteArea = _val;
    return S_OK;
}

HRESULT MeshQuality::setBodyDemoteVolume(const FloatP_t &_val) {
    if(_val < 0) 
        return E_FAIL;
    bodyDemoteVolume = _val;
    return S_OK;
}

HRESULT MeshQuality::setEdgeSplitDist(const FloatP_t &_val) {
    if(_val <= 0)
        return E_FAIL;
    edgeSplitDist = _val;
    return S_OK;
}

HRESULT MeshQuality::setReconnectLength(const FloatP_t &_val) {
    if(_val < 0)
        return E_FAIL;
    reconnectLength = _val;
    return S_OK;
}

HRESULT MeshQuality::setReconnectHysteresis(const FloatP_t &_val) {
    if(_val < 0)
        return E_FAIL;
    reconnectHysteresis = _val;
    return S_OK;
}

HRESULT MeshQuality::setReconnectEnergyGate(const bool &_val) {
    reconnectEnergyGate = _val;
    return S_OK;
}

HRESULT MeshQuality::setStockQualityOps(const bool &_val) {
    stockQualityOps = _val;
    return S_OK;
}

HRESULT MeshQuality::setReconnectInterval(const unsigned int &_val) {
    reconnectInterval = _val < 1 ? 1 : _val;
    return S_OK;
}

std::string MeshQuality::analyzeIReconnection(const unsigned int &v10Id, const unsigned int &v11Id) const {
    Mesh *mesh = Mesh::get();
    if(!mesh) return "{\"valid\": false, \"reason\": \"no mesh\"}";
    Vertex *v10 = mesh->getVertex(v10Id);
    Vertex *v11 = mesh->getVertex(v11Id);
    if(!v10 || !v11) return "{\"valid\": false, \"reason\": \"missing vertex\"}";
    RNR_IConfig cfg;
    if(!rnr_iNeighbourhood(v10, v11, cfg)) return "{\"valid\": false}";
    return rnr_iConfigJson(cfg, rnr_iToHVeto(cfg));
}

std::string MeshQuality::analyzeHReconnection(const unsigned int &triId) const {
    Mesh *mesh = Mesh::get();
    if(!mesh) return "{\"valid\": false, \"reason\": \"no mesh\"}";
    Surface *tri = mesh->getSurface(triId);
    if(!tri) return "{\"valid\": false, \"reason\": \"missing surface\"}";
    RNR_HConfig cfg;
    if(!rnr_hNeighbourhood(tri, cfg)) return "{\"valid\": false}";
    return rnr_hConfigJson(cfg, rnr_hToIVeto(cfg));
}

std::string MeshQuality::findReconnectionCandidates() const {
    std::stringstream ss;
    ss << "[";
    Mesh *mesh = Mesh::get();
    if(mesh && reconnectLength > 0) {
        bool first = true;
        std::vector<RNR_IConfig> iCands;
        rnr_findShortEdges(mesh, reconnectLength, iCands);
        for(auto &cfg : iCands) {
            if(!first) ss << ", ";
            ss << rnr_iConfigJson(cfg, rnr_iToHVeto(cfg));
            first = false;
        }
        std::vector<RNR_HConfig> hCands;
        rnr_findSmallTriangles(mesh, reconnectLength, hCands);
        for(auto &cfg : hCands) {
            if(!first) ss << ", ";
            ss << rnr_hConfigJson(cfg, rnr_hToIVeto(cfg));
            first = false;
        }
    }
    ss << "]";
    return ss.str();
}

std::string MeshQuality::forceReconnectIToH(const unsigned int &v10Id, const unsigned int &v11Id) const {
    Mesh *mesh = Mesh::get();
    if(!mesh) return rnr_reconnectResultJson(false, "no mesh", -1, {});
    if(reconnectLength <= 0) return rnr_reconnectResultJson(false, "reconnectLength must be > 0", -1, {});

    Vertex *v10 = mesh->getVertex(v10Id);
    Vertex *v11 = mesh->getVertex(v11Id);
    if(!v10 || !v11) return rnr_reconnectResultJson(false, "missing vertex", -1, {});

    RNR_IConfig cfg;
    if(!rnr_iNeighbourhood(v10, v11, cfg))
        return rnr_reconnectResultJson(false, "no valid I-configuration", -1, {});
    std::string veto = rnr_iToHVeto(cfg);
    if(!veto.empty())
        return rnr_reconnectResultJson(false, veto, -1, {});

    ReconnectionOperation op(mesh, cfg, reconnectLength, reconnectHysteresis, reconnectEnergyGate, false);
    op.prep();
    if(!op.check())
        return rnr_reconnectResultJson(false, "operation no longer valid", -1, {});

    if(mesh->ensureAvailableVertices(op.numNewVertices()) != S_OK ||
       mesh->ensureAvailableSurfaces(op.numNewSurfaces()) != S_OK ||
       mesh->ensureAvailableBodies(op.numNewBodies()) != S_OK)
        return rnr_reconnectResultJson(false, "failed to reserve mesh storage", -1, {});

    op.implement();
    return rnr_reconnectResultJson(op.lastOk, op.lastReason, op.lastNewSurfaceId, op.lastNewVertexIds);
}

std::string MeshQuality::forceReconnectHToI(const unsigned int &triId) const {
    Mesh *mesh = Mesh::get();
    if(!mesh) return rnr_reconnectResultJson(false, "no mesh", -1, {});
    if(reconnectLength <= 0) return rnr_reconnectResultJson(false, "reconnectLength must be > 0", -1, {});

    Surface *tri = mesh->getSurface(triId);
    if(!tri) return rnr_reconnectResultJson(false, "missing triangle surface", -1, {});

    RNR_HConfig cfg;
    if(!rnr_hNeighbourhood(tri, cfg))
        return rnr_reconnectResultJson(false, "no valid H-configuration", -1, {});
    std::string veto = rnr_hToIVeto(cfg);
    if(!veto.empty())
        return rnr_reconnectResultJson(false, veto, -1, {});

    ReconnectionOperation op(mesh, cfg, reconnectLength, reconnectHysteresis, reconnectEnergyGate, false);
    op.prep();
    if(!op.check())
        return rnr_reconnectResultJson(false, "operation no longer valid", -1, {});

    if(mesh->ensureAvailableVertices(op.numNewVertices()) != S_OK ||
       mesh->ensureAvailableSurfaces(op.numNewSurfaces()) != S_OK ||
       mesh->ensureAvailableBodies(op.numNewBodies()) != S_OK)
        return rnr_reconnectResultJson(false, "failed to reserve mesh storage", -1, {});

    op.implement();
    return rnr_reconnectResultJson(op.lastOk, op.lastReason, op.lastNewSurfaceId, op.lastNewVertexIds);
}

HRESULT MeshQuality::setCollision2D(const bool &_collision2D) {
    collision2D = _collision2D;
    return S_OK;
}

HRESULT MeshQuality::excludeVertex(const unsigned int &id) {
    excludedVertices.insert(id);
    return S_OK;
}

HRESULT MeshQuality::excludeSurface(const unsigned int &id) {
    excludedSurfaces.insert(id);
    return S_OK;
}

HRESULT MeshQuality::excludeBody(const unsigned int &id) {
    excludedBodies.insert(id);
    return S_OK;
}

HRESULT MeshQuality::includeVertex(const unsigned int &id) {
    excludedVertices.erase(id);
    return S_OK;
}

HRESULT MeshQuality::includeSurface(const unsigned int &id) {
    excludedSurfaces.erase(id);
    return S_OK;
}

HRESULT MeshQuality::includeBody(const unsigned int &id) {
    excludedBodies.erase(id);
    return S_OK;
}

namespace TissueForge::io {


    template <>
    HRESULT toFile(const TissueForge::models::vertex::MeshQuality &dataElement, const MetaData &metaData, IOElement &fileElement) {

        TF_IOTOEASY(fileElement, metaData, "vertexMergeDist", dataElement.getVertexMergeDistance());
        TF_IOTOEASY(fileElement, metaData, "surfaceDemoteArea", dataElement.getSurfaceDemoteArea());
        TF_IOTOEASY(fileElement, metaData, "bodyDemoteVolume", dataElement.getBodyDemoteVolume());
        TF_IOTOEASY(fileElement, metaData, "edgeSplitDist", dataElement.getEdgeSplitDist());
        TF_IOTOEASY(fileElement, metaData, "reconnectLength", dataElement.getReconnectLength());
        TF_IOTOEASY(fileElement, metaData, "reconnectHysteresis", dataElement.getReconnectHysteresis());
        TF_IOTOEASY(fileElement, metaData, "reconnectEnergyGate", dataElement.getReconnectEnergyGate());
        TF_IOTOEASY(fileElement, metaData, "stockQualityOps", dataElement.getStockQualityOps());
        TF_IOTOEASY(fileElement, metaData, "reconnectInterval", dataElement.getReconnectInterval());
        TF_IOTOEASY(fileElement, metaData, "collision2D", dataElement.getCollision2D());
        TF_IOTOEASY(fileElement, metaData, "excludedVertices", dataElement.getExcludedVertices());
        TF_IOTOEASY(fileElement, metaData, "excludedSurfaces", dataElement.getExcludedSurfaces());
        TF_IOTOEASY(fileElement, metaData, "excludedBodies", dataElement.getExcludedBodies());

        fileElement.get()->type = "MeshQuality";

        return S_OK;
    }

    template <>
    HRESULT fromFile(const IOElement &fileElement, const MetaData &metaData, TissueForge::models::vertex::MeshQuality *dataElement) {
        
        FloatP_t vertexMergeDist;
        TF_IOFROMEASY(fileElement, metaData, "vertexMergeDist", &vertexMergeDist);
        dataElement->setVertexMergeDistance(vertexMergeDist);

        FloatP_t surfaceDemoteArea;
        TF_IOFROMEASY(fileElement, metaData, "surfaceDemoteArea", &surfaceDemoteArea);
        dataElement->setSurfaceDemoteArea(surfaceDemoteArea);

        FloatP_t bodyDemoteVolume;
        TF_IOFROMEASY(fileElement, metaData, "bodyDemoteVolume", &bodyDemoteVolume);
        dataElement->setBodyDemoteVolume(bodyDemoteVolume);

        FloatP_t edgeSplitDist;
        TF_IOFROMEASY(fileElement, metaData, "edgeSplitDist", &edgeSplitDist);
        dataElement->setEdgeSplitDist(edgeSplitDist);

        // Reconnection knobs are read optionally (default if absent) so meshes saved by builds
        // before the native RNR port still deserialize.
        {
            ::TissueForge::io::IOChildMap _rcChildren = ::TissueForge::io::IOElement::children(fileElement);

            FloatP_t reconnectLength = 0.0;
            auto _rlItr = _rcChildren.find("reconnectLength");
            if(_rlItr != _rcChildren.end())
                ::TissueForge::io::fromFile(_rlItr->second, metaData, &reconnectLength);
            dataElement->setReconnectLength(reconnectLength);

            FloatP_t reconnectHysteresis = 0.0;
            auto _rhItr = _rcChildren.find("reconnectHysteresis");
            if(_rhItr != _rcChildren.end())
                ::TissueForge::io::fromFile(_rhItr->second, metaData, &reconnectHysteresis);
            dataElement->setReconnectHysteresis(reconnectHysteresis);

            bool reconnectEnergyGate = false;
            auto _rgItr = _rcChildren.find("reconnectEnergyGate");
            if(_rgItr != _rcChildren.end())
                ::TissueForge::io::fromFile(_rgItr->second, metaData, &reconnectEnergyGate);
            dataElement->setReconnectEnergyGate(reconnectEnergyGate);

            bool stockQualityOps = true;
            auto _sqItr = _rcChildren.find("stockQualityOps");
            if(_sqItr != _rcChildren.end())
                ::TissueForge::io::fromFile(_sqItr->second, metaData, &stockQualityOps);
            dataElement->setStockQualityOps(stockQualityOps);

            unsigned int reconnectInterval = 1;
            auto _riItr = _rcChildren.find("reconnectInterval");
            if(_riItr != _rcChildren.end())
                ::TissueForge::io::fromFile(_riItr->second, metaData, &reconnectInterval);
            dataElement->setReconnectInterval(reconnectInterval);
        }

        bool collision2D;
        TF_IOFROMEASY(fileElement, metaData, "collision2D", &collision2D);
        dataElement->setCollision2D(collision2D);

        if(FIO::hasImport() && TissueForge::models::vertex::io::VertexSolverFIOModule::hasImport()) {
            std::unordered_set<unsigned int> excludedVertices;
            TF_IOFROMEASY(fileElement, metaData, "excludedVertices", &excludedVertices);
            for(auto &oldId : excludedVertices) {
                auto id_itr = TissueForge::models::vertex::io::VertexSolverFIOModule::importSummary->vertexIdMap.find(oldId);
                if(id_itr != TissueForge::models::vertex::io::VertexSolverFIOModule::importSummary->vertexIdMap.end()) 
                    dataElement->excludeVertex(id_itr->second);
            }

            std::unordered_set<unsigned int> excludedSurfaces;
            TF_IOFROMEASY(fileElement, metaData, "excludedSurfaces", &excludedSurfaces);
            for(auto &oldId : excludedSurfaces) {
                auto id_itr = TissueForge::models::vertex::io::VertexSolverFIOModule::importSummary->surfaceIdMap.find(oldId);
                if(id_itr != TissueForge::models::vertex::io::VertexSolverFIOModule::importSummary->surfaceIdMap.end()) 
                    dataElement->excludeSurface(id_itr->second);
            }

            std::unordered_set<unsigned int> excludedBodies;
            TF_IOFROMEASY(fileElement, metaData, "excludedBodies", &excludedBodies);
            for(auto &oldId : excludedBodies) {
                auto id_itr = TissueForge::models::vertex::io::VertexSolverFIOModule::importSummary->bodyIdMap.find(oldId);
                if(id_itr != TissueForge::models::vertex::io::VertexSolverFIOModule::importSummary->bodyIdMap.end()) 
                    dataElement->excludeBody(id_itr->second);
            }
        }

        return S_OK;
    }
}

std::string TissueForge::models::vertex::MeshQuality::toString() {
    return TissueForge::io::toString(*this);
}

MeshQuality MeshQuality::fromString(const std::string &s) {
    return TissueForge::io::fromString<MeshQuality>(s);
}

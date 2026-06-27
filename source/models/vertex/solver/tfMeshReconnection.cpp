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
 * @file tfMeshReconnection.cpp
 *
 * Native 3D T1 / reversible network reconnection (RNR): the Okuda I<->H face<->edge swap. Holds the
 * read-only neighborhood walk + Condition-4 vetoes, the ReconnectionOperation (Okuda Appendix-1
 * placement + surface surgery), and the MeshQuality scheduler entry point + diagnostic methods.
 * Extracted verbatim from tfMeshQuality.cpp (no behavior change); the MeshQuality scheduler in that
 * file drives this pass via MeshReconnection_buildOperations().
 */

#include "tfMeshReconnection.h"

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
#include <cmath>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>


using namespace TissueForge;
using namespace TissueForge::models::vertex;


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

#ifdef TF_VERTEX_RNR_DEBUG
// --- JSON emit for the gated diagnostic/debug entry points (see tfMeshQuality.h) -------------

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

#endif // TF_VERTEX_RNR_DEBUG

// --- placement math (Okuda 2013 Appendix 1) -------------------------------------------------

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


/////////////////////////////////
// MeshQuality scheduler bridge //
/////////////////////////////////

// The reconnection-candidate scan, extracted from MeshQuality_constructOperationsReconnection so the
// RNR-specific logic (scan + Condition-4 veto + slot de-dup + op creation) lives with the rest of the
// RNR code. The generic dependency-chain build + serial execution stay in tfMeshQuality.cpp.
HRESULT TissueForge::models::vertex::MeshReconnection_buildOperations(
    Mesh *mesh,
    const std::vector<bool> &passMask,
    const FloatP_t &reconnectLength,
    const FloatP_t &reconnectHysteresis,
    const bool &reconnectEnergyGate,
    std::vector<MeshQualityOperation*> &ops
) {
    // ops is indexed by mesh-object id (here surface id, since both triggers live on a Surface: an
    // I->H short edge is a consecutive vertex pair on a surface, an H->I feature IS a triangular
    // surface). targets reference this same id space, so the dependency chains serialize ops with
    // overlapping touched surfaces.
    ops = std::vector<MeshQualityOperation*>(mesh->sizeSurfaces(), 0);

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

    return S_OK;
}


#ifdef TF_VERTEX_RNR_DEBUG

////////////////////////////////////////////////////
// MeshQuality RNR diagnostics + debug entry points //
////////////////////////////////////////////////////
//
// Gated behind TF_VERTEX_RNR_DEBUG (default ON in development builds; a knobs-only production build
// sets it OFF). These return JSON strings for Python-side introspection/testing; the
// force_reconnect_* pair bypasses the candidate scan and scheduling and can corrupt the mesh if
// misused, so it is not part of the stable API.

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

#endif // TF_VERTEX_RNR_DEBUG

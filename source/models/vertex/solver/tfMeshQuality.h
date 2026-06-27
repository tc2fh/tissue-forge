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
 * @file tfMeshQuality.h
 * 
 */

#ifndef _MODELS_VERTEX_SOLVER_TFMESHQUALITY_H_
#define _MODELS_VERTEX_SOLVER_TFMESHQUALITY_H_


#include "tfMeshObj.h"

#include <tf_port.h>

#include <mutex>
#include <set>
#include <unordered_set>
#include <vector>


namespace TissueForge::models::vertex {


    class Mesh;


    /**
     * @brief An operation that modifies the topology of a mesh to improve its quality
     */
    struct MeshQualityOperation {

        enum Flag : unsigned int {
            None    = 0, 
            Active  = 1 << 0, 
            Custom  = 1 << 1
        };

        unsigned int flags;

        /**
         * @brief Target mesh objects.
         * 
         * Used to identify dependencies between operations
         */
        std::vector<int> targets;

        /** Upstream operations, if any */
        std::set<MeshQualityOperation*> prev;

        /** Downstream operations, if any */
        std::set<MeshQualityOperation*> next;

        /** Lock, for safe modification during concurrent work */
        std::mutex lock;

        MeshQualityOperation(Mesh *_mesh);

        virtual ~MeshQualityOperation() {};

        /**
         * @brief Add an operation to the list of next operations
         * 
         * If the operation is already upstream of this operation, 
         * then the call is ignored
         * 
         * @param _next the operation
         */
        HRESULT appendNext(MeshQualityOperation *_next);

        /**
         * @brief Remove an operation to the list of next operations
         * 
         * @param _next the operation
         */
        HRESULT removeNext(MeshQualityOperation *_next);

        /**
         * @brief Compute all upstream operations
         */
        std::set<MeshQualityOperation*> upstreams() const;

        /**
         * @brief Compute all downstream operations
         */
        std::set<MeshQualityOperation*> downstreams() const;

        /**
         * @brief Compute all upstream operations that have no dependencies
         */
        std::set<MeshQualityOperation*> headOperations() const;

        /**
         * @brief Validate this operation
         */
        virtual HRESULT validate() { return S_OK; }

        /**
         * @brief Check whether this operation is still valid
         * 
         * @return true if this operation is still valid
         */
        virtual bool check() { return !(flags & Flag::None); };

        /**
         * @brief Do all prep, checks and planning for this operation
         */
        virtual void prep() {}

        /**
         * @brief Returns how many vertices will be created by this operation
         */
        virtual size_t numNewVertices() const { return 0; };

        /**
         * @brief Returns how many surfaces will be created by this operation
         */
        virtual size_t numNewSurfaces() const { return 0; };

        /**
         * @brief Returns how many bodies will be created by this operation
         */
        virtual size_t numNewBodies() const { return 0; };

        /**
         * @brief Implement this operation
         * 
         * @return ids of affected children, if any
         */
        virtual std::vector<int> implement() { return {}; }

    protected:

        Mesh *mesh;
    };


    /**
     * @brief Custom mesh quality operation.
     * 
     * todo: implement support for custom mesh quality operations
     */
    struct CAPI_EXPORT CustomQualityOperation : MeshQualityOperation {

        typedef bool (*OperationCheck)();
        typedef void (*OperationPrep)();
        typedef std::vector<int> (*OperationFunction)(std::vector<int>);

        OperationCheck *opCheck;
        OperationPrep *opPrep;
        OperationFunction *opFunc;

        CustomQualityOperation(Mesh *_mesh, OperationFunction *_opFunc, OperationCheck *_opCheck=NULL, OperationPrep *_opPrep=NULL) : 
            MeshQualityOperation(_mesh), 
            opFunc{_opFunc}, 
            opCheck{_opCheck}, 
            opPrep{_opPrep}
        {
            flags = Flag::Active | Flag::Custom;
        };

        virtual ~CustomQualityOperation() {
            if(opCheck) {
                delete opCheck;
                opCheck = 0;
            }
            delete opFunc;
            opFunc = 0;
        };

        bool check() override {
            if(opCheck) return (*opCheck)();
            return true;
        }

        void prep() override {
            if(opPrep) (*opPrep)();
        }

        std::vector<int> implement() override { return (*opFunc)(targets); }
    };


    /**
     * @brief An object that schedules topological operations on a mesh to maintain its quality
     * 
     */
    class CAPI_EXPORT MeshQuality {

        /** 
         * Vertex merge criterion.
         * 
         * Two vertices are merged if their distance is less than this value. 
         */
        FloatP_t vertexMergeDist;

        /**
         * Surface demotion criterion.
         * 
         * A surface becomes a vertex if its area is less than this value.
         */
        FloatP_t surfaceDemoteArea;

        /**
         * Body demotion criterion.
         * 
         * A body becomes a vertex if its volume is less than this value.
         */
        FloatP_t bodyDemoteVolume;

        /**
         * Initial length of an edge created by splitting a vertex.
         */
        FloatP_t edgeSplitDist;

        /**
         * Reconnection trigger length (Okuda Condition 2, Delta_l_th).
         *
         * Length below which the native 3D T1 / reversible network reconnection (RNR, the
         * Okuda I<->H face<->edge swap) fires: a short interior edge is reconnected (I->H)
         * when its length is below this value, and a small triangular face is reconnected
         * (H->I) when its longest edge is below it. This is an ABSOLUTE length (physical
         * units, like an edge length), NOT a box-fraction coefficient like vertexMergeDist.
         * 0 disables reconnection (the default), so the reconnection pass is a no-op until set.
         */
        FloatP_t reconnectLength;

        /**
         * Reconnection placement hysteresis.
         *
         * Features created by a reconnection are sized at reconnectLength*(1+reconnectHysteresis),
         * so a fresh feature sits above the trigger and does not immediately reconnect back
         * (an anti-thrash gap). 0 = faithful (infinitesimal features, as in Okuda).
         */
        FloatP_t reconnectHysteresis;

        /**
         * Optional energy gate on reconnection (OFF by default).
         *
         * When true, a reconnection that raises the local heterotypic interfacial energy is
         * rejected (greedy/Metropolis-at-T=0). This is a DEPARTURE from Okuda's purely
         * geometric Condition-2 trigger and is an instability driver; the faithful default
         * is false.
         */
        bool reconnectEnergyGate;

        /**
         * Master flag for whether the stock (non-RNR) quality operations are enabled.
         *
         * These are TissueForge's pre-existing vertex split/merge, surface/body demote, and 2D
         * collision repair passes. Enabled by default to preserve historical MeshQuality behavior.
         * When false, doQuality() runs only the native RNR (Okuda I<->H) reconnection pass. The
         * three per-pass flags below give finer control: a pass runs only if BOTH this master flag
         * and its per-pass flag are true. This lets the native RNR coexist with the stock passes
         * without an all-or-nothing switch (e.g. keep vertex/surface ops while excluding body
         * demotion, which is not finite-block-safe -- see stockBodyOps).
         */
        bool stockQualityOps;

        /** Whether the stock vertex pass (merge/split) runs (when stockQualityOps is also true). */
        bool stockVertexOps;

        /** Whether the stock surface pass (surface demote + 2D collision) runs (when stockQualityOps
         *  is also true). The collision sub-step is additionally gated by collision2D. */
        bool stockSurfaceOps;

        /**
         * Whether the stock body pass (body -> vertex demotion) runs (when stockQualityOps is also
         * true). Default true preserves historical behavior, BUT the body-demotion collapse is not
         * robust on finite (free-surface) blocks: when a cell's volume falls below bodyDemoteVolume
         * the degenerate body->vertex collapse (Vertex::replace) can cascade and crash. A vertex
         * model that drives topology change through the native RNR reconnection instead of stock
         * collapses should set this false to coexist safely with the stock vertex/surface passes.
         */
        bool stockBodyOps;

        /**
         * Reconnection throttle interval (the 3DVertVor oracle's dtr, in doQuality calls).
         *
         * The native reconnection pass runs only on every reconnectInterval-th doQuality() call
         * (i.e. every reconnectInterval-th integration step). 1 = every step (the original
         * behavior). Larger values give the mesh time to relax between reconnections, which breaks
         * the post-reconnection overshoot storm (the 3DVertVor oracle reconnects every dtr = 10*dt
         * rather than every step; see rnr/PORTING_NOTES.md sections 6c/6d). Trade-off: a larger
         * interval is more stable but reconnects (and therefore sorts) more slowly. Values < 1 are
         * treated as 1.
         */
        unsigned int reconnectInterval;

        /**
         * Flag for whether doing 2D collisions
         */
        bool collision2D;

        /** Flag for whether currently doing work */
        bool _working;

        /** doQuality() call counter; gates the reconnection pass by reconnectInterval (transient,
         * not serialized). */
        unsigned int reconnectCounter;

        std::unordered_set<unsigned int> excludedVertices;

        std::unordered_set<unsigned int> excludedSurfaces;

        std::unordered_set<unsigned int> excludedBodies;

    public:

        MeshQuality(
            const FloatP_t &vertexMergeDistCf=0.0001,
            const FloatP_t &surfaceDemoteAreaCf=0.0001,
            const FloatP_t &bodyDemoteVolumeCf=0.0001,
            const FloatP_t &_edgeSplitDistCf=2.0,
            const FloatP_t &_reconnectLength=0.0
        );

        /**
         * @brief Get a summary string
         */
        std::string str() const;

        /**
         * @brief Get a JSON string representation
         */
        std::string toString();

        /**
         * @brief Create an instance from a JSON string representation
         * 
         * @param s JSON string representation
         */
        static MeshQuality fromString(const std::string &s);

        /**
         * @brief Perform quality operations work
         */
        HRESULT doQuality();

        /**
         * @brief Test whether quality operations are being done
         * 
         * @return true if quality operations are being done
         */
        const bool working() const { return _working; }

        /**
         * @brief Get the distance below which two vertices are scheduled for merging
         */
        FloatP_t getVertexMergeDistance() const { return vertexMergeDist; };

        /**
         * @brief Set the distance below which two vertices are scheduled for merging
         * 
         * @param _val distance
         */
        HRESULT setVertexMergeDistance(const FloatP_t &_val);

        /**
         * @brief Get the area below which a surface is scheduled to become a vertex
         */
        FloatP_t getSurfaceDemoteArea() const { return surfaceDemoteArea; };

        /**
         * @brief Set the area below which a surface is scheduled to become a vertex
         * 
         * @param _val area
         */
        HRESULT setSurfaceDemoteArea(const FloatP_t &_val);

        /**
         * @brief Get the volume below which a body is scheduled to become a vertex
         */
        FloatP_t getBodyDemoteVolume() const { return bodyDemoteVolume; }

        /**
         * @brief Set the volume below which a body is scheduled to become a vertex
         * 
         * @param _val volume
         */
        HRESULT setBodyDemoteVolume(const FloatP_t &_val);

        /**
         * @brief Get the distance at which two vertices are seperated when a vertex is split
         */
        FloatP_t getEdgeSplitDist() const { return edgeSplitDist; };

        /**
         * @brief Set the distance at which two vertices are seperated when a vertex is split
         * 
         * @param _val distance
         */
        HRESULT setEdgeSplitDist(const FloatP_t &_val);

        /**
         * @brief Get the reconnection trigger length (Okuda Condition 2, Delta_l_th)
         */
        FloatP_t getReconnectLength() const { return reconnectLength; };

        /**
         * @brief Set the reconnection trigger length (Okuda Condition 2, Delta_l_th)
         *
         * @param _val length (absolute; 0 disables reconnection)
         */
        HRESULT setReconnectLength(const FloatP_t &_val);

        /**
         * @brief Get the reconnection placement hysteresis
         */
        FloatP_t getReconnectHysteresis() const { return reconnectHysteresis; };

        /**
         * @brief Set the reconnection placement hysteresis
         *
         * @param _val hysteresis (>= 0)
         */
        HRESULT setReconnectHysteresis(const FloatP_t &_val);

        /**
         * @brief Get whether the optional reconnection energy gate is enabled
         */
        bool getReconnectEnergyGate() const { return reconnectEnergyGate; };

        /**
         * @brief Set whether the optional reconnection energy gate is enabled
         *
         * @param _val flag (a DEPARTURE from Okuda's geometric trigger; default false)
         */
        HRESULT setReconnectEnergyGate(const bool &_val);

        /**
         * @brief Get whether stock non-RNR quality operations are enabled
         */
        bool getStockQualityOps() const { return stockQualityOps; };

        /**
         * @brief Set whether stock non-RNR quality operations are enabled
         *
         * When false, doQuality() skips the legacy vertex/surface/body/collision passes and
         * runs only the native reconnection pass (if reconnectLength > 0).
         *
         * @param _val flag
         */
        HRESULT setStockQualityOps(const bool &_val);

        /** @brief Get whether the stock vertex pass (merge/split) runs (with stockQualityOps). */
        bool getStockVertexOps() const { return stockVertexOps; };

        /** @brief Set whether the stock vertex pass (merge/split) runs (with stockQualityOps). */
        HRESULT setStockVertexOps(const bool &_val);

        /** @brief Get whether the stock surface pass (demote + 2D collision) runs (with stockQualityOps). */
        bool getStockSurfaceOps() const { return stockSurfaceOps; };

        /** @brief Set whether the stock surface pass (demote + 2D collision) runs (with stockQualityOps). */
        HRESULT setStockSurfaceOps(const bool &_val);

        /** @brief Get whether the stock body pass (body->vertex demotion) runs (with stockQualityOps). */
        bool getStockBodyOps() const { return stockBodyOps; };

        /**
         * @brief Set whether the stock body pass (body->vertex demotion) runs (with stockQualityOps).
         *
         * Default true preserves historical behavior; set false to coexist safely with the native
         * RNR pass on finite blocks (the body-demotion collapse is not finite-block-safe).
         *
         * @param _val flag
         */
        HRESULT setStockBodyOps(const bool &_val);

        /**
         * @brief Get the reconnection throttle interval (doQuality calls between reconnection passes)
         */
        unsigned int getReconnectInterval() const { return reconnectInterval; };

        /**
         * @brief Set the reconnection throttle interval (the oracle's dtr, in doQuality calls)
         *
         * The reconnection pass runs only every _val-th doQuality() call; larger values relax the
         * mesh between reconnections and break the post-reconnection overshoot storm. Values < 1
         * are treated as 1 (run every step).
         *
         * @param _val interval (>= 1)
         */
        HRESULT setReconnectInterval(const unsigned int &_val);

#ifdef TF_VERTEX_RNR_DEBUG
        // --- RNR diagnostic / debug entry points (NOT part of the stable API) ----------------
        //
        // Exposed only when built with TF_VERTEX_RNR_DEBUG (the default in development builds; a
        // knobs-only production build sets it OFF and omits these). They return JSON object/array
        // strings for Python-side introspection and testing. The forceReconnect* pair bypasses the
        // candidate scan and all scheduling/safety, so it can corrupt the mesh if misused.

        /**
         * @brief Diagnostic (read-only): analyze the I->H reconnection neighborhood of the short
         *        edge (v10Id, v11Id) on the current mesh. Returns a JSON object string
         *        {valid, kind, v10_id, v11_id, cap_top_id, cap_bot_id, side_cell_ids, length,
         *        legal, veto_reason}. Does not mutate the mesh.
         */
        std::string analyzeIReconnection(const unsigned int &v10Id, const unsigned int &v11Id) const;

        /**
         * @brief Diagnostic (read-only): analyze the H->I reconnection neighborhood of the
         *        triangular surface triId. Returns a JSON object string; does not mutate the mesh.
         */
        std::string analyzeHReconnection(const unsigned int &triId) const;

        /**
         * @brief Diagnostic (read-only): JSON array of every reconnection candidate the native scan
         *        finds at the current reconnectLength (Okuda Condition 2) -- both I->H short edges
         *        and H->I small triangles -- each tagged legal + veto_reason. Returns "[]" when
         *        reconnectLength <= 0. Same scanners doQuality uses; does not mutate the mesh.
         */
        std::string findReconnectionCandidates() const;

        /**
         * @brief Debug entry point: force one native I->H reconnection on the current mesh,
         *        bypassing the scan and stock quality passes. Uses reconnectLength as the Okuda
         *        placement length. Returns a JSON object string {ok, reason, new_surface_id,
         *        new_vertex_ids}. Mutates the mesh; not part of the stable API.
         */
        std::string forceReconnectIToH(const unsigned int &v10Id, const unsigned int &v11Id) const;

        /**
         * @brief Debug entry point: force one native H->I reconnection on the current mesh,
         *        bypassing the scan and stock quality passes. Uses reconnectLength as the Okuda
         *        placement length. Returns a JSON object string {ok, reason, new_surface_id,
         *        new_vertex_ids}. Mutates the mesh; not part of the stable API.
         */
        std::string forceReconnectHToI(const unsigned int &triId) const;
#endif // TF_VERTEX_RNR_DEBUG

        /**
         * @brief Get whether 2D collisions are implemented
         *
         * @return true if 2D collisions are implemented
         */
        bool getCollision2D() const { return collision2D; }

        /**
         * @brief Set whether 2D collisions are implemented
         * 
         * @param _collision2D flag indicating whether 2D collisions are implemented
         */
        HRESULT setCollision2D(const bool &_collision2D);

        /**
         * @brief Exclude a vertex from quality operations
         * 
         * @param id vertex id
         */
        HRESULT excludeVertex(const unsigned int &id);

        /**
         * @brief Exclude a surface from quality operations
         * 
         * @param id surface id
         */
        HRESULT excludeSurface(const unsigned int &id);

        /**
         * @brief Exclude a body from quality operations
         * 
         * @param id body id
         */
        HRESULT excludeBody(const unsigned int &id);

        /**
         * @brief Include a vertex from quality operations. 
         * 
         * Unless otherwise specified, all vertices are subject to operations.
         * 
         * @param id vertex id
         */
        HRESULT includeVertex(const unsigned int &id);

        /**
         * @brief Include a surface from quality operations. 
         * 
         * Unless otherwise specified, all surfaces are subject to operations.
         * 
         * @param id surface id
         */
        HRESULT includeSurface(const unsigned int &id);

        /**
         * @brief Include a body from quality operations. 
         * 
         * Unless otherwise specified, all bodies are subject to operations.
         * 
         * @param id body id
         */
        HRESULT includeBody(const unsigned int &id);

        /**
         * @brief Get the current set of excluded vertices
         */
        std::unordered_set<unsigned int> getExcludedVertices() const { return excludedVertices; }

        /**
         * @brief Get the current set of excluded surfaces
         */
        std::unordered_set<unsigned int> getExcludedSurfaces() const { return excludedSurfaces; }

        /**
         * @brief Get the current set of excluded bodies
         */
        std::unordered_set<unsigned int> getExcludedBodies() const { return excludedBodies; }
    };
}


#endif // _MODELS_VERTEX_SOLVER_TFMESHQUALITY_H_

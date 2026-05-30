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

%{

#include <models/vertex/solver/tfMeshQuality.h>
%}

%ignore TissueForge::models::vertex::MeshQualityOperation;
%ignore TissueForge::models::vertex::CustomQualityOperation;

%rename(do_quality) TissueForge::models::vertex::MeshQuality::doQuality;
%rename(exclude_vertex) TissueForge::models::vertex::MeshQuality::excludeVertex;
%rename(exclude_surface) TissueForge::models::vertex::MeshQuality::excludeSurface;
%rename(exclude_body) TissueForge::models::vertex::MeshQuality::excludeBody;
%rename(include_vertex) TissueForge::models::vertex::MeshQuality::includeVertex;
%rename(include_surface) TissueForge::models::vertex::MeshQuality::includeSurface;
%rename(include_body) TissueForge::models::vertex::MeshQuality::includeBody;

%rename(_vertex_solver_Quality) TissueForge::models::vertex::MeshQuality;

%include <models/vertex/solver/tfMeshQuality.h>

%extend TissueForge::models::vertex::MeshQuality {
    %pythoncode %{
        @property
        def vertex_merge_distance(self) -> float:
            """Distance below which two vertices are scheduled for merging"""
            return self.getVertexMergeDistance()

        @vertex_merge_distance.setter
        def vertex_merge_distance(self, _val: float):
            self.setVertexMergeDistance(_val)

        @property
        def surface_demote_area(self) -> float:
            """Area below which a surface is scheduled to become a vertex"""
            return self.getSurfaceDemoteArea()

        @surface_demote_area.setter
        def surface_demote_area(self, _val: float):
            self.setSurfaceDemoteArea(_val)

        @property
        def body_demote_volume(self) -> float:
            """Volume below which a body is scheduled to become a vertex"""
            return self.getBodyDemoteVolume()

        @body_demote_volume.setter
        def body_demote_volume(self, _val: float):
            self.setBodyDemoteVolume(_val)

        @property
        def edge_split_distance(self) -> float:
            """Distance at which two vertices are seperated when a vertex is split"""
            return self.getEdgeSplitDist()

        @edge_split_distance.setter
        def edge_split_distance(self, _val: float):
            self.setEdgeSplitDist(_val)

        @property
        def reconnect_length(self) -> float:
            """Length below which a short interior edge / small triangular face is reconnected
            (native 3D T1 / Okuda I<->H RNR; Okuda Condition 2, Delta_l_th). Absolute length;
            0 disables reconnection (the default)."""
            return self.getReconnectLength()

        @reconnect_length.setter
        def reconnect_length(self, _val: float):
            self.setReconnectLength(_val)

        @property
        def reconnect_hysteresis(self) -> float:
            """Placement hysteresis: features created by a reconnection are sized at
            reconnect_length*(1+reconnect_hysteresis) (anti-thrash gap). 0 = faithful."""
            return self.getReconnectHysteresis()

        @reconnect_hysteresis.setter
        def reconnect_hysteresis(self, _val: float):
            self.setReconnectHysteresis(_val)

        @property
        def reconnect_energy_gate(self) -> bool:
            """Whether to reject reconnections that raise local heterotypic energy (a DEPARTURE
            from Okuda's geometric trigger; an instability driver). Default False."""
            return self.getReconnectEnergyGate()

        @reconnect_energy_gate.setter
        def reconnect_energy_gate(self, _val: bool):
            self.setReconnectEnergyGate(_val)

        def analyze_i_reconnection(self, v10_id: int, v11_id: int) -> dict:
            """Diagnostic (read-only): the native I->H neighborhood walk + Condition-4 veto for
            the short edge (v10_id, v11_id) on the current mesh -- the C++ port of
            topology.i_neighbourhood + conditions.i_to_h_veto. Returns a dict: keys `valid`
            (bool) and, when valid, `kind` ('I'), `v10_id`, `v11_id`, `triangle_id` (-1),
            `trigger_surface_id`, `cap_top_id`, `cap_bot_id`, `side_cell_ids` (list), `length`,
            `legal` (bool), `veto_reason` (str, '' if legal). For cross-checking against the
            Python prototype (Phase-B gate); does not mutate the mesh."""
            import json
            return json.loads(self.analyzeIReconnection(v10_id, v11_id))

        def analyze_h_reconnection(self, triangle_id: int) -> dict:
            """Diagnostic (read-only): the native H->I neighborhood walk + Condition-4 veto for
            the triangular surface `triangle_id` (port of h_neighbourhood + h_to_i_veto). Returns
            a dict like analyze_i_reconnection but with `kind` ('H'), `triangle_id`, and `length`
            = the triangle's max edge. Does not mutate the mesh."""
            import json
            return json.loads(self.analyzeHReconnection(triangle_id))

        def find_reconnection_candidates(self) -> list:
            """Diagnostic (read-only): every reconnection candidate the native scan finds at the
            current `reconnect_length` (Okuda Condition 2) -- both I->H short edges and H->I small
            triangles -- each as a dict tagged `legal` + `veto_reason`. Empty list when
            reconnect_length <= 0. Same scanners doQuality uses; does not mutate the mesh."""
            import json
            return json.loads(self.findReconnectionCandidates())

        def force_reconnect_i_to_h(self, v10_id: int, v11_id: int) -> dict:
            """Debug/test mutate entry point: force one native I->H reconnection on the current
            mesh, using `reconnect_length` as Okuda Delta_l_th. Bypasses the scan and stock
            quality passes; returns `{ok, reason, new_surface_id, new_vertex_ids}`."""
            import json
            return json.loads(self.forceReconnectIToH(v10_id, v11_id))

        def force_reconnect_h_to_i(self, triangle_id: int) -> dict:
            """Debug/test mutate entry point: force one native H->I reconnection on the current
            mesh, using `reconnect_length` as Okuda Delta_l_th. Bypasses the scan and stock
            quality passes; returns `{ok, reason, new_surface_id, new_vertex_ids}`."""
            import json
            return json.loads(self.forceReconnectHToI(triangle_id))

        @property
        def collision_2d(self) -> bool:
            """Whether 2D collisions are implemented"""
            return self.getCollision2D()

        @collision_2d.setter
        def collision_2d(self, _collision_2d: bool):
            self.setCollision2D(_collision_2d)

        def __str__(self) -> str:
            return self.str()
    %}
}

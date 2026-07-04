# ******************************************************************************
# This file is part of Tissue Forge.
# Copyright (c) 2022-2024 T.J. Sego and Tien Comlekoglu
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# ******************************************************************************

"""
Three-dimensional cell sorting in the vertex model, driven by a native 3D T1.

Two cell types carry a heterotypic interfacial tension. Unlike the 2D ``cell_sorting`` example --
which reaches neighbor exchange through the 2D T1 transition (``vertex_merge`` + ``edge_split``) --
a 3D vertex model requires a genuine 3D T1: the reversible network reconnection (RNR) / I<->H
operation of

    Okuda et al. (2013) "Reversible network reconnection model for simulating large
    deformation in dynamic tissue morphogenesis." Biomech. Model. Mechanobiol. 12:627-644.

This example demonstrates that operation running composably inside ``tf.step()``: it is provided
natively by the ``MeshQuality`` reconnection pass (enabled with the single knob
``reconnect_length``), running alongside the safe stock quality passes -- no global kill-switch.
Cells are driven by an active self-propulsion motility drive (``MeshSolver.set_motility``); the
heterotypic tension biases the reconnections, the mechanism that sorts the tissue.

The mesh is a cubic lattice of cells (``create_plpd_mesh``); an interior edge is shared by four
cells -- the textbook Okuda [I] configuration -- so RNR triggers as the tension deforms the mesh.
The run reports the heterotypic-interface fraction (the sorting order parameter) and the running
reconnection count. Note this small, free-surface block illustrates the *machinery*; pronounced
macroscopic demixing emerges in a larger, periodic bulk run over many more steps (see
Zhang & Schwarz, Phys. Rev. Research 4:043148, whose Fig. 1E/1F this capability reproduces).

Runs headless (no rendering required). Set ``WITH_WINDOW = True`` to watch it live.
"""

import numpy as np
import tissue_forge as tf
from tissue_forge.models.vertex import solver as tfv
from tissue_forge.models.vertex.solver.mesh_types import BodyTypeSpec, SurfaceTypeSpec

WITH_WINDOW = False

# Interfacial tension between UNLIKE cells (homotypic = 0). Larger -> stronger sorting drive.
sigma = 0.5
# RNR trigger: reconnect an edge once it shrinks below this length (Okuda's Delta l_th).
reconnect_length = 0.2
# Active self-propulsion speed and director rotational-diffusion rate (the sorting drive).
motility_v0 = 0.2
motility_dr = 1.0
n_per_axis = 4          # 4 x 4 x 4 = 64 cells
cell_len = 1.0
n_steps = 2000
seed = 5

tf.init(windowless=not WITH_WINDOW, dim=[20., 20., 20.], cutoff=3.0, dt=0.01)
tfv.init()

mesh = tfv.MeshSolver.get().get_mesh()
mesh.quality = None


class Interface(SurfaceTypeSpec):
    pass


class TypeA(BodyTypeSpec):
    volume_lam = 5.0
    volume_val = 1.0
    surface_area_lam = 1.0
    surface_area_val = 6.0
    adhesion = {"TypeA": 0.0, "TypeB": sigma}


class TypeB(BodyTypeSpec):
    volume_lam = 5.0
    volume_val = 1.0
    surface_area_lam = 1.0
    surface_area_val = 6.0
    adhesion = {"TypeA": sigma, "TypeB": 0.0}


stype, btype_a, btype_b = Interface.get(), TypeA.get(), TypeB.get()
BodyTypeSpec.bind_adhesion([TypeA, TypeB])

# Cubic lattice of cells (all TypeA to start).
grid = tfv.create_plpd_mesh(btype_a, stype, tf.FVector3(3., 3., 3.),
                            n_per_axis, n_per_axis, n_per_axis,
                            cell_len, cell_len, cell_len)
bodies = [b for plane in grid for row in plane for b in row]

# Randomly relabel half the cells TypeB (the mixed initial condition).
rng = np.random.default_rng(seed)
for b in bodies:
    if rng.random() < 0.5:
        b.become(btype_b)

# Enable the native RNR reconnection alongside the safe stock quality passes. The stock
# body-demote pass is excluded (per-pass flag) because it is not robust on this finite,
# free-surface block; the vertex/surface passes stay on.
q = tfv.Quality()
q.stock_quality_operations = True
q.stock_body_operations = False
q.reconnect_length = reconnect_length
q.reconnect_hysteresis = 0.2
mesh.quality = q

# Native active-motility drive: per-cell directors + a per-vertex active force, all in C++.
tfv.MeshSolver.set_motility(motility_v0, motility_dr, seed + 1)


def heterotypic_fraction():
    """Fraction of cell-cell interfaces that separate UNLIKE cells (falls as cells sort)."""
    het = tot = 0
    for b in bodies:
        b_type = b.type().id
        for nb in b.connected_bodies:
            if nb.id <= b.id:
                continue
            tot += 1
            if nb.type().id != b_type:
                het += 1
    return het / tot if tot else 0.0


_recon = {"count": 0, "nv": mesh.num_vertices}


def track_reconnections():
    nv = mesh.num_vertices
    if nv != _recon["nv"]:
        _recon["count"] += abs(nv - _recon["nv"])
        _recon["nv"] = nv


def report(step):
    vols = [b.volume for b in bodies]
    print(f"step {step:5d}: heterotypic_fraction={heterotypic_fraction():.3f} "
          f"reconnections={_recon['count']:5d} "
          f"min_vol={min(vols):.3f} max_vol={max(vols):.3f}", flush=True)


if WITH_WINDOW:
    tf.run()
else:
    report(0)
    for step in range(1, n_steps + 1):
        tf.step()
        track_reconnections()
        if step % 200 == 0:
            report(step)

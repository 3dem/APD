#!/usr/bin/env python

import sys
import os
import math
from collections import defaultdict, deque
from Bio.PDB import PDBParser, MMCIFParser, NeighborSearch, is_aa

# contacts.py
#
# This script:
#   1) Identifies protofilaments based on CA-CA stacking.
#   2) Finds the "middle chain" (chain_offset == 0) in each protofilament.
#   3) For residues in the middle chains, finds contacts between:
#        - side-chain atoms (non-backbone, non-hydrogen), and
#        - CA atoms of GLY residues,
#      within sidechain_cutoff.
#   4) Applies an "outer side" filter per residue and per contact segment:
#        - For each residue i in a protofilament, a local plane is defined
#          using:
#              N(i) in a reference chain of that protofilament,
#              N(i) in a neighboring amyloid layer,
#              C(i) in the same reference chain.
#        - For each protofilament, every chain picks a single neighbor chain:
#              the chain in the same protofilament with the smallest
#              absolute difference in chain_offset (not equal to itself).
#        - The plane normal is stored together with N(i) as a reference
#          point for that residue.
#        - The sidechain "direction" is approximated by the position of
#          CB(i).
#        - A contact segment between two atoms is kept only if both
#          endpoints lie on the same side of this plane as CB(i), for
#          each residue where this plane can be defined; otherwise it is
#          discarded.
#        - If the plane or CB is missing, that residue does not veto the
#          contact (the contact is not discarded on that basis).
#   5) Excludes local sequence neighbors along the protofilament axis:
#        - For each protofilament, each chain gets a "sequence index"
#          per residue in PDB order. The first residue of the chain gets
#          index 0. For each subsequent residue:
#              * Let delta_resnum be the absolute difference in residue
#                numbers between this residue and the previous one.
#              * If delta_resnum == 1, the sequence index increases by 1.
#              * If delta_resnum != 1, and both residues have CA atoms
#                and their CA-CA distance is < max_backbone_gap_distance,
#                the sequence index also increases by 1 (non-consecutive
#                numbering is tolerated because the backbone is continuous).
#              * Otherwise the sequence index increases by delta_resnum.
#        - For each chain, an integer axial offset is computed relative
#          to the middle chain in that protofilament using overlapping
#          residue numbers and these sequence indices.
#        - The axial index of a residue is:
#              axial_index = sequence_index_in_chain + chain_axial_offset
#        - Intra-protofilament contacts (same pf) are excluded when the
#          absolute difference between these axial indices is less than
#          min_sequence_separation.
#   6) For each residue pair, stores:
#        - shortest atom-atom distance among considered atoms (distance)
#        - coordinates of the closest atom pair.
#        - contact is stored only once per unordered residue pair when
#          both residues correspond to middle-chain residues in their
#          protofilaments, using an ordering on (pf, residue_number).
#        - multi-layer contacts (different chain offsets for residue2)
#          are collapsed to the shortest distance.
#   7) For each residue in the middle chain of a protofilament, computes
#      a handedness value based on the left/right position of the side
#      chain in a CA-based local frame.
#   8) Computes a helical rise per protofilament from the average CA-CA
#      Z-distance between identical residues in layers 0 and 1.
#   9) For each stored contact, computes a real-valued ca_offset for
#      residue2 as a fractional number of rises based on the vertical
#      (Z) separation between the middle-layer copies of residue2 and
#      residue1:
#
#          ca_offset = (z_mid_partner(residue2) - z_mid_src(residue1)) / rise
#
#      where z_mid_src(residue1) is the CA z-coordinate of residue1 in
#      the middle chain of its protofilament, and z_mid_partner(residue2)
#      is the CA z-coordinate of residue2 in the middle chain of its
#      protofilament. The rise used is the average of the available
#      per-protofilament rises for the two protofilaments involved in
#      the contact (or the single available rise if only one is defined).
#
# Output:
#   <root>_contacts.csv with columns:
#     pf,
#     residue,resname,source_chain,handedness,
#     ca_x,ca_y,ca_z,
#     cb_x,cb_y,cb_z,
#     partner_pf,partner_chain_offset,ca_offset,partner_residue,partner_resname,partner_chain,
#     contact_type,distance
#
#   <root>_contacts.bild
#

# Parameters

# Side-chain contact cutoff (Angstrom), used as relaxed threshold in compare.py
sidechain_cutoff = 6.5

# Minimum sequence separation (in axial index units) for
# intra-protofilament contacts. Residues whose axial indices
# differ by less than this value in the same protofilament are
# excluded from intra-protofilament contacts.
min_sequence_separation = 3

# Parameters for protofilament detection based on CA-CA stacking
stacking_ca_cutoff = 5.5
min_overlap_residues = 5
min_stack_fraction = 0.8

# BILD cylinder radius for contacts
cylinder_radius = 0.4

# Tolerance used in the "outer side" test. If the dot product involved
# in classifying "side of plane" has absolute value below this, we treat
# it as ambiguous and do not reject the contact on that basis.
outer_side_eps = 1.0e-3

# Maximum CA-CA distance to allow treating non-consecutive residue
# numbering as if it were consecutive for the sequence index.
max_backbone_gap_distance = 4.2  # Angstrom

# 0) Command line and file names

if len(sys.argv) != 2:
    print("Usage: python contacts.py <input.pdb|input.cif>")
    sys.exit(1)

infile = sys.argv[1]
if not os.path.isfile(infile):
    print("Error: file not found:", infile)
    sys.exit(1)

root = os.path.splitext(infile)[0]
csv_file = root + "_contacts.csv"
bild_file = root + "_contacts.bild"

# 1) Parse structure

ext = os.path.splitext(infile)[1].lower()
if ext in [".cif", ".mmcif"]:
    parser = MMCIFParser(QUIET=True)
else:
    parser = PDBParser(QUIET=True)

structure = parser.get_structure("prot", infile)

# 2) Work on first model

model = next(structure.get_models())

# 3) Collect chains and CA atoms for protofilament detection

protein_chains = []
ca_maps = {}  # chain_id -> {resnum: CA atom}

for chain in model:
    ca_map = {}
    has_protein = False
    for res in chain:
        if not is_aa(res, standard=True):
            continue
        has_protein = True
        if "CA" in res:
            resnum = res.id[1]
            ca_map[resnum] = res["CA"]
    if has_protein and ca_map:
        protein_chains.append(chain)
        ca_maps[chain.id] = ca_map

if not protein_chains:
    raise RuntimeError("No protein chains with CA atoms found.")

chain_ids = [c.id for c in protein_chains]
chain_by_id = {c.id: c for c in protein_chains}


def chain_center(cid):
    xs, ys, zs = [], [], []
    for r, atom in ca_maps[cid].items():
        xs.append(float(atom.coord[0]))
        ys.append(float(atom.coord[1]))
        zs.append(float(atom.coord[2]))
    n = float(len(xs))
    return (sum(xs) / n, sum(ys) / n, sum(zs) / n)


def ca_distance(a1, a2):
    dv = a1.coord - a2.coord
    return math.sqrt(float(dv[0] ** 2 + dv[1] ** 2 + dv[2] ** 2))


# 4) Build graph of stacked chains using CA-CA distances

from collections import defaultdict as _dd

stack_graph = _dd(set)

for i in range(len(chain_ids)):
    for j in range(i + 1, len(chain_ids)):
        c1 = chain_ids[i]
        c2 = chain_ids[j]
        cm1 = ca_maps[c1]
        cm2 = ca_maps[c2]

        common = set(cm1.keys()) & set(cm2.keys())
        if len(common) < min_overlap_residues:
            continue

        close = 0
        tot = 0
        for r in common:
            if ca_distance(cm1[r], cm2[r]) <= stacking_ca_cutoff:
                close += 1
            tot += 1

        if tot > 0 and (close / float(tot)) >= min_stack_fraction:
            stack_graph[c1].add(c2)
            stack_graph[c2].add(c1)

# 5) Identify protofilaments and order chains

def bfs(start):
    visited = set([start])
    q = deque([start])
    while q:
        u = q.popleft()
        for v in stack_graph[u]:
            if v not in visited:
                visited.add(v)
                q.append(v)
    return visited

# Ensure all chains appear in stack_graph
for cid in chain_ids:
    if cid not in stack_graph:
        stack_graph[cid] = set()

unseen = set(chain_ids)
protofilaments = []

while unseen:
    start = next(iter(unseen))
    comp = bfs(start)
    unseen -= comp

    deg = {c: len(stack_graph[c]) for c in comp}
    ends = [c for c in comp if deg[c] == 1]
    if ends:
        start_c = ends[0]
    else:
        start_c = next(iter(comp))

    ordered = []
    visited = set([start_c])
    cur = start_c

    while True:
        ordered.append(cur)
        nxt = [n for n in stack_graph[cur] if n not in visited]
        if not nxt:
            break
        cur = nxt[0]
        visited.add(cur)

    protofilaments.append(ordered)

# 6) Align protofilaments to same direction

ref_axis = None

for pf in protofilaments:
    if len(pf) < 2:
        continue
    a = chain_center(pf[0])
    b = chain_center(pf[-1])
    ax = (b[0] - a[0], b[1] - a[1], b[2] - a[2])
    norm = math.sqrt(sum(x * x for x in ax))
    if norm > 1.0e-6:
        ref_axis = (ax[0] / norm, ax[1] / norm, ax[2] / norm)
        break

if ref_axis is not None:
    for k, pf in enumerate(protofilaments):
        if len(pf) < 2:
            continue
        a = chain_center(pf[0])
        b = chain_center(pf[-1])
        ax = (b[0] - a[0], b[1] - a[1], b[2] - a[2])
        norm = math.sqrt(sum(x * x for x in ax))
        if norm < 1.0e-6:
            continue
        axn = (ax[0] / norm, ax[1] / norm, ax[2] / norm)
        dot = axn[0] * ref_axis[0] + axn[1] * ref_axis[1] + axn[2] * ref_axis[2]
        if dot < 0:
            protofilaments[k] = list(reversed(pf))

# 7) Assign chain offsets (middle = 0) and identify middle chains

chain_to_pf = {}
chain_offset = {}
pf_to_middle_chain_id = {}

for pf_index, pf in enumerate(protofilaments):
    mid = len(pf) // 2
    for pos, cid in enumerate(pf):
        chain_to_pf[cid] = pf_index
        chain_offset[cid] = pos - mid
    if pf:
        pf_to_middle_chain_id[pf_index] = pf[mid]

# Check that each protofilament has at least two amyloid layers
for pf_index, pf in enumerate(protofilaments):
    if len(pf) < 2:
        raise RuntimeError(
            "Protofilament %d has only one amyloid layer. "
            "Need at least two stacked chains to define handedness." % pf_index
        )

middle_chain_ids = []
for pf_index in sorted(pf_to_middle_chain_id.keys()):
    middle_chain_ids.append(pf_to_middle_chain_id[pf_index])

middle_chains = [chain_by_id[cid] for cid in middle_chain_ids]

print("Protofilaments:")
for idx, pf in enumerate(protofilaments):
    m = len(pf) // 2
    s = " ".join("%s(%d)" % (cid, pos - m) for pos, cid in enumerate(pf))
    print("PF", idx, ":", s)
print()

# Check that stacking is approximately along the Z axis.
print("Checking that protofilaments are approximately aligned with Z:")
for pf_index, pf in enumerate(protofilaments):
    if not pf:
        continue
    mid_cid = pf_to_middle_chain_id.get(pf_index, None)
    if mid_cid is None:
        continue

    neighbor_cid = None
    # Prefer offset +1 over -1
    for cid in pf:
        if chain_offset.get(cid) == 1:
            neighbor_cid = cid
            break
    if neighbor_cid is None:
        for cid in pf:
            if chain_offset.get(cid) == -1:
                neighbor_cid = cid
                break
    if neighbor_cid is None:
        continue

    c0 = chain_center(mid_cid)
    c1 = chain_center(neighbor_cid)
    dx = c1[0] - c0[0]
    dy = c1[1] - c0[1]
    dz = c1[2] - c0[2]
    radial = math.sqrt(dx * dx + dy * dy)
    axial = abs(dz)

    if axial < 1.0:
        print(
            "  Warning: PF %d has very small Z separation between chains "
            "(dz=%.3f). Assumed Z axis may be invalid." % (pf_index, dz)
        )
        continue

    ratio = radial / axial
    print(
        "  PF %d: delta center (dx,dy,dz)=(%.3f,%.3f,%.3f), radial/axial=%.3f"
        % (pf_index, dx, dy, dz, ratio)
    )
    if ratio > 0.2:
        print(
            "    Warning: PF %d may not be aligned with Z (radial/axial > 0.2)."
            % pf_index
        )
print()

# 8) Build residue lists, contact atoms, CA and CB positions

protein_residues = []
chain_residues = defaultdict(list)

sidechain_atoms_all = []           # all contact atoms (sidechain heavy atoms + GLY CA)
res_sidechain_atoms = defaultdict(list)

res_handedness = {}                # residue -> "R", "L", or "NA"
res_ca_coords = {}                 # residue -> (x, y, z) or None
res_cb_coords = {}                 # residue -> (x, y, z) or None (GLY or missing CB -> None)

# For the outer-side plane: residue -> {"normal": (nx,ny,nz), "point": (x,y,z)}
# The "point" is N(i) in the chain used as reference for that residue.
res_plane_data = {}

mainchain_names = {"N", "CA", "C", "O"}


def atom_distance(a, b):
    dv = a.coord - b.coord
    return math.sqrt(float(dv[0] ** 2 + dv[1] ** 2 + dv[2] ** 2))


def get_ca_atom(res):
    try:
        return res["CA"]
    except KeyError:
        return None


for chain in protein_chains:
    for res in chain:
        if not is_aa(res, standard=True):
            continue

        protein_residues.append(res)
        chain_residues[chain.id].append(res)

        ca_atom = get_ca_atom(res)
        if ca_atom is not None:
            res_ca_coords[res] = (
                float(ca_atom.coord[0]),
                float(ca_atom.coord[1]),
                float(ca_atom.coord[2]),
            )
        else:
            res_ca_coords[res] = None

        resname_upper = res.get_resname().upper()
        cb_coord_tuple = None
        if resname_upper != "GLY" and "CB" in res:
            cb_atom = res["CB"]
            cb_coord_tuple = (
                float(cb_atom.coord[0]),
                float(cb_atom.coord[1]),
                float(cb_atom.coord[2]),
            )
        res_cb_coords[res] = cb_coord_tuple

        for a in res:
            aname = a.get_name().strip()

            # Ignore hydrogens (any atom whose name starts with H)
            if aname.startswith("H"):
                continue

            # Ignore backbone atoms, except allow CA of GLY
            if aname in mainchain_names:
                if not (resname_upper == "GLY" and aname == "CA"):
                    continue

            sidechain_atoms_all.append(a)
            res_sidechain_atoms[res].append(a)

if not sidechain_atoms_all:
    raise RuntimeError("No sidechain heavy atoms found for contact detection.")

# Build mapping of residue numbers in middle chains for each protofilament.
pf_mid_resnums = {}
for pf_index, pf in enumerate(protofilaments):
    mid_cid = pf_to_middle_chain_id.get(pf_index, None)
    if mid_cid is None:
        continue
    pf_mid_resnums[pf_index] = set(res.id[1] for res in chain_residues[mid_cid])

# Estimate the helical rise per protofilament as the average Z-distance between CA atoms
# between identical residues in layers 0 and 1 (middle chain
# offset 0 and chain offset +1 if present, otherwise -1).
pf_rise = {}
print("Estimating helical rise from CA Z-distances between layers 0 and 1:")
for pf_index, pf in enumerate(protofilaments):
    if not pf:
        continue
    mid_cid = pf_to_middle_chain_id.get(pf_index, None)
    if mid_cid is None:
        continue

    # Choose layer 1 chain: prefer offset +1, otherwise -1
    layer1_cid = None
    for cid in pf:
        if chain_offset.get(cid) == 1:
            layer1_cid = cid
            break
    if layer1_cid is None:
        for cid in pf:
            if chain_offset.get(cid) == -1:
                layer1_cid = cid
                break

    if layer1_cid is None:
        print(
            "  PF %d: cannot find chain with offset +1 or -1; "
            "rise will be NA." % pf_index
        )
        continue

    mid_res_list = chain_residues[mid_cid]
    layer1_res_list = chain_residues[layer1_cid]
    mid_by_resnum = {res.id[1]: res for res in mid_res_list}
    layer1_by_resnum = {res.id[1]: res for res in layer1_res_list}

    common_resnums = sorted(set(mid_by_resnum.keys()) & set(layer1_by_resnum.keys()))
    distances = []
    for rnum in common_resnums:
        res0 = mid_by_resnum[rnum]
        res1 = layer1_by_resnum[rnum]
        ca0 = res_ca_coords.get(res0)
        ca1 = res_ca_coords.get(res1)
        if ca0 is None or ca1 is None:
            continue
        dz = abs(ca1[2] - ca0[2])
        distances.append(dz)

    if not distances:
        print(
            "  PF %d: no overlapping residues with CA atoms between "
            "layers 0 and 1; rise will be NA." % pf_index
        )
        continue

    rise_val = sum(distances) / float(len(distances))
    pf_rise[pf_index] = rise_val
    print(
        "  PF %d: rise = %.3f A (from %d CA-CA pairs)"
        % (pf_index, rise_val, len(distances))
    )
print()

# 8a) Sequence index per residue (per chain, PDB order)

res_chain_seq_index = {}

for cid, res_list in chain_residues.items():
    prev_res = None
    prev_seq_index = None
    for res in res_list:
        resnum = res.id[1]
        if prev_res is None:
            seq_index = 0
        else:
            resnum_prev = prev_res.id[1]
            delta_resnum = abs(resnum - resnum_prev)
            if delta_resnum < 1:
                delta_resnum = 1
            jump = delta_resnum

            if delta_resnum != 1:
                ca_prev = get_ca_atom(prev_res)
                ca_curr = get_ca_atom(res)
                if ca_prev is not None and ca_curr is not None:
                    d_ca = ca_distance(ca_prev, ca_curr)
                    if d_ca < max_backbone_gap_distance:
                        # Treat as consecutive despite numbering jump
                        jump = 1

            seq_index = prev_seq_index + jump

        res_chain_seq_index[res] = seq_index
        prev_res = res
        prev_seq_index = seq_index

# 8b) Vector utilities

def vec_sub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def vec_dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def vec_cross(a, b):
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def vec_norm(a):
    return math.sqrt(vec_dot(a, a))

# 8c) Axial index per residue (per protofilament)

chain_axial_offset = {}
res_axial_index = {}


def median_int(values):
    vals = sorted(values)
    n = len(vals)
    if n == 0:
        return 0
    mid = n // 2
    if n % 2 == 1:
        return int(round(vals[mid]))
    else:
        return int(round(0.5 * (vals[mid - 1] + vals[mid])))


for pf_index, pf in enumerate(protofilaments):
    if not pf:
        continue

    mid_cid = pf_to_middle_chain_id.get(pf_index, None)
    if mid_cid is None:
        mid_cid = pf[len(pf) // 2]

    mid_res_list = chain_residues[mid_cid]
    mid_resnum_to_seq = {
        res.id[1]: res_chain_seq_index[res]
        for res in mid_res_list
        if res in res_chain_seq_index
    }

    chain_axial_offset[(pf_index, mid_cid)] = 0

    for cid in pf:
        if cid == mid_cid:
            continue
        res_list = chain_residues[cid]
        resnum_to_seq = {
            res.id[1]: res_chain_seq_index[res]
            for res in res_list
            if res in res_chain_seq_index
        }
        common_resnums = set(mid_resnum_to_seq.keys()) & set(resnum_to_seq.keys())
        if common_resnums:
            offsets = []
            for rnum in common_resnums:
                i_mid = mid_resnum_to_seq[rnum]
                i_chain = resnum_to_seq[rnum]
                offsets.append(i_mid - i_chain)
            offset_c = median_int(offsets)
        else:
            len_mid = len(mid_res_list)
            len_c = len(res_list)
            offset_c = int(round(0.5 * (len_mid - len_c)))
        chain_axial_offset[(pf_index, cid)] = offset_c

for pf_index, pf in enumerate(protofilaments):
    for cid in pf:
        offset = chain_axial_offset.get((pf_index, cid), 0)
        for res in chain_residues[cid]:
            seq_idx = res_chain_seq_index.get(res, None)
            if seq_idx is not None:
                res_axial_index[res] = seq_idx + offset

# 8d) Outer-side planes

def choose_neighbor_chain_for_plane(pf_chain_ids, cid):
    off_self = chain_offset[cid]
    best_c = None
    best_delta = None
    for c in pf_chain_ids:
        if c == cid:
            continue
        off_other = chain_offset[c]
        delta = abs(off_other - off_self)
        if best_delta is None or delta < best_delta:
            best_delta = delta
            best_c = c
    return best_c


def compute_outer_side_planes_for_pf(pf_chain_ids):
    for cid in pf_chain_ids:
        neighbor_cid = choose_neighbor_chain_for_plane(pf_chain_ids, cid)
        if neighbor_cid is None:
            continue

        res_list = chain_residues[cid]
        neighbor_res_list = chain_residues[neighbor_cid]
        neighbor_by_resnum = {res.id[1]: res for res in neighbor_res_list}

        for res_i in res_list:
            resnum_i = res_i.id[1]
            res_neig_i = neighbor_by_resnum.get(resnum_i, None)
            if res_neig_i is None:
                continue

            try:
                N_mid = res_i["N"]
                C_mid = res_i["C"]
                N_neig = res_neig_i["N"]
            except KeyError:
                continue

            P_mid_N = (
                float(N_mid.coord[0]),
                float(N_mid.coord[1]),
                float(N_mid.coord[2]),
            )
            P_neig_N = (
                float(N_neig.coord[0]),
                float(N_neig.coord[1]),
                float(N_neig.coord[2]),
            )
            P_mid_C = (
                float(C_mid.coord[0]),
                float(C_mid.coord[1]),
                float(C_mid.coord[2]),
            )

            v1 = vec_sub(P_neig_N, P_mid_N)
            v2 = vec_sub(P_mid_C, P_mid_N)
            n_vec = vec_cross(v1, v2)
            n_norm = vec_norm(n_vec)
            if n_norm < 1.0e-8:
                continue
            n_unit = (n_vec[0] / n_norm, n_vec[1] / n_norm, n_vec[2] / n_norm)
            res_plane_data[res_i] = {"normal": n_unit, "point": P_mid_N}

for pf in protofilaments:
    if not pf:
        continue
    compute_outer_side_planes_for_pf(pf)

# 8e) Handedness

def compute_handedness_for_middle_chain(pf_chain_ids):
    eps = 1.0e-6

    mid_cid = None
    for cid in pf_chain_ids:
        if chain_offset.get(cid, None) == 0:
            mid_cid = cid
            break
    if mid_cid is None:
        mid_cid = pf_chain_ids[len(pf_chain_ids) // 2]

    neighbor_cid = None
    for cid in pf_chain_ids:
        if chain_offset.get(cid, None) == 1:
            neighbor_cid = cid
            break
    if neighbor_cid is None:
        for cid in pf_chain_ids:
            if chain_offset.get(cid, None) == -1:
                neighbor_cid = cid
                break
    if neighbor_cid is None:
        return

    mid_res_list = chain_residues[mid_cid]
    neighbor_res_list = chain_residues[neighbor_cid]
    neighbor_by_resnum = {res.id[1]: res for res in neighbor_res_list}

    n_mid = len(mid_res_list)

    for idx, res_i in enumerate(mid_res_list):
        res_handedness[res_i] = "NA"

        resname = res_i.get_resname().upper()

        if idx == 0 or idx == n_mid - 1:
            res_handedness[res_i] = "NA"
            continue
        if resname == "GLY":
            res_handedness[res_i] = "NA"
            continue

        try:
            ca_i = res_i["CA"]
            cb_i = res_i["CB"]
        except KeyError:
            res_handedness[res_i] = "NA"
            continue

        res_prev_mid = mid_res_list[idx - 1]
        res_next_mid = mid_res_list[idx + 1]

        ca_prev_mid = get_ca_atom(res_prev_mid)
        ca_next_mid = get_ca_atom(res_next_mid)
        if ca_prev_mid is None or ca_next_mid is None:
            res_handedness[res_i] = "NA"
            continue

        resnum_prev = res_prev_mid.id[1]
        res_prev_neig = neighbor_by_resnum.get(resnum_prev, None)
        if res_prev_neig is None:
            res_handedness[res_i] = "NA"
            continue
        ca_prev_neig = get_ca_atom(res_prev_neig)
        if ca_prev_neig is None:
            res_handedness[res_i] = "NA"
            continue

        P_prev_mid = (
            float(ca_prev_mid.coord[0]),
            float(ca_prev_mid.coord[1]),
            float(ca_prev_mid.coord[2]),
        )
        P_prev_neig = (
            float(ca_prev_neig.coord[0]),
            float(ca_prev_neig.coord[1]),
            float(ca_prev_neig.coord[2]),
        )
        P_next_mid = (
            float(ca_next_mid.coord[0]),
            float(ca_next_mid.coord[1]),
            float(ca_next_mid.coord[2]),
        )

        z_vec = vec_sub(P_next_mid, P_prev_mid)
        z_norm = vec_norm(z_vec)
        if z_norm < 1.0e-8:
            res_handedness[res_i] = "NA"
            continue
        z = (z_vec[0] / z_norm, z_vec[1] / z_norm, z_vec[2] / z_norm)

        u_raw = vec_sub(P_prev_neig, P_prev_mid)
        proj = vec_dot(u_raw, z)
        u_vec = (
            u_raw[0] - proj * z[0],
            u_raw[1] - proj * z[1],
            u_raw[2] - proj * z[2],
        )
        u_norm = vec_norm(u_vec)
        if u_norm < 1.0e-8:
            res_handedness[res_i] = "NA"
            continue
        y = (u_vec[0] / u_norm, u_vec[1] / u_norm, u_vec[2] / u_norm)

        x_vec = vec_cross(z, y)
        x_norm = vec_norm(x_vec)
        if x_norm < 1.0e-8:
            res_handedness[res_i] = "NA"
            continue
        x = (x_vec[0] / x_norm, x_vec[1] / x_norm, x_vec[2] / x_norm)

        ca_i_coord = (
            float(ca_i.coord[0]),
            float(ca_i.coord[1]),
            float(ca_i.coord[2]),
        )
        cb_i_coord = (
            float(cb_i.coord[0]),
            float(cb_i.coord[1]),
            float(cb_i.coord[2]),
        )
        s_vec = vec_sub(cb_i_coord, ca_i_coord)
        s_norm = vec_norm(s_vec)
        if s_norm < 1.0e-8:
            res_handedness[res_i] = "NA"
            continue
        s = (s_vec[0] / s_norm, s_vec[1] / s_norm, s_vec[2] / s_norm)

        proj_sx = vec_dot(s, x)
        if proj_sx > eps:
            res_handedness[res_i] = "R"
        elif proj_sx < -eps:
            res_handedness[res_i] = "L"
        else:
            res_handedness[res_i] = "NA"

for pf in protofilaments:
    if not pf:
        continue
    compute_handedness_for_middle_chain(pf)

# 8f) Outer-side segment filter

def sidechain_vector(res):
    rn = res.get_resname().upper()
    if rn == "GLY":
        return None
    if "CA" not in res or "CB" not in res:
        return None
    ca = res["CA"]
    cb = res["CB"]
    vx = float(cb.coord[0] - ca.coord[0])
    vy = float(cb.coord[1] - ca.coord[1])
    vz = float(cb.coord[2] - ca.coord[2])
    n2 = vx * vx + vy * vy + vz * vz
    if n2 < 1.0e-8:
        return None
    n = math.sqrt(n2)
    return (vx / n, vy / n, vz / n)


def segment_on_outer_side_for_residue(res, p1, p2):
    pdata = res_plane_data.get(res, None)
    svec = sidechain_vector(res)
    if pdata is None or svec is None:
        return True

    nvec = pdata["normal"]
    plane_point = pdata["point"]

    try:
        cb = res["CB"]
    except KeyError:
        return True

    cb_coord = (
        float(cb.coord[0]),
        float(cb.coord[1]),
        float(cb.coord[2]),
    )

    vs = (
        cb_coord[0] - plane_point[0],
        cb_coord[1] - plane_point[1],
        cb_coord[2] - plane_point[2],
    )
    s_sc = nvec[0] * vs[0] + nvec[1] * vs[1] + nvec[2] * vs[2]
    if abs(s_sc) < outer_side_eps:
        return True

    v1 = (
        p1[0] - plane_point[0],
        p1[1] - plane_point[1],
        p1[2] - plane_point[2],
    )
    v2 = (
        p2[0] - plane_point[0],
        p2[1] - plane_point[1],
        p2[2] - plane_point[2],
    )
    s1 = nvec[0] * v1[0] + nvec[1] * v1[1] + nvec[2] * v1[2]
    s2 = nvec[0] * v2[0] + nvec[1] * v2[1] + nvec[2] * v2[2]

    if abs(s1) < outer_side_eps or abs(s2) < outer_side_eps:
        other = s2 if abs(s1) < outer_side_eps else s1
        if abs(other) < outer_side_eps:
            return True
        return (s_sc * other) > 0.0

    if (s_sc * s1) <= 0.0:
        return False
    if (s_sc * s2) <= 0.0:
        return False

    return True

# 8g) Real-valued ca_offset (fractional number of rises, middle-layer reference)

def compute_ca_offset_for_pair(pf_src, resnum_src, pf_partner, resnum_partner):
    """
    Compute ca_offset for a contact between residue1 (pf_src, resnum_src)
    and residue2 (pf_partner, resnum_partner) as a fractional number of rises,
    based on the vertical separation between their middle-layer CA atoms.

    ca_offset = (z_mid_partner(residue2) - z_mid_src(residue1)) / rise

    The rise used is the average of the available per-protofilament rises
    for pf_src and pf_partner (if only one is defined and positive, that
    value is used).
    """
    rise_src = pf_rise.get(pf_src, None)
    rise_partner = pf_rise.get(pf_partner, None)

    rise_candidates = []
    if rise_src is not None and rise_src > 0.0:
        rise_candidates.append(rise_src)
    if rise_partner is not None and rise_partner > 0.0:
        rise_candidates.append(rise_partner)

    if not rise_candidates:
        return None

    rise_val = sum(rise_candidates) / float(len(rise_candidates))

    mid_src_cid = pf_to_middle_chain_id.get(pf_src, None)
    mid_partner_cid = pf_to_middle_chain_id.get(pf_partner, None)
    if mid_src_cid is None or mid_partner_cid is None:
        return None

    res_mid_src = None
    for res in chain_residues.get(mid_src_cid, []):
        if res.id[1] == resnum_src:
            res_mid_src = res
            break
    if res_mid_src is None:
        return None

    res_mid_partner = None
    for res in chain_residues.get(mid_partner_cid, []):
        if res.id[1] == resnum_partner:
            res_mid_partner = res
            break
    if res_mid_partner is None:
        return None

    ca_src = res_ca_coords.get(res_mid_src)
    ca_partner = res_ca_coords.get(res_mid_partner)
    if ca_src is None or ca_partner is None:
        return None

    z_src = ca_src[2]
    z_partner = ca_partner[2]

    delta_z = z_partner - z_src
    offset_float = delta_z / rise_val

    return offset_float

# 9) Collect contacts

contacts = {}

ns_sc = NeighborSearch(sidechain_atoms_all)

middle_chain_set = set(middle_chains)

middle_residues_all = []
for chain in middle_chains:
    for res in chain_residues[chain.id]:
        middle_residues_all.append(res)

middle_res_index = {res: idx for idx, res in enumerate(middle_residues_all)}

chain_res_by_num = {}
for cid, res_list in chain_residues.items():
    for res in res_list:
        chain_res_by_num[(cid, res.id[1])] = res

for chain_obj in middle_chains:
    cid_src = chain_obj.id
    pf_src = chain_to_pf[cid_src]

    for res in chain_residues[cid_src]:
        sc_list_i = res_sidechain_atoms.get(res, [])
        if not sc_list_i:
            continue

        axial_i = res_axial_index.get(res, None)
        resnum_i = res.id[1]

        for a in sc_list_i:
            for nb in ns_sc.search(a.coord, sidechain_cutoff):
                other_res = nb.get_parent()
                if other_res is res:
                    continue
                if not is_aa(other_res, standard=True):
                    continue

                sc_list_j = res_sidechain_atoms.get(other_res, [])
                if not sc_list_j:
                    continue

                chain_j = other_res.get_parent()
                cid_j = chain_j.id
                pf_j = chain_to_pf[cid_j]
                axial_j = res_axial_index.get(other_res, None)
                resnum_j = other_res.id[1]

                d = atom_distance(a, nb)
                if d > sidechain_cutoff:
                    continue

                if (
                    resnum_i in pf_mid_resnums.get(pf_src, set())
                    and resnum_j in pf_mid_resnums.get(pf_j, set())
                ):
                    if (pf_j, resnum_j) <= (pf_src, resnum_i):
                        continue

                if pf_src == pf_j:
                    if axial_i is not None and axial_j is not None:
                        if abs(axial_i - axial_j) < min_sequence_separation:
                            continue
                    ctype = "intra"
                else:
                    ctype = "inter"

                p1 = (
                    float(a.coord[0]),
                    float(a.coord[1]),
                    float(a.coord[2]),
                )
                p2 = (
                    float(nb.coord[0]),
                    float(nb.coord[1]),
                    float(nb.coord[2]),
                )

                if not segment_on_outer_side_for_residue(res, p1, p2):
                    continue
                if not segment_on_outer_side_for_residue(other_res, p1, p2):
                    continue

                key = (pf_src, resnum_i, pf_j, resnum_j, ctype)
                partner_chain_offset = chain_offset[cid_j]

                old = contacts.get(key, None)
                if old is None or d < old["distance"]:
                    contacts[key] = {
                        "distance": d,
                        "p1": p1,
                        "p2": p2,
                        "source_chain_id": cid_src,
                        "partner_chain_id": cid_j,
                        "partner_chain_offset": partner_chain_offset,
                    }

print("Total contact pairs (after outer-side segment filter):", len(contacts))

# 10) Write contacts.csv

with open(csv_file, "w") as fout:
    fout.write(
        "pf,"
        "residue,resname,source_chain,handedness,"
        "ca_x,ca_y,ca_z,"
        "cb_x,cb_y,cb_z,"
        "partner_pf,partner_chain_offset,ca_offset,partner_residue,partner_resname,partner_chain,"
        "contact_type,distance\n"
    )

    contacts_by_source = defaultdict(list)
    for (pf_src, resnum_src, pf_partner, resnum_partner, ctype), info in contacts.items():
        contacts_by_source[(pf_src, resnum_src)].append(
            (pf_partner, resnum_partner, ctype, info)
        )

    for pf_index, pf in enumerate(protofilaments):
        if not pf:
            continue

        mid_cid = pf_to_middle_chain_id.get(pf_index, None)
        if mid_cid is None:
            mid_cid = pf[len(pf) // 2]
        mid_chain = chain_by_id[mid_cid]

        for res in chain_residues[mid_cid]:
            resnum_i = res.id[1]
            resname_i = res.get_resname()
            source_chain_id = mid_chain.id
            handed = res_handedness.get(res, "NA")
            ca_src = res_ca_coords.get(res, None)
            cb_src = res_cb_coords.get(res, None)

            if ca_src is None:
                src_ca_x = src_ca_y = src_ca_z = "NA"
            else:
                src_ca_x = "%.3f" % ca_src[0]
                src_ca_y = "%.3f" % ca_src[1]
                src_ca_z = "%.3f" % ca_src[2]

            if cb_src is None:
                src_cb_x = src_cb_y = src_cb_z = "NA"
            else:
                src_cb_x = "%.3f" % cb_src[0]
                src_cb_y = "%.3f" % cb_src[1]
                src_cb_z = "%.3f" % cb_src[2]

            key_src = (pf_index, resnum_i)
            if key_src not in contacts_by_source:
                fout.write(
                    "%d,%d,%s,%s,%s,"
                    "%s,%s,%s,"
                    "%s,%s,%s,"
                    "NA,NA,NA,NA,NA,NA,"
                    "none,NA\n"
                    % (
                        pf_index,
                        resnum_i,
                        resname_i,
                        source_chain_id,
                        handed,
                        src_ca_x,
                        src_ca_y,
                        src_ca_z,
                        src_cb_x,
                        src_cb_y,
                        src_cb_z,
                    )
                )
                continue

            for (pf_partner, resnum_partner, ctype, info) in contacts_by_source[key_src]:
                partner_chain_id = info["partner_chain_id"]
                partner_chain_offset = info["partner_chain_offset"]

                partner_res = chain_res_by_num.get(
                    (partner_chain_id, resnum_partner), None
                )
                if partner_res is not None:
                    partner_resname = partner_res.get_resname()
                else:
                    partner_resname = "NA"

                ca_offset_val = compute_ca_offset_for_pair(
                    pf_index, resnum_i, pf_partner, resnum_partner
                )
                if ca_offset_val is None:
                    ca_offset_str = "NA"
                else:
                    ca_offset_str = "%.3f" % ca_offset_val

                dist_str = "%.3f" % info["distance"]

                fout.write(
                    "%d,%d,%s,%s,%s,"
                    "%s,%s,%s,"
                    "%s,%s,%s,"
                    "%d,%d,%s,%d,%s,%s,"
                    "%s,%s\n"
                    % (
                        pf_index,
                        resnum_i,
                        resname_i,
                        source_chain_id,
                        handed,
                        src_ca_x,
                        src_ca_y,
                        src_ca_z,
                        src_cb_x,
                        src_cb_y,
                        src_cb_z,
                        pf_partner,
                        partner_chain_offset,
                        ca_offset_str,
                        resnum_partner,
                        partner_resname,
                        partner_chain_id,
                        ctype,
                        dist_str,
                    )
                )

print("CSV file written to:", csv_file)

# 11) Write BILD file with all contacts

with open(bild_file, "w") as f:
    for (pf_src, resnum_src, pf_partner, resnum_partner, ctype), info in contacts.items():
        p1 = info["p1"]
        p2 = info["p2"]

        if ctype == "intra":
            f.write(".color 0 1 0\n")
        else:
            f.write(".color 1 1 0\n")

        f.write(
            ".cylinder %.3f %.3f %.3f %.3f %.3f %.3f %.3f\n"
            % (
                p1[0],
                p1[1],
                p1[2],
                p2[0],
                p2[1],
                p2[2],
                cylinder_radius,
            )
        )

print("BILD file written to:", bild_file)
print("Done.")

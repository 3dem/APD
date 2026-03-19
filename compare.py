#!/usr/bin/env python

import sys
import os
import csv
import math
from collections import defaultdict

# compare.py
#
# Compare amyloid protofilament contacts between two structures using
# the CSV output of contacts.py.
#
# CSV format (from contacts.py):
#   pf,
#   residue,resname,source_chain,handedness,
#   ca_x,ca_y,ca_z,
#   cb_x,cb_y,cb_z,
#   partner_pf,partner_chain_offset,ca_offset,partner_residue,partner_resname,partner_chain,
#   contact_type,distance
#
# Key behavior:
#   - Contacts are classified as common / unique based on strict vs relaxed cutoffs.
#   - right-left handedness can be globally flipped for structure 2 if it reduces mismatches.
#   - Z-offset (ca_offset) comparisons:
#       * Intra-protofilament contacts: use rounded delta = round(co1 - co2),
#         with optional global sign flip for structure 2 (up-side-down check).
#       * Inter-protofilament contacts: in addition to the optional sign flip,
#         a global integer shift in {-2,-1,0,1,2} can be applied to structure 2
#         to minimize the number of inter contacts with non-zero rounded delta.
#   - SVG coloring of common contacts uses abs(rounded delta) buckets:
#       0, 1, 2, >=3
#     using the intra rule for intra contacts and the inter rule (with shift) for inter contacts.
#   - Summary:
#       * Protofilament blocks: Z-offsets are counted only from common intra contacts.
#       * Interface blocks (pf i-j): Z-offsets are counted from common inter contacts,
#         using the best inter shift; XY and XYZ distances are reported.
#

# Distance thresholds (Angstrom)
strict_cutoff_nonterminal = 4.5
strict_cutoff_terminal = 3.5
relaxed_cutoff = 6.5

# Backbone CA-CA maximum distance to draw black lines
max_backbone_ca_distance = 4.2

# SVG parameters
svg_scale = 10.0
svg_margin = 5.0
svg_panel_spacing = 10.0
svg_circle_radius = 1.2  # in Angstrom, then scaled

# Allow global flip of handedness in structure 2
allow_handedness_global_flip = True

# Allow global sign flip of Z-offsets (ca_offset) in structure 2 (up-side-down check)
allow_zoffset_global_flip = True

# Allow global integer shift of inter-protofilament Z-offsets for structure 2
# Tested shifts are in inter_zoffset_shift_candidates.
allow_inter_zoffset_global_shift = True
inter_zoffset_shift_candidates = [-2, -1, 0, 1, 2]

one_letter_code = {
    "ARG": "R", "HIS": "H", "LYS": "K", "ASP": "D", "GLU": "E",
    "SER": "S", "THR": "T", "ASN": "N", "GLN": "Q", "CYS": "C",
    "GLY": "G", "PRO": "P", "ALA": "A", "VAL": "V", "ILE": "I",
    "LEU": "L", "MET": "M", "PHE": "F", "TYR": "Y", "TRP": "W",
}


def parse_float_or_none(s):
    s = s.strip()
    if s in ("", "NA", "NaN", "nan"):
        return None
    return float(s)


def parse_int_or_none(s):
    s = s.strip()
    if s in ("", "NA"):
        return None
    return int(s)


def normalize_handed_label(h):
    """Normalize handedness labels; treat blank/NA-like values as missing."""
    h = str(h).strip()
    if h.upper() in ("", "NA", "N/A", "NONE", "NULL"):
        return None
    return h


def flip_handed_label(h):
    """Flip handedness label for global comparison."""
    h = normalize_handed_label(h)
    if h == "R":
        return "L"
    if h == "L":
        return "R"
    if h == "1":
        return "-1"
    if h == "-1":
        return "1"
    return h


def load_contacts_csv(csv_path, residue_offset=0):
    """Load residues and contacts from a contacts.py CSV file.

    residue_offset:
      Integer offset applied to ALL residue numbers in this file:
        new_resnum      = resnum_raw + residue_offset
        new_partner_res = partner_res_raw + residue_offset

    Returns:
      residues: dict[(pf, residue)] -> {
          "resname": str,
          "handed": str,
          "ca": (x, y, z) or None,
          "cb": (x, y, z) or None,
          "chain_id": str,
      }

      contacts: dict[(pf, res, partner_pf, partner_res, ctype)] -> {
          "distance": float,
          "terminal": bool,
          "ca_offset": float or None,
      }
    """
    residues = {}
    contacts = {}
    chain_resnums = defaultdict(set)
    raw_contacts = []

    with open(csv_path, "r", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            pf = parse_int_or_none(row["pf"])
            if pf is None:
                continue

            resnum_raw = parse_int_or_none(row["residue"])
            if resnum_raw is None:
                continue
            resnum = resnum_raw + residue_offset

            resname = row["resname"].strip()
            handed = row["handedness"].strip()
            source_chain = row["source_chain"].strip()

            ca_x = parse_float_or_none(row["ca_x"])
            ca_y = parse_float_or_none(row["ca_y"])
            ca_z = parse_float_or_none(row["ca_z"])
            if ca_x is not None and ca_y is not None and ca_z is not None:
                ca_coord = (ca_x, ca_y, ca_z)
            else:
                ca_coord = None

            cb_x = parse_float_or_none(row["cb_x"])
            cb_y = parse_float_or_none(row["cb_y"])
            cb_z = parse_float_or_none(row["cb_z"])
            if cb_x is not None and cb_y is not None and cb_z is not None:
                cb_coord = (cb_x, cb_y, cb_z)
            else:
                cb_coord = None

            rkey = (pf, resnum)
            if rkey not in residues:
                residues[rkey] = {
                    "resname": resname,
                    "handed": handed,
                    "ca": ca_coord,
                    "cb": cb_coord,
                    "chain_id": source_chain,
                }
            else:
                r = residues[rkey]
                if r["ca"] is None and ca_coord is not None:
                    r["ca"] = ca_coord
                if r["cb"] is None and cb_coord is not None:
                    r["cb"] = cb_coord
                if not r["resname"]:
                    r["resname"] = resname
                if not r["handed"]:
                    r["handed"] = handed
                if not r["chain_id"]:
                    r["chain_id"] = source_chain

            chain_resnums[source_chain].add(resnum)

            ctype = row["contact_type"].strip()
            partner_pf = parse_int_or_none(row["partner_pf"])
            partner_res_raw = parse_int_or_none(row["partner_residue"])
            partner_chain = row["partner_chain"].strip()

            partner_offset_str = row.get("partner_chain_offset", "").strip()
            if partner_offset_str in ("", "NA"):
                partner_offset = None
            else:
                partner_offset = int(partner_offset_str)

            ca_offset = parse_float_or_none(row.get("ca_offset", "").strip())

            if partner_pf is not None and partner_res_raw is not None and partner_chain:
                partner_res = partner_res_raw + residue_offset
                chain_resnums[partner_chain].add(partner_res)
            else:
                partner_pf = None
                partner_res = None

            if ctype == "none":
                continue

            dist = parse_float_or_none(row["distance"])
            if dist is None:
                continue

            if partner_pf is None or partner_res is None:
                continue

            raw_contacts.append({
                "pf": pf,
                "resnum": resnum,
                "chain_id": source_chain,
                "partner_pf": partner_pf,
                "partner_res": partner_res,
                "partner_chain": partner_chain,
                "partner_chain_offset": partner_offset,
                "ctype": ctype,
                "distance": dist,
                "ca_offset": ca_offset,
            })

    terminal_positions = set()
    for cid, resnums in chain_resnums.items():
        if not resnums:
            continue
        rmin = min(resnums)
        rmax = max(resnums)
        terminal_positions.add((cid, rmin))
        terminal_positions.add((cid, rmax))

    for rc in raw_contacts:
        key = (rc["pf"], rc["resnum"], rc["partner_pf"], rc["partner_res"], rc["ctype"])

        is_terminal = (
            (rc["chain_id"], rc["resnum"]) in terminal_positions
            or (rc["partner_chain"], rc["partner_res"]) in terminal_positions
        )

        info = contacts.get(key)
        if info is None:
            contacts[key] = {
                "distance": rc["distance"],
                "terminal": is_terminal,
                "ca_offset": rc["ca_offset"],
            }
        else:
            if rc["distance"] < info["distance"]:
                info["distance"] = rc["distance"]
                info["terminal"] = is_terminal
                info["ca_offset"] = rc["ca_offset"]
            elif rc["distance"] == info["distance"]:
                if info["ca_offset"] is None and rc["ca_offset"] is not None:
                    info["ca_offset"] = rc["ca_offset"]
                if not info["terminal"] and is_terminal:
                    info["terminal"] = True

    return residues, contacts


def classify_contacts(contacts1, contacts2):
    """Classify contacts into common and unique sets."""
    keys1 = set(contacts1.keys())
    keys2 = set(contacts2.keys())
    all_keys = keys1 | keys2

    common_keys = set()
    unique1_keys = set()
    unique2_keys = set()

    # IMPORTANT (inter-/intra-PF consistency):
    # The project definition for contact classification uses a single strict
    # threshold (4.5 Å) and a relaxed threshold (6.5 Å):
    #   * Common: strict in one structure AND relaxed in the other.
    #   * Unique: strict in one structure AND absent even at relaxed in the other.
    # This must apply identically for intra- and inter-protofilament contacts.
    #
    # Previous versions used a special "terminal" strict cutoff (3.5 Å) and
    # additionally *excluded* terminal contacts from the unique sets. In the CSVs
    # from contacts.py, "terminal" is only inferred from the min/max residue
    # numbers that appear in the contact table, which can incorrectly label
    # internal residues (often glycines) as terminal when only a contact subset
    # is present. That caused missing inter-PF contacts involving glycines.
    #
    # We therefore ignore the inferred terminal flag for classification.
    def get_info(info):
        if info is None:
            return None, False, False
        d = info["distance"]
        s = d <= strict_cutoff_nonterminal
        r = d <= relaxed_cutoff
        return d, s, r

    for key in all_keys:
        info1 = contacts1.get(key)
        info2 = contacts2.get(key)

        _, s1, r1 = get_info(info1)
        _, s2, r2 = get_info(info2)

        if info1 is not None and info2 is not None:
            if (s1 and r2) or (s2 and r1):
                common_keys.add(key)
                continue

        if info1 is not None and s1 and not (info2 is not None and r2):
            unique1_keys.add(key)

        if info2 is not None and s2 and not (info1 is not None and r1):
            unique2_keys.add(key)

    return common_keys, unique1_keys, unique2_keys


def project_point(point3d, flip_x, rot_rad):
    """Project 3D point (x, y, z) to 2D (x, y) with optional flip and rotation."""
    x, y, z = point3d
    if flip_x:
        x = -x
    cos_a = math.cos(rot_rad)
    sin_a = math.sin(rot_rad)
    x2 = cos_a * x - sin_a * y
    y2 = sin_a * x + cos_a * y
    return (x2, y2)


def compute_projected_positions(residues, residue_keys, flip_x, rot_rad):
    """Compute projected CA and CB positions for a given set of residues (keys are (pf,res))."""
    ca_proj = {}
    cb_proj = {}

    minx = float("inf")
    maxx = float("-inf")
    miny = float("inf")
    maxy = float("-inf")

    for key in residue_keys:
        r = residues.get(key)
        if r is None:
            continue
        ca = r.get("ca")
        cb = r.get("cb")
        if ca is None:
            continue

        x_ca, y_ca = project_point(ca, flip_x, rot_rad)
        ca_proj[key] = (x_ca, y_ca)
        minx = min(minx, x_ca)
        maxx = max(maxx, x_ca)
        miny = min(miny, y_ca)
        maxy = max(maxy, y_ca)

        if cb is not None:
            x_cb, y_cb = project_point(cb, flip_x, rot_rad)
        else:
            x_cb, y_cb = x_ca, y_ca

        cb_proj[key] = (x_cb, y_cb)
        minx = min(minx, x_cb)
        maxx = max(maxx, x_cb)
        miny = min(miny, y_cb)
        maxy = max(maxy, y_cb)

    if minx == float("inf"):
        return {}, {}, (0.0, 0.0, 0.0, 0.0)

    return ca_proj, cb_proj, (minx, maxx, miny, maxy)


def get_ca_offset(info):
    if info is None:
        return None
    return info.get("ca_offset", None)


def zoffset_bucket(co1, co2):
    if co1 is None or co2 is None:
        return 0
    return int(round(co1 - co2))


def compute_zoffset_bucket_intra(info1, info2, flip_sign_struct2):
    co1 = get_ca_offset(info1)
    co2 = get_ca_offset(info2)
    if flip_sign_struct2 and co2 is not None:
        co2 = -co2
    return zoffset_bucket(co1, co2)


def compute_zoffset_bucket_inter(info1, info2, flip_sign_struct2, inter_shift):
    co1 = get_ca_offset(info1)
    co2 = get_ca_offset(info2)
    if flip_sign_struct2 and co2 is not None:
        co2 = -co2
    if co2 is not None:
        co2 = co2 + float(inter_shift)
    return zoffset_bucket(co1, co2)


def choose_best_inter_shift(common_keys, contacts1, contacts2, flip_sign_struct2, candidates):
    """Choose inter_shift minimizing the number of common inter contacts with non-zero rounded delta.

    Tie-break:
      1) prefer shift == 0
      2) then prefer smallest abs(shift)
      3) then prefer negative (stable deterministic)
    """
    counts = []
    for sh in candidates:
        nonzero = 0
        total = 0
        for key in common_keys:
            pf, res, partner_pf, partner_res, ctype = key
            if ctype != "inter":
                continue
            total += 1
            b = compute_zoffset_bucket_inter(contacts1.get(key), contacts2.get(key), flip_sign_struct2, sh)
            if abs(b) >= 1:
                nonzero += 1
        counts.append((nonzero, total, sh))

    # If no inter contacts, keep shift 0
    any_inter = any(tot > 0 for _, tot, _ in counts)
    if not any_inter:
        return 0, 0, 0

    # Rank by objective then tie-break
    def rank_key(tup):
        nonzero, total, sh = tup
        return (
            nonzero,
            0 if sh == 0 else 1,
            abs(sh),
            0 if sh < 0 else 1,
        )

    counts_sorted = sorted(counts, key=rank_key)
    best = counts_sorted[0]
    return best[2], best[0], best[1]


def write_svg(
    svg_path,
    residues1,
    residues2,
    contacts1,
    contacts2,
    common_keys,
    unique1_keys,
    unique2_keys,
    mutated_residues,
    handed_diff_residues,
    common_residues,
    flip1,
    rot1_deg,
    flip2,
    rot2_deg,
    flip_zoffset_sign_struct2,
    inter_zoffset_shift,
    pf_stats_lines,
):
    """Write SVG with two panels (structure 1 and 2), legend, and stats."""
    rot1_rad = math.radians(rot1_deg)
    rot2_rad = math.radians(rot2_deg)

    base_name = os.path.basename(svg_path)
    name_no_ext = os.path.splitext(base_name)[0]
    parts = name_no_ext.split("_vs_")
    if len(parts) >= 2:
        label1 = parts[0]
        label2 = parts[1]
    else:
        label1 = name_no_ext
        label2 = name_no_ext
    if label1.endswith("_contacts"):
        label1 = label1[:-9]
    if label2.endswith("_contacts"):
        label2 = label2[:-9]

    keys1_all = set(residues1.keys())
    keys2_all = set(residues2.keys())
    unique_residues1 = keys1_all - set(common_residues)
    unique_residues2 = keys2_all - set(common_residues)

    draw_residues1 = set(common_residues) | unique_residues1
    draw_residues2 = set(common_residues) | unique_residues2

    ca1_proj, cb1_proj, bounds1 = compute_projected_positions(residues1, draw_residues1, flip1, rot1_rad)
    ca2_proj, cb2_proj, bounds2 = compute_projected_positions(residues2, draw_residues2, flip2, rot2_rad)

    if not ca1_proj and not ca2_proj:
        print("No drawable residues for SVG; skipping SVG output.")
        return

    minx1, maxx1, miny1, maxy1 = bounds1
    minx2, maxx2, miny2, maxy2 = bounds2

    width1 = 0.0 if minx1 == maxx1 else (maxx1 - minx1)
    height1 = 0.0 if miny1 == maxy1 else (maxy1 - miny1)
    height2 = 0.0 if miny2 == maxy2 else (maxy2 - miny2)

    offset1_x = svg_margin - minx1
    offset2_x = svg_margin + width1 + svg_panel_spacing + svg_margin - minx2

    half1 = height1 / 2.0
    half2 = height2 / 2.0
    max_half = max(half1, half2)
    global_center_y = svg_margin + max_half

    center1 = miny1 + half1 if height1 > 0 else 0.0
    center2 = miny2 + half2 if height2 > 0 else 0.0

    offset1_y = global_center_y - center1
    offset2_y = global_center_y - center2

    circle_r_world = svg_circle_radius
    structures_bottom_world = global_center_y + max_half + circle_r_world

    circle_r_px = svg_circle_radius * svg_scale
    row_step_px = circle_r_px * 2.5
    row_step_world = row_step_px / svg_scale

    label_y_world = structures_bottom_world + svg_margin * 0.5
    legend_y_world = label_y_world + svg_margin * 1.0

    # Estimate space for legend + summary
    legend_start_world = legend_y_world
    y_world = legend_start_world + 3.0 * row_step_world
    y_world += 1.2 * row_step_world
    y_world += 10.0 * 0.8 * row_step_world
    legend_bottom_world = y_world

    stats_header_y_world = legend_bottom_world + 1.2 * row_step_world
    n_pf_lines = len(pf_stats_lines)
    stats_bottom_world = (
        stats_header_y_world
        + 0.7 * row_step_world
        + max(0, n_pf_lines - 1) * 0.6 * row_step_world
        + 0.5 * row_step_world
    )

    bottom_world = stats_bottom_world + svg_margin
    rightmost = offset2_x + maxx2
    total_width_world = rightmost + svg_margin

    total_width = total_width_world * svg_scale
    total_height = bottom_world * svg_scale

    def world_to_svg(x, y, off_x, off_y):
        xs = (x + off_x) * svg_scale
        ys = (y + off_y) * svg_scale
        return xs, ys

    ca1_svg = {}
    circ1_svg = {}
    for key, (x, y) in ca1_proj.items():
        ca1_svg[key] = world_to_svg(x, y, offset1_x, offset1_y)
    for key, (x, y) in cb1_proj.items():
        circ1_svg[key] = world_to_svg(x, y, offset1_x, offset1_y)

    ca2_svg = {}
    circ2_svg = {}
    for key, (x, y) in ca2_proj.items():
        ca2_svg[key] = world_to_svg(x, y, offset2_x, offset2_y)
    for key, (x, y) in cb2_proj.items():
        circ2_svg[key] = world_to_svg(x, y, offset2_x, offset2_y)

    pf_indices = sorted({pf for (pf, res) in (draw_residues1 | draw_residues2)})

    def backbone_ca_distance(residues_dict, key1, key2):
        r1 = residues_dict.get(key1)
        r2 = residues_dict.get(key2)
        if r1 is None or r2 is None:
            return None
        ca1 = r1.get("ca")
        ca2 = r2.get("ca")
        if ca1 is None or ca2 is None:
            return None
        dx = ca1[0] - ca2[0]
        dy = ca1[1] - ca2[1]
        dz = ca1[2] - ca2[2]
        return math.sqrt(dx * dx + dy * dy + dz * dz)

    with open(svg_path, "w", encoding="ascii") as f:
        f.write(
            '<?xml version="1.0" standalone="no"?>\n'
            '<svg viewBox="0 0 %.1f %.1f" version="1.1" '
            'xmlns="http://www.w3.org/2000/svg">\n' % (total_width, total_height)
        )

        def draw_backbone(residues, ca_svg, residue_keys, unique_residue_set):
            for pf in pf_indices:
                resnums = sorted([r for (p, r) in residue_keys if p == pf])
                last = None
                for r in resnums:
                    keyr = (pf, r)
                    if keyr not in ca_svg:
                        last = None
                        continue
                    if last is not None and last in ca_svg:
                        d_ca = backbone_ca_distance(residues, last, keyr)
                        if d_ca is not None and d_ca <= max_backbone_ca_distance:
                            x1, y1 = ca_svg[last]
                            x2, y2 = ca_svg[keyr]
                            f.write(
                                '  <line x1="%.1f" y1="%.1f" '
                                'x2="%.1f" y2="%.1f" '
                                'stroke="%s" stroke-width="6" stroke-linecap="butt" />\n'
                                % (x1, y1, x2, y2, ("#D3D3D3" if (last in unique_residue_set or keyr in unique_residue_set) else "black"))
                            )
                    last = keyr

        def draw_unique(unique_keys, circ_svg):
            """Draw unique contacts that involve only common residues (legacy behavior)."""
            for key in sorted(unique_keys):
                pf, res, partner_pf, partner_res, ctype = key
                k1 = (pf, res)
                k2 = (partner_pf, partner_res)

                # Contacts involving residues that exist only in one structure are drawn separately
                if k1 not in common_residues or k2 not in common_residues:
                    continue

                if k1 not in circ_svg or k2 not in circ_svg:
                    continue
                x1, y1 = circ_svg[k1]
                x2, y2 = circ_svg[k2]
                if ctype == "intra":
                    color = "#FF5342"
                else:
                    color = "#FF9463"
                f.write(
                    '  <line x1="%.1f" y1="%.1f" '
                    'x2="%.1f" y2="%.1f" '
                    'stroke="%s" stroke-width="6" stroke-linecap="butt" />\n'
                    % (x1, y1, x2, y2, color)
                )

        def draw_unique_residue_contacts(residues_unique_set, contacts_dict, circ_svg):
            """Draw strict contacts that involve at least one residue unique to this structure."""
            for key, info in contacts_dict.items():
                pf, res, partner_pf, partner_res, ctype = key
                k1 = (pf, res)
                k2 = (partner_pf, partner_res)

                if k1 not in residues_unique_set and k2 not in residues_unique_set:
                    continue

                # Only draw strict contacts. Use the same strict definition as
                # classify_contacts() (4.5 Å), ignoring the inferred terminal flag.
                if info is None:
                    continue
                d = info.get("distance", None)
                if d is None:
                    continue
                if d > strict_cutoff_nonterminal:
                    continue

                if k1 not in circ_svg or k2 not in circ_svg:
                    continue

                x1, y1 = circ_svg[k1]
                x2, y2 = circ_svg[k2]
                f.write(
                    '  <line x1="%.1f" y1="%.1f" '
                    'x2="%.1f" y2="%.1f" '
                    'stroke="#F2C14E" stroke-width="3" stroke-linecap="butt" />\n'
                    % (x1, y1, x2, y2)
                )

        def common_bucket_for_key(key):
            info1 = contacts1.get(key)
            info2 = contacts2.get(key)
            ctype = key[4]
            if ctype == "inter":
                return compute_zoffset_bucket_inter(info1, info2, flip_zoffset_sign_struct2, inter_zoffset_shift)
            return compute_zoffset_bucket_intra(info1, info2, flip_zoffset_sign_struct2)

        def draw_common(circ_svg):
            for key in sorted(common_keys):
                pf, res, partner_pf, partner_res, ctype = key
                k1 = (pf, res)
                k2 = (partner_pf, partner_res)
                if k1 not in circ_svg or k2 not in circ_svg:
                    continue

                bucket = common_bucket_for_key(key)
                ab = abs(bucket)
                stroke_width = 6
                if ab == 0:
                    color = "#D3D3D3"
                    stroke_width = 3
                elif ab == 1:
                    color = "#4C7D8E"
                elif ab == 2:
                    color = "#174D7C"
                else:
                    color = "#1C284E"

                x1, y1 = circ_svg[k1]
                x2, y2 = circ_svg[k2]
                f.write(
                    '  <line x1="%.1f" y1="%.1f" '
                    'x2="%.1f" y2="%.1f" '
                    'stroke="%s" stroke-width="%s" stroke-linecap="butt" />\n'
                    % (x1, y1, x2, y2, color, stroke_width)
                )

        def draw_residues(residues, circ_svg, residue_keys, unique_residue_set):
            for keyr in sorted(residue_keys):
                if keyr not in circ_svg:
                    continue
                x, y = circ_svg[keyr]
                rinfo = residues.get(keyr)
                if rinfo is None:
                    continue
                resname = rinfo["resname"].upper()
                letter = one_letter_code.get(resname, "?")

                is_unique = keyr in unique_residue_set
                is_mutated = keyr in mutated_residues
                is_handed_flip = keyr in handed_diff_residues

                if is_unique:
                    # Unique residues: white fill, grey outline, grey label
                    fill_color = "white"
                    text_color = "#D3D3D3"
                else:
                    if is_mutated:
                        fill_color = "#FF5342"
                    elif is_handed_flip:
                        fill_color = "#FF9463"
                    else:
                        fill_color = "white"
                    text_color = "black"

                stroke_color = "#D3D3D3" if is_unique else "black"
                if (not is_unique) and is_mutated and is_handed_flip:
                    f.write(
                        '  <path d="M %.1f %.1f A %.1f %.1f 0 0 1 %.1f %.1f L %.1f %.1f Z" '
                        'fill="#FF5342" stroke="none" />\n'
                        % (x, y - circle_r_px, circle_r_px, circle_r_px, x, y + circle_r_px, x, y)
                    )
                    f.write(
                        '  <path d="M %.1f %.1f A %.1f %.1f 0 0 0 %.1f %.1f L %.1f %.1f Z" '
                        'fill="#FF9463" stroke="none" />\n'
                        % (x, y - circle_r_px, circle_r_px, circle_r_px, x, y + circle_r_px, x, y)
                    )
                    f.write(
                        '  <circle cx="%.1f" cy="%.1f" r="%.1f" '
                        'fill="none" stroke="%s" stroke-width="2" />\n'
                        % (x, y, circle_r_px, stroke_color)
                    )
                else:
                    f.write(
                        '  <circle cx="%.1f" cy="%.1f" r="%.1f" '
                        'fill="%s" stroke="%s" stroke-width="2" />\n'
                        % (x, y, circle_r_px, fill_color, stroke_color)
                    )
                label_y = y + circle_r_px * 0.5
                f.write(
                    '  <text x="%.1f" y="%.1f" '
                    'fill="%s" font-family="sans-serif" '
                    'font-size="17.6" text-anchor="middle">%s</text>\n'
                    % (x, label_y, text_color, letter)
                )

        # Structure 1
        f.write('<g id="structure1">\n')
        draw_backbone(residues1, ca1_svg, draw_residues1, unique_residues1)
        draw_common(circ1_svg)
        draw_unique_residue_contacts(unique_residues1, contacts1, circ1_svg)
        draw_unique(unique1_keys, circ1_svg)
        draw_residues(residues1, circ1_svg, draw_residues1, unique_residues1)
        f.write('</g>\n')

        # Structure 2
        f.write('<g id="structure2">\n')
        draw_backbone(residues2, ca2_svg, draw_residues2, unique_residues2)
        draw_common(circ2_svg)
        draw_unique_residue_contacts(unique_residues2, contacts2, circ2_svg)
        draw_unique(unique2_keys, circ2_svg)
        draw_residues(residues2, circ2_svg, draw_residues2, unique_residues2)
        f.write('</g>\n')

        # Structure labels
        center1_x_world = (minx1 + maxx1) / 2.0 + offset1_x
        center2_x_world = (minx2 + maxx2) / 2.0 + offset2_x
        center1_x = center1_x_world * svg_scale
        center2_x = center2_x_world * svg_scale
        label_y = label_y_world * svg_scale

        f.write('<g id="structure_labels">\n')
        f.write(
            '  <text x="%.1f" y="%.1f" text-anchor="middle" '
            'font-family="sans-serif" font-size="17.6" font-weight="bold" '
            'fill="black">%s</text>\n'
            % (center1_x, label_y, label1)
        )
        f.write(
            '  <text x="%.1f" y="%.1f" text-anchor="middle" '
            'font-family="sans-serif" font-size="17.6" font-weight="bold" '
            'fill="black">%s</text>\n'
            % (center2_x, label_y, label2)
        )
        f.write('</g>\n')

        # Legend (slightly larger font)
        legend_x = svg_margin * svg_scale
        legend_y = legend_y_world * svg_scale
        legend_text_offset = circle_r_px * 3.0
        legend_text_dy = circle_r_px * 0.3

        f.write('<g id="legend">\n')
        f.write(
            '  <text x="%.1f" y="%.1f" text-anchor="start" '
            'font-family="sans-serif" font-size="19.2" font-weight="bold" '
            'fill="black">Legend</text>\n'
            % (legend_x, legend_y)
        )
        row_y = legend_y + row_step_px * 0.9

        def legend_circle_row(fill_color, stroke_color, label):
            nonlocal row_y
            cx = legend_x
            cy = row_y
            f.write(
                '  <circle cx="%.1f" cy="%.1f" r="%.1f" '
                'fill="%s" stroke="%s" stroke-width="2" />\n'
                % (cx, cy, circle_r_px, fill_color, stroke_color)
            )
            text_y = cy + legend_text_dy
            f.write(
                '  <text x="%.1f" y="%.1f" text-anchor="start" '
                'font-family="sans-serif" font-size="19.2" '
                'fill="black">%s</text>\n'
                % (legend_x + legend_text_offset, text_y, label)
            )
            row_y += row_step_px

        def legend_line_row(color, width, label):
            nonlocal row_y
            cy = row_y
            line_center_x = legend_x
            line_x1 = line_center_x - circle_r_px * 1.5
            line_x2 = line_center_x + circle_r_px * 1.5
            f.write(
                '  <line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" '
                'stroke="%s" stroke-width="%.1f" stroke-linecap="butt" />\n'
                % (line_x1, cy, line_x2, cy, color, width)
            )
            text_y = cy + legend_text_dy
            f.write(
                '  <text x="%.1f" y="%.1f" text-anchor="start" '
                'font-family="sans-serif" font-size="19.2" '
                'fill="black">%s</text>\n'
                % (legend_x + legend_text_offset, text_y, label)
            )
            row_y += row_step_px

        legend_circle_row("white", "black", "common residue")
        legend_circle_row("white", "#D3D3D3", "extra residue")
        legend_circle_row("#FF5342", "black", "mutated residue")
        legend_circle_row("#FF9463", "black", "right-left flip")

        legend_line_row("black", 6.0, "backbone (CA-CA)")
        legend_line_row("#FF5342", 6.0, "unique intra contact")
        legend_line_row("#FF9463", 6.0, "unique inter contact")
        legend_line_row("#F2C14E", 6.0, "extra contact")
        legend_line_row("#D3D3D3", 3.0, "common contact (Z-offset=0)")
        legend_line_row("#4C7D8E", 6.0, "common contact (Z-offset=1)")
        legend_line_row("#174D7C", 6.0, "common contact (Z-offset=2)")
        legend_line_row("#1C284E", 6.0, "common contact (Z-offset>=3)")
        f.write("</g>\n")

        # Summary with aligned value column using tspans and dynamic value column
        summary_x = legend_x + legend_text_offset * 10.0
        summary_y = legend_y

        f.write('<g id="pf_stats">\n')
        f.write(
            '  <text x="%.1f" y="%.1f" text-anchor="start" '
            'font-family="sans-serif" font-size="19.2" font-weight="bold" '
            'fill="black">Summary</text>\n'
            % (summary_x, summary_y)
        )

        # Blank line between Summary header and first block
        summary_y += row_step_px * 1.3

        base_world = summary_x / svg_scale
        indent_dx_world = 1.5

        font_size_px = 19.2
        char_width_px = font_size_px * 0.6
        char_width_world = char_width_px / svg_scale

        def escape_text(s):
            return s.replace("&", "&amp;").replace("<", "&lt;")

        processed_lines = []
        max_label_end_world = base_world

        for indent_level, line, is_bold in pf_stats_lines:
            if not line:
                processed_lines.append((indent_level, "", is_bold, None))
                continue

            colon_pos = line.find(":")
            if colon_pos != -1:
                left = line[: colon_pos + 1]
                right = line[colon_pos + 1:].strip()
                if right:
                    label = left
                    value = right
                else:
                    label = line
                    value = None
            else:
                label = line
                value = None

            label_x_world = base_world + indent_dx_world * float(indent_level)
            label_end_world = label_x_world + char_width_world * float(len(label))
            if label_end_world > max_label_end_world:
                max_label_end_world = label_end_world

            processed_lines.append((indent_level, label, is_bold, value))

        value_col_world = max_label_end_world + char_width_world * 3.0

        for indent_level, label, is_bold, value in processed_lines:
            if label == "" and value is None:
                summary_y += row_step_px * 0.7
                continue

            label_x_world = base_world + indent_dx_world * float(indent_level)
            label_x = label_x_world * svg_scale
            value_x = value_col_world * svg_scale
            y = summary_y

            font_weight = "bold" if is_bold else "normal"
            label_esc = escape_text(label)
            value_esc = escape_text(value) if value is not None else None

            f.write(
                '  <text x="%.1f" y="%.1f" text-anchor="start" '
                'font-family="sans-serif" font-size="19.2" font-weight="%s" '
                'fill="black">' % (label_x, y, font_weight)
            )
            f.write('<tspan x="%.1f" y="%.1f">%s</tspan>' % (label_x, y, label_esc))
            if value_esc is not None:
                f.write('<tspan x="%.1f" y="%.1f">%s</tspan>' % (value_x, y, value_esc))
            f.write('</text>\n')

            summary_y += row_step_px * 0.6

        f.write("</g>\n")
        f.write("</svg>\n")

    print("SVG file written to:", svg_path)


def main():
    if len(sys.argv) < 3:
        print(
            "Usage: python compare.py <contacts1.csv> <contacts2.csv> "
            "[--offset2 INT] [--flip1] [--rot1 DEG] [--flip2] [--rot2 DEG]"
        )
        sys.exit(1)

    csv1 = sys.argv[1]
    csv2 = sys.argv[2]

    if not os.path.isfile(csv1):
        print("Error: file not found:", csv1)
        sys.exit(1)
    if not os.path.isfile(csv2):
        print("Error: file not found:", csv2)
        sys.exit(1)

    flip1 = False
    flip2 = False
    rot1_deg = 0.0
    rot2_deg = 0.0
    offset2 = 0

    i = 3
    while i < len(sys.argv):
        arg = sys.argv[i]
        if arg == "--offset2":
            if i + 1 >= len(sys.argv):
                print("Error: --offset2 requires an integer value.")
                sys.exit(1)
            offset2 = int(sys.argv[i + 1])
            i += 2
        elif arg == "--flip1":
            flip1 = True
            i += 1
        elif arg == "--flip2":
            flip2 = True
            i += 1
        elif arg == "--rot1":
            if i + 1 >= len(sys.argv):
                print("Error: --rot1 requires a degree value.")
                sys.exit(1)
            rot1_deg = float(sys.argv[i + 1])
            i += 2
        elif arg == "--rot2":
            if i + 1 >= len(sys.argv):
                print("Error: --rot2 requires a degree value.")
                sys.exit(1)
            rot2_deg = float(sys.argv[i + 1])
            i += 2
        else:
            print("Unknown option:", arg)
            print(
                "Usage: python compare.py <contacts1.csv> <contacts2.csv> "
                "[--offset2 INT] [--flip1] [--rot1 DEG] [--flip2] [--rot2 DEG]"
            )
            sys.exit(1)

    root1 = os.path.splitext(csv1)[0]
    root2 = os.path.splitext(csv2)[0]
    svg_path = "%s_vs_%s.svg" % (root1, root2)
    log_path = "%s_vs_%s.log" % (root1, root2)

    orig_stdout = sys.stdout
    log_file = open(log_path, "w")

    class Tee(object):
        def __init__(self, *streams):
            self.streams = streams

        def write(self, data):
            for s in self.streams:
                s.write(data)

        def flush(self):
            for s in self.streams:
                s.flush()

    sys.stdout = Tee(orig_stdout, log_file)

    try:
        print("Loading contacts from:", csv1)
        residues1, contacts1 = load_contacts_csv(csv1, residue_offset=0)
        print("  residues:", len(residues1), "contacts:", len(contacts1))

        print("Loading contacts from:", csv2, "(offset2 = %d)" % offset2)
        residues2, contacts2 = load_contacts_csv(csv2, residue_offset=offset2)
        print("  residues:", len(residues2), "contacts:", len(contacts2))

        keys1 = set(residues1.keys())
        keys2 = set(residues2.keys())
        common_residues = keys1 & keys2
        print("Common residues:", len(common_residues))
        if not common_residues:
            print("No common residues between the two structures after applying offset2.")
            return

        print(
            "Classifying contacts with strict_nonterminal=%.2f, strict_terminal=%.2f, relaxed=%.2f"
            % (strict_cutoff_nonterminal, strict_cutoff_terminal, relaxed_cutoff)
        )
        common_keys, unique1_keys, unique2_keys = classify_contacts(contacts1, contacts2)
        print("Common contacts:", len(common_keys))
        print("Unique contacts in structure 1:", len(unique1_keys))
        print("Unique contacts in structure 2:", len(unique2_keys))

        pf_indices = sorted({pf for (pf, r) in common_residues})

        # Handedness global flip
        apply_handed_flip = False
        if allow_handedness_global_flip:
            direct_mismatches = 0
            flipped_mismatches = 0
            for key in common_residues:
                h1 = normalize_handed_label(residues1[key]["handed"])
                h2 = normalize_handed_label(residues2[key]["handed"])
                if h1 is None or h2 is None:
                    continue
                if h1 != h2:
                    direct_mismatches += 1
                if h1 != flip_handed_label(h2):
                    flipped_mismatches += 1
            print("Handedness mismatches without flip (ignoring blank/NA):", direct_mismatches)
            print("Handedness mismatches with global flip (ignoring blank/NA):", flipped_mismatches)
            if flipped_mismatches < direct_mismatches:
                apply_handed_flip = True
                print("Using GLOBAL FLIP of handedness for structure 2.")
            else:
                print("Using direct handedness comparison (no global flip).")

        mutated = set()
        handed_diff_all = set()
        for key in common_residues:
            r1 = residues1[key]
            r2 = residues2[key]
            if r1["resname"] != r2["resname"]:
                mutated.add(key)

            h1 = normalize_handed_label(r1["handed"])
            h2 = normalize_handed_label(r2["handed"])
            if h1 is None or h2 is None:
                continue
            if apply_handed_flip:
                h2 = flip_handed_label(h2)
            if h1 != h2:
                handed_diff_all.add(key)

        # Exclude PF-end residues from handedness differences
        ends_to_exclude = set()
        for pf in pf_indices:
            resnums = sorted([r for (p, r) in common_residues if p == pf])
            if not resnums:
                continue
            ends_to_exclude.add((pf, resnums[0]))
            ends_to_exclude.add((pf, resnums[-1]))
        handed_diff = handed_diff_all - ends_to_exclude

        print("Mutated residues:", len(mutated))
        print("Handedness differences (excluding PF-end residues):", len(handed_diff))

        residues_with_intra_contact_change = set()
        residues_with_inter_contact_change = set()

        def add_residues_for_key(key, target_set):
            pf, res, partner_pf, partner_res, ctype = key
            r1k = (pf, res)
            r2k = (partner_pf, partner_res)
            if r1k in common_residues:
                target_set.add(r1k)
            if r2k in common_residues:
                target_set.add(r2k)

        for key in unique1_keys | unique2_keys:
            if key[4] == "intra":
                add_residues_for_key(key, residues_with_intra_contact_change)
            else:
                add_residues_for_key(key, residues_with_inter_contact_change)

        print("Residues with intra-protofilament contact changes:", len(residues_with_intra_contact_change))
        print("Residues with inter-protofilament contact changes:", len(residues_with_inter_contact_change))

        # Choose global sign flip of Z-offsets (structure 2) to reduce non-zero deltas overall (intra+inter)
        apply_zoffset_sign_flip = False
        if allow_zoffset_global_flip:
            count_direct = 0
            count_flipped = 0
            for key in common_keys:
                info1 = contacts1.get(key)
                info2 = contacts2.get(key)
                if key[4] == "inter":
                    b_direct = compute_zoffset_bucket_inter(info1, info2, False, 0)
                    b_flip = compute_zoffset_bucket_inter(info1, info2, True, 0)
                else:
                    b_direct = compute_zoffset_bucket_intra(info1, info2, False)
                    b_flip = compute_zoffset_bucket_intra(info1, info2, True)
                if abs(b_direct) >= 1:
                    count_direct += 1
                if abs(b_flip) >= 1:
                    count_flipped += 1
            print("Common contacts with non-zero rounded Z-offset (direct):", count_direct)
            print("Common contacts with non-zero rounded Z-offset (sign-flipped):", count_flipped)
            if count_flipped < count_direct:
                apply_zoffset_sign_flip = True
                print("Using GLOBAL SIGN FLIP of Z-offsets for structure 2.")
            else:
                print("Using direct Z-offset comparison (no global sign flip).")

        # Choose best inter shift to minimize different inter Z-offsets
        inter_shift = 0
        if allow_inter_zoffset_global_shift:
            inter_shift, best_nonzero, best_total = choose_best_inter_shift(
                common_keys, contacts1, contacts2, apply_zoffset_sign_flip, inter_zoffset_shift_candidates
            )
            if best_total > 0:
                print(
                    "Best inter Z-offset shift for structure 2:",
                    inter_shift,
                    "(non-zero deltas:",
                    best_nonzero,
                    "of",
                    best_total,
                    ")"
                )

        # Per-protofilament (intra only) Z-offset residue set
        residues_with_nonzero_intra_offset = set()
        for key in common_keys:
            if key[4] != "intra":
                continue
            b = compute_zoffset_bucket_intra(contacts1.get(key), contacts2.get(key), apply_zoffset_sign_flip)
            if abs(b) >= 1:
                residues_with_nonzero_intra_offset.add((key[0], key[1]))
                residues_with_nonzero_intra_offset.add((key[2], key[3]))

        # Per-interface (inter) non-zero Z-offset inter-contact count helper
        def inter_contact_nonzero_for_key(key):
            b = compute_zoffset_bucket_inter(contacts1.get(key), contacts2.get(key), apply_zoffset_sign_flip, inter_shift)
            return abs(b) >= 1

        # Build Summary lines (indent_level, line, is_bold)
        pf_stats_lines = []

        # Protofilament blocks
        for pf in pf_indices:
            res_pf_common = {(p, r) for (p, r) in common_residues if p == pf}
            y = len(res_pf_common)
            if y == 0:
                continue

            res_pf_left = {(p, r) for (p, r) in residues1.keys() if p == pf}
            res_pf_right = {(p, r) for (p, r) in residues2.keys() if p == pf}
            aa = len(res_pf_left)
            bb = len(res_pf_right)
            cc = y

            mut_pf = mutated & res_pf_common
            hand_pf = handed_diff & res_pf_common
            intra_pf = residues_with_intra_contact_change & res_pf_common
            zoff_pf = residues_with_nonzero_intra_offset & res_pf_common

            d = len(zoff_pf)
            e = len(mut_pf)
            f_ = len(hand_pf)
            g = len(intra_pf)

            frac_z = 100 * float(d) / float(y)
            frac_mut = 100 * float(e) / float(y)
            frac_hand = 100 * float(f_) / float(y)
            frac_intra = 100 * float(g) / float(y)

            xy_set = set()
            # Don't count mutations in the overall distances!
            #xy_set |= mut_pf
            xy_set |= hand_pf
            xy_set |= intra_pf

            xyz_set = set(xy_set)
            xyz_set |= zoff_pf

            frac_xy = float(100 * len(xy_set)) / float(y)
            frac_xyz = float(100 * len(xyz_set)) / float(y)

            # Adjusted APD to account for residues ordered in only one structure
            # Nuniq = max(#unique residues in left, #unique residues in right)
            nuniq_left = max(0, aa - cc)
            nuniq_right = max(0, bb - cc)
            Nuniq = max(nuniq_left, nuniq_right)
            adj_frac_xy = float(100 * (len(xy_set) + Nuniq)) / float(y + Nuniq)
            adj_frac_xyz = float(100 * (len(xyz_set) + Nuniq)) / float(y + Nuniq)

            header = "Protofilament %d" % pf

            # Console
            print(header)
            print("  - residues in structure 1: %d" % aa)
            print("  - residues in structure 2: %d" % bb)
            print("  - residues in common: %d" % cc)
            print("  - max(extra residues): %d" % Nuniq)
            print("  For common residues:")
            print("  - right-left flips: %d" % f_)
            print("  - uniq intra-contacts: %d" % g)
            print("  - Z-offsets: %s" % d)
            print("  - XY-differences: %d" % len(xy_set))
            print("  - XYZ-differences: %d" % len(xyz_set))
            print("  -----------------------------------------------")
            print("  - PF%d Amyloid Packing Difference (XY): %.0f%%" % (pf, adj_frac_xy) )
            print("  - PF%d Amyloid Packing Difference (XYZ): %.0f%%" % (pf, adj_frac_xyz) )
            print("")

            # SVG Summary
            pf_stats_lines.append((0, header, True))
            pf_stats_lines.append((1, "- residues in structure 1: %d" % aa, False))
            pf_stats_lines.append((1, "- residues in structure 2: %d" % bb, False))
            pf_stats_lines.append((1, "- residues in common: %d" % cc, False))
            pf_stats_lines.append((1, "- max(extra residues): %d" % Nuniq, False))
            pf_stats_lines.append((1, "For common residues:", False))
            pf_stats_lines.append((1, "- right-left flips: %d" % f_, False))
            pf_stats_lines.append((1, "- uniq intra-contacts: %d" % g, False))
            pf_stats_lines.append((1, "- Z-offsets: %d" % d, False))
            pf_stats_lines.append((1, "- XY-differences: %d" % len(xy_set), False))
            pf_stats_lines.append((1, "- XYZ-differences: %d" % len(xyz_set), False))
            pf_stats_lines.append((1, "-----------------------------------------------", False))
            pf_stats_lines.append((1, "- PF%d Amyloid Packing Difference (XY): %.0f%%" % (pf, adj_frac_xy), True))
            pf_stats_lines.append((1, "- PF%d Amyloid Packing Difference (XYZ): %.0f%%" % (pf, adj_frac_xyz), True))
            pf_stats_lines.append((0, "", False))

        # Interface blocks
        for idx_i in range(len(pf_indices)):
            for idx_j in range(idx_i + 1, len(pf_indices)):
                pf_i = pf_indices[idx_i]
                pf_j = pf_indices[idx_j]

                # Residues involved in inter contacts for this interface in left/right.
                # IMPORTANT: interface membership is defined by STRICT contacts only
                # (<= 4.5 Å), consistent with the contact classification rules.
                inter_res_left = set()
                for key, info in contacts1.items():
                    if key[4] != "inter":
                        continue
                    if {key[0], key[2]} != {pf_i, pf_j}:
                        continue
                    if info is None or info.get("distance") is None:
                        continue
                    if info["distance"] > strict_cutoff_nonterminal:
                        continue
                    inter_res_left.add((key[0], key[1]))
                    inter_res_left.add((key[2], key[3]))

                inter_res_right = set()
                for key, info in contacts2.items():
                    if key[4] != "inter":
                        continue
                    if {key[0], key[2]} != {pf_i, pf_j}:
                        continue
                    if info is None or info.get("distance") is None:
                        continue
                    if info["distance"] > strict_cutoff_nonterminal:
                        continue
                    inter_res_right.add((key[0], key[1]))
                    inter_res_right.add((key[2], key[3]))

                a = len(inter_res_left)
                b = len(inter_res_right)
                inter_res_common = inter_res_left & inter_res_right
                c = len(inter_res_common)

                # Extra residues for this interface (residues participating in inter-PF contacts
                # that are present only in one structure)
                extra_left = max(0, a - c)
                extra_right = max(0, b - c)
                Nextra = max(extra_left, extra_right)

                # d_contact: number of common inter contacts with non-zero Z-offset delta (after best shift)
                d_contact = 0
                zoff_res_inter = set()
                for key in common_keys:
                    if key[4] != "inter":
                        continue
                    if {key[0], key[2]} != {pf_i, pf_j}:
                        continue
                    if inter_contact_nonzero_for_key(key):
                        d_contact += 1
                        zoff_res_inter.add((key[0], key[1]))
                        zoff_res_inter.add((key[2], key[3]))

                # e_res: right-left flipped residues involved in inter contacts and present in common interface residue set
                inout_res_inter = handed_diff & inter_res_common
                e_res = len(inout_res_inter)

                # f_res: number of COMMON interface residues that participate in at least
                # one UNIQUE inter-PF contact for this interface.
                uniq_res_inter = set()
                for key in (unique1_keys | unique2_keys):
                    if key[4] != "inter":
                        continue
                    if {key[0], key[2]} != {pf_i, pf_j}:
                        continue
                    uniq_res_inter.add((key[0], key[1]))
                    uniq_res_inter.add((key[2], key[3]))
                f_res = len(uniq_res_inter & inter_res_common)

                header_if = "Protofilament interface %d-%d" % (pf_i, pf_j)

                if c == 0:
                    z_line = "- Z-offsets (%d/%d): NA" % (d_contact, c)
                    inout_line = "- right-left flips (%d/%d): NA" % (e_res, c)
                    uniq_line = "- uniq inter-contacts (%d/%d): NA" % (f_res, c)
                    xy_line = "-  PF%d-PF%d Amyloid Packing Difference (XY): NA" % (pf_i, pf_j)
                    xyz_line = "-  PF%d-PF%d Amyloid Packing Difference (XYZ): NA" % (pf_i, pf_j)
                else:
                    frac_z = 100 * float(d_contact) / float(c)
                    frac_inout = 100 * float(e_res) / float(c)
                    frac_uniq = 100 * float(f_res) / float(c)

                    xy_set = set()
                    xy_set |= inout_res_inter
                    xy_set |= (uniq_res_inter & inter_res_common)

                    xyz_set = set(xy_set)
                    xyz_set |= (zoff_res_inter & inter_res_common)

                    frac_xy = 100 * float(len(xy_set)) / float(c)
                    frac_xyz = 100 * float(len(xyz_set)) / float(c)

                    # Adjusted APD to account for extra interface residues present in only one structure
                    adj_frac_xy = 100 * float(len(xy_set) + Nextra) / float(c + Nextra)
                    adj_frac_xyz = 100 * float(len(xyz_set) + Nextra) / float(c + Nextra)

                    z_line = "- Z-offsets (%d/%d): %.1f%%" % (d_contact, c, frac_z)
                    inout_line = "- right-left flips (%d/%d): %.1f%%" % (e_res, c, frac_inout)
                    uniq_line = "- uniq inter-contacts (%d/%d): %.1f%%" % (f_res, c, frac_uniq)
                    xy_line =  "- PF%d-PF%d Amyloid Packing Difference (XY): %.1f%%" % (pf_i, pf_j, adj_frac_xy)
                    xyz_line = "- PF%d-PF%d Amyloid Packing Difference (XYZ): %.1f%%" % (pf_i, pf_j, adj_frac_xyz)

                # Console
                print(header_if)
                print("  - residues in left:", a)
                print("  - residues in right:", b)
                print("  - residues in common:", c)
                print("  - max(extra residues):", Nextra)
                print("  For common residues:")
                print("  " + z_line)
                print("  " + inout_line)
                print("  " + uniq_line)
                print("  ---------------------------------------------")
                print("  " + xy_line)
                print("  " + xyz_line)
                print("")

                # SVG Summary
                pf_stats_lines.append((0, header_if, True))
                pf_stats_lines.append((1, "- residues in left: %d" % a, False))
                pf_stats_lines.append((1, "- residues in right: %d" % b, False))
                pf_stats_lines.append((1, "- residues in common: %d" % c, False))
                pf_stats_lines.append((1, "- max(extra residues): %d" % Nextra, False))
                pf_stats_lines.append((1, "For common residues:", False))
                pf_stats_lines.append((1, z_line, False))
                pf_stats_lines.append((1, inout_line, False))
                pf_stats_lines.append((1, uniq_line, False))
                pf_stats_lines.append((1, "---------------------------------------------", False))
                pf_stats_lines.append((1, xy_line, True))
                pf_stats_lines.append((1, xyz_line, True))
                pf_stats_lines.append((0, "", False))

        print("Writing SVG:", svg_path)
        write_svg(
            svg_path,
            residues1,
            residues2,
            contacts1,
            contacts2,
            common_keys,
            unique1_keys,
            unique2_keys,
            mutated,
            handed_diff,
            common_residues,
            flip1,
            rot1_deg,
            flip2,
            rot2_deg,
            apply_zoffset_sign_flip,
            inter_shift,
            pf_stats_lines,
        )

        print("Done.")

    finally:
        sys.stdout = orig_stdout
        log_file.close()
        print("Log file written to:", log_path)


if __name__ == "__main__":
    main()

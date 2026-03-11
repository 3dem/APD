#!/usr/bin/env python3
"""
helix.py

Generate helical symmetry copies for an amyloid structure using Biopython.

Changes in this version:
  - Only one positional argument: input PDB filename.
  - Output is always written as CIF: <inputname>_helix.cif
  - When box size and pixel size are given, the computed XY coordinates
    of the helical axis are printed to the screen.
  - The output file contains only:
        - the original chains selected for copying, and
        - the symmetry-generated copies of those chains.
    All other chains from the input model are removed.
"""

import argparse
import math
import sys
import os
from copy import deepcopy

from Bio.PDB import PDBParser, MMCIFParser, PDBIO, MMCIFIO


def parse_args():
    parser = argparse.ArgumentParser(
        description="Apply helical symmetry to selected chains of an amyloid structure."
    )

    parser.add_argument("input", help="Input structure file (PDB format).")

    parser.add_argument(
        "--chains",
        required=True,
        help="Comma-separated list of chain IDs to replicate, for example A,C,D",
    )
    parser.add_argument(
        "--twist",
        type=float,
        required=True,
        help="Helical twist per subunit in degrees (RELION convention)",
    )
    parser.add_argument(
        "--rise",
        type=float,
        required=True,
        help="Helical rise per subunit in Angstroms",
    )
    parser.add_argument(
        "--n_copies",
        default=10,
        type=int,
        required=False,
        help=(
            "Number of copies to generate in each direction along the helical axis "
            "(total copies = original + 2 * n_copies)"
        ),
    )
    parser.add_argument(
        "--output",
        default=None,
        type=str,
        required=False,
        help=(
            "Name of the output file, add _helix.cif to input if none provided"
        ),
    )

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--axis_xy",
        nargs=2,
        metavar=("AXIS_X", "AXIS_Y"),
        type=float,
        help="XY position (Angstroms) of the helical axis",
    )
    group.add_argument(
        "--box_and_pixel",
        nargs=2,
        metavar=("BOX_SIZE_PIX", "PIXEL_SIZE_A"),
        help=(
            "Box size (pixels) and pixel size (Angstroms). "
            "Helical axis is placed at the center of pixel int(boxsize/2)."
        ),
    )

    return parser.parse_args()


def compute_axis_xy(args):
    """
    Compute XY position of the helical axis.

    Either:
      - directly from --axis_xy, or
      - from RELION-style box size and pixel size:
            center_pixel = int(boxsize / 2)
            axis_x = axis_y = (center_pixel + 0.5 )* pixel_size

    When box_and_pixel is used, the computed axis coordinates are printed.
    """
    if args.axis_xy is not None:
        ax, ay = args.axis_xy
        return ax, ay

    # RELION-style box size and pixel size
    box_size_pix = int(args.box_and_pixel[0])
    pixel_size = float(args.box_and_pixel[1])

    center_pixel = int(box_size_pix / 2)
    axis_pos = (center_pixel + 0.5) * pixel_size
    axis_x = axis_pos
    axis_y = axis_pos

    print(
        "Computed helical axis position from box size and pixel size: "
        "axis_x = {:.3f} A, axis_y = {:.3f} A".format(axis_x, axis_y)
    )

    return axis_x, axis_y


def generate_chain_id_pool(existing_ids):
    """
    Generate a pool of chain IDs (single-character) not used in existing_ids.
    Uses A-Z, a-z, 0-9.
    """
    candidates = (
        [chr(i) for i in range(ord("A"), ord("Z") + 1)]
        + [chr(i) for i in range(ord("a"), ord("z") + 1)]
        + [str(i) for i in range(10)]
    )
    return [c for c in candidates if c not in existing_ids]


def apply_helical_transform_to_chain(chain, n_step, twist_deg, rise, axis_x, axis_y):
    """
    Apply helical symmetry transformation to a chain in place.

    Transformation:
      - rotation by (n_step * twist_deg) around the z-axis passing through (axis_x, axis_y)
      - translation by (n_step * rise) along z
    """
    angle_rad = math.radians(n_step * twist_deg)
    cos_a = math.cos(angle_rad)
    sin_a = math.sin(angle_rad)
    dz = n_step * rise

    for atom in chain.get_atoms():
        x, y, z = atom.coord

        # Translate so axis is at origin
        x_shift = x - axis_x
        y_shift = y - axis_y

        # Rotate around z-axis
        x_rot = x_shift * cos_a - y_shift * sin_a
        y_rot = x_shift * sin_a + y_shift * cos_a
        z_rot = z + dz

        # Translate back
        atom.set_coord((x_rot + axis_x, y_rot + axis_y, z_rot))


def main():
    args = parse_args()

    if args.output is None:
        output_path = os.path.splitext(os.path.basename(args.input))[0] + "_helix.cif"
    else:
        output_path = args.output

    # Parse input structure

    ext = os.path.splitext(args.input)[1].lower()
    if ext in [".cif", ".mmcif"]:
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)

    structure = parser.get_structure("input_structure", args.input)

    try:
        model = next(structure.get_models())
    except StopIteration:
        print("No models found in the input structure.", file=sys.stderr)
        sys.exit(1)

    selected_chain_ids = [c.strip() for c in args.chains.split(",") if c.strip()]

    # Check existing chains
    existing_chain_ids = {chain.id for chain in model.get_chains()}

    for cid in selected_chain_ids:
        if cid not in existing_chain_ids:
            print(
                "ERROR: Chain {} not found in the input structure.".format(cid),
                file=sys.stderr,
            )
            sys.exit(1)

    if args.n_copies < 1:
        print("ERROR: --n_copies must be at least 1.", file=sys.stderr)
        sys.exit(1)

    axis_x, axis_y = compute_axis_xy(args)

    # Extract selected chains
    selected_chains = {chain.id: chain for chain in model.get_chains()
                       if chain.id in selected_chain_ids}

    # Remove all chains from the model
    for chain in list(model.get_chains()):
        model.detach_child(chain.id)

    # Re-add selected chains in the order given by --chains
    for cid in selected_chain_ids:
        model.add(selected_chains[cid])

    # Rebuild existing_chain_ids after removal
    existing_chain_ids = {chain.id for chain in model.get_chains()}

    # Determine chain IDs for new copies
    total_needed = len(selected_chain_ids) * 2 * args.n_copies
    chain_pool = generate_chain_id_pool(existing_chain_ids)

    if len(chain_pool) < total_needed:
        print(
            "ERROR: Not enough chain IDs available ({} needed, {} available).".format(
                total_needed, len(chain_pool)
            ),
            file=sys.stderr,
        )
        sys.exit(1)

    # Map original chain IDs (only selected ones remain) to Chain objects
    chain_map = {chain.id: chain for chain in model.get_chains()}

    twist_deg = args.twist
    rise = args.rise

    # Create symmetry copies
    for step in range(1, args.n_copies + 1):
        for direction in (+1, -1):
            n_step = step * direction

            for cid in selected_chain_ids:
                original_chain = chain_map[cid]
                new_chain = deepcopy(original_chain)

                new_id = chain_pool.pop(0)
                new_chain.id = new_id

                apply_helical_transform_to_chain(
                    new_chain, n_step, twist_deg, rise, axis_x, axis_y
                )

                model.add(new_chain)

    # Save PDB/CIF output (structure contains only selected chains and their copies)
    ext = os.path.splitext(output_path)[1].lower()

    if ext == ".pdb":
        io = PDBIO()
    elif ext in (".cif", ".mmcif"):
        io = MMCIFIO()
    else:
        raise ValueError(
            f"Unsupported output format '{ext}'. "
            "Use .pdb, .cif, or .mmcif"
        )

    io.set_structure(structure)
    io.save(output_path)

    total_new_chains = len(selected_chain_ids) * 2 * args.n_copies
    print(
        "Wrote helical symmetry expanded structure to {} "
        "with {} new chains ({} copies in +Z and {} in -Z for each selected chain).".format(
            output_path, total_new_chains, args.n_copies, args.n_copies
        )
    )


if __name__ == "__main__":
    main()

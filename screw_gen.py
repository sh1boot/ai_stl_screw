#!/usr/bin/env python3
"""Generate a matched screw and threaded-hole (nut) STL pair.

Thread profile is flat-flanked (trapezoidal) with configurable crest and root
flat widths.  Each ring vertex samples the thread profile at its own helix
phase, producing a true helical ridge rather than a stack of rotationally-
symmetric discs.

Usage:
    python3 screw_gen.py [options]
    python3 screw_gen.py --help
"""

import argparse
import math
import struct
import sys


# ---------------------------------------------------------------------------
# Binary STL helpers
# ---------------------------------------------------------------------------

def _normal(v0, v1, v2):
    ax, ay, az = v1[0]-v0[0], v1[1]-v0[1], v1[2]-v0[2]
    bx, by, bz = v2[0]-v0[0], v2[1]-v0[1], v2[2]-v0[2]
    nx = ay*bz - az*by
    ny = az*bx - ax*bz
    nz = ax*by - ay*bx
    mag = math.sqrt(nx*nx + ny*ny + nz*nz) or 1.0
    return nx/mag, ny/mag, nz/mag


def write_stl(path, triangles):
    """Write binary STL; triangles is a list of (v0, v1, v2) 3-tuples."""
    buf = bytearray(80)
    buf += struct.pack('<I', len(triangles))
    for v0, v1, v2 in triangles:
        n = _normal(v0, v1, v2)
        buf += struct.pack('<3f', *n)
        for v in (v0, v1, v2):
            buf += struct.pack('<3f', *v)
        buf += b'\x00\x00'
    with open(path, 'wb') as f:
        f.write(buf)


# ---------------------------------------------------------------------------
# Primitive mesh helpers
# ---------------------------------------------------------------------------

def quad_tris(v00, v01, v10, v11):
    """Two CCW triangles for a quad: (v00,v01) bottom edge, (v10,v11) top edge."""
    return [(v00, v10, v11), (v00, v11, v01)]


def fan_cap(centre, ring, reverse=False):
    """Triangle fan from centre to ring, closing a circular opening."""
    tris = []
    n = len(ring)
    for i in range(n):
        j = (i + 1) % n
        if reverse:
            tris.append((centre, ring[j], ring[i]))
        else:
            tris.append((centre, ring[i], ring[j]))
    return tris


# ---------------------------------------------------------------------------
# Thread profile
# ---------------------------------------------------------------------------

def thread_radius(helix_phase, r_inner, r_outer, crest_frac, root_frac):
    """
    Thread surface radius at a given fractional position along one pitch [0, 1).

    Profile (symmetric, flat-flanked):
      root flat  : r_inner  for root_frac/2 at each end of the pitch
      crest flat : r_outer  for crest_frac at the centre of the pitch
      linear flanks in between
    """
    root_half = root_frac  / 2.0
    flank     = 0.5 - root_half - crest_frac / 2.0

    if flank < -1e-9:
        raise ValueError("crest_frac + root_frac must be <= 1.0")
    flank = max(flank, 0.0)

    p = helix_phase % 1.0
    if p > 0.5:
        p = 1.0 - p   # symmetric profile

    if p <= root_half:
        return r_inner
    elif flank > 0 and p <= root_half + flank:
        return r_inner + (r_outer - r_inner) * (p - root_half) / flank
    else:
        return r_outer


# ---------------------------------------------------------------------------
# Helix ring
# ---------------------------------------------------------------------------

def _helix_ring(z, ri, ro, pitch, crest_frac, root_frac, n_sides):
    """
    One ring of n_sides vertices at height z.

    Each vertex i sits at angle  ang_base + 2π·i/n  where
    ang_base = 2π·z/pitch (the helix phase at z).  Its radius is
    thread_radius evaluated at helix_phase = (z/pitch + i/n) % 1,
    i.e. at the fractional pitch position that corresponds to that
    vertex's angular position on the helix.

    This makes only a narrow arc of each ring reach r_outer (the crest),
    and that arc rotates continuously with z — producing a true helix
    instead of a stack of rotationally-symmetric discs.
    """
    ang_base  = 2.0 * math.pi * z / pitch
    z_frac    = (z / pitch) % 1.0
    verts = []
    for i in range(n_sides):
        hp = (z_frac + i / n_sides) % 1.0
        r  = thread_radius(hp, ri, ro, crest_frac, root_frac)
        a  = ang_base + 2.0 * math.pi * i / n_sides
        verts.append((r * math.cos(a), r * math.sin(a), z))
    return verts


# ---------------------------------------------------------------------------
# Screw geometry
# ---------------------------------------------------------------------------

def build_screw(r_inner, r_outer, pitch, crest_frac, root_frac,
                length, taper_mm, tolerance, n_sides, n_thread_steps):
    """
    Return triangle list for a solid threaded screw shaft.

    Oriented along +Z: tip at z=0, head end at z=length.
    taper_mm  : radial reduction at the tip vs the base.
    tolerance : subtracted from both radii (shrinks screw for clearance).
    """
    ri_base = r_inner - tolerance
    ro_base = r_outer - tolerance
    total_steps = max(1, int(math.ceil(length / pitch))) * n_thread_steps

    rings = []
    for s in range(total_steps + 1):
        z = length * s / total_steps
        t  = (1.0 - z / length) if length > 0 else 0.0
        ri = max(0.0, ri_base - taper_mm * t)
        ro = max(ri,  ro_base - taper_mm * t)
        rings.append(_helix_ring(z, ri, ro, pitch, crest_frac, root_frac, n_sides))

    tris = []

    # Lateral surface: reversed winding → outward-facing normals for a solid body.
    for s in range(total_steps):
        r0, r1 = rings[s], rings[s + 1]
        for i in range(n_sides):
            j = (i + 1) % n_sides
            tris += quad_tris(r0[j], r0[i], r1[j], r1[i])

    # End caps
    tris += fan_cap((0.0, 0.0, 0.0),    rings[0],  reverse=True)   # bottom, normal −Z
    tris += fan_cap((0.0, 0.0, length), rings[-1], reverse=False)  # top,    normal +Z

    return tris


# ---------------------------------------------------------------------------
# Nut (threaded hole in a cuboid) geometry
# ---------------------------------------------------------------------------

def build_nut(r_inner, r_outer, pitch, crest_frac, root_frac,
              length, taper_mm, tolerance, n_sides, n_thread_steps, box_margin):
    """
    Return triangle list for a rectangular block with a matching threaded bore.

    The bore is the exact negative of the screw: same helix, tolerances added,
    taper reversed so the larger bore end receives the tapered screw tip.
    box_margin: extra wall thickness beyond r_outer + tolerance on each side.
    """
    ri_base = r_inner + tolerance
    ro_base = r_outer + tolerance
    total_steps = max(1, int(math.ceil(length / pitch))) * n_thread_steps

    bore_rings = []
    for s in range(total_steps + 1):
        z = length * s / total_steps
        t  = (1.0 - z / length) if length > 0 else 0.0   # 1 at tip, 0 at base
        ri = ri_base + taper_mm * t                        # bore widens toward tip
        ro = max(ri, ro_base + taper_mm * t)
        bore_rings.append(
            _helix_ring(z, ri, ro, pitch, crest_frac, root_frac, n_sides))

    tris = []

    # Bore surface: non-reversed winding → inward-facing normals (correct for a
    # cavity whose material boundary faces toward the bore axis).
    for s in range(total_steps):
        r0, r1 = bore_rings[s], bore_rings[s + 1]
        for i in range(n_sides):
            j = (i + 1) % n_sides
            tris += quad_tris(r0[i], r0[j], r1[i], r1[j])

    # The bore surface open ends (bore_rings[0] and bore_rings[-1]) are closed
    # by the square-face holes below — no separate bore end caps needed.

    # ---- Outer cuboid ----
    box = ro_base + box_margin
    z0, z1 = 0.0, length
    corners_xy = [(-box, -box), (box, -box), (box, box), (-box, box)]

    def cv(ci, z):
        x, y = corners_xy[ci]
        return (x, y, z)

    # Four side faces, normals pointing outward
    for i in range(4):
        j = (i + 1) % 4
        v00 = cv(i, z0); v01 = cv(j, z0)
        v10 = cv(i, z1); v11 = cv(j, z1)
        tris += quad_tris(v01, v00, v11, v10)

    # Top and bottom faces: square with a helical-profile hole
    tris += _square_with_hole([cv(i, z1) for i in range(4)], bore_rings[-1], reverse=False)
    tris += _square_with_hole([cv(i, z0) for i in range(4)], bore_rings[0],  reverse=True)

    return tris


def _square_with_hole(box_ring, hole_ring, reverse):
    """
    Triangulate a square face with a hole using radial sector decomposition.

    box_ring : 4 vertices of the square, CCW when viewed from the face normal.
    hole_ring: n vertices of the hole boundary, CCW from the same viewpoint.
    reverse  : flip all triangle windings (for the bottom face).

    For each hole edge, project its endpoints radially onto the box boundary,
    collect any box corners in the same angular sector, and fan-triangulate
    from the inner hole vertex.  No fan diagonal can cross the hole.
    """
    n   = len(hole_ring)
    z   = box_ring[0][2]
    box = max(abs(box_ring[0][0]), abs(box_ring[0][1]))

    def _proj(px, py):
        m = max(abs(px), abs(py))
        if m < 1e-12:
            return (box, 0.0, z)
        s = box / m
        return (px * s, py * s, z)

    corner_angles = [math.atan2(c[1], c[0]) % (2 * math.pi) for c in box_ring]

    def _in_sector(a, a0, a1):
        a = a % (2 * math.pi)
        if a1 > a0:
            return a0 < a < a1
        return a > a0 or a < a1

    def _add(tris, t):
        tris.append((t[0], t[2], t[1]) if reverse else t)

    tris = []
    for j in range(n):
        h0 = hole_ring[j]
        h1 = hole_ring[(j + 1) % n]
        b0 = _proj(h0[0], h0[1])
        b1 = _proj(h1[0], h1[1])
        a0 = math.atan2(h0[1], h0[0]) % (2 * math.pi)
        a1 = math.atan2(h1[1], h1[0]) % (2 * math.pi)

        mid_corners = [box_ring[ci] for ci, ca in enumerate(corner_angles)
                       if _in_sector(ca, a0, a1)]

        poly = [h0, b0] + mid_corners + [b1, h1]
        p0 = poly[0]
        for i in range(1, len(poly) - 1):
            _add(tris, (p0, poly[i], poly[i + 1]))

    return tris


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Generate matched screw and threaded-hole STL files.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    p.add_argument("--outer-radius", type=float, default=5.0, metavar="MM",
                   help="Thread outer (crest) radius")
    p.add_argument("--inner-radius", type=float, default=3.5, metavar="MM",
                   help="Thread inner (root) radius")
    p.add_argument("--pitch",        type=float, default=1.0, metavar="MM",
                   help="Thread pitch (axial distance per full turn)")
    p.add_argument("--crest-frac",   type=float, default=0.25, metavar="F",
                   help="Fraction of pitch that is crest flat [0,1)")
    p.add_argument("--root-frac",    type=float, default=0.25, metavar="F",
                   help="Fraction of pitch that is root flat [0,1)")
    p.add_argument("--length",       type=float, default=20.0, metavar="MM",
                   help="Axial length of screw / threaded bore")
    p.add_argument("--taper",        type=float, default=0.0,  metavar="MM",
                   help="Radial reduction at tip vs base (0 = parallel)")
    p.add_argument("--tolerance",    type=float, default=0.2,  metavar="MM",
                   help="Radial clearance: screw shrunk and hole grown by this")
    p.add_argument("--sides",        type=int,   default=64,   metavar="N",
                   help="Polygon approximation sides (>=6)")
    p.add_argument("--thread-steps", type=int,   default=32,   metavar="N",
                   help="Axial sample steps per thread pitch (>=4)")
    p.add_argument("--box-margin",   type=float, default=3.0,  metavar="MM",
                   help="Wall thickness around bore in the nut block")
    p.add_argument("--screw-out",    type=str,   default="screw.stl",
                   help="Output path for screw STL")
    p.add_argument("--nut-out",      type=str,   default="nut.stl",
                   help="Output path for nut STL")

    return p.parse_args()


def main():
    args = parse_args()

    if args.inner_radius <= 0:
        sys.exit("error: --inner-radius must be positive")
    if args.outer_radius <= args.inner_radius:
        sys.exit("error: --outer-radius must be greater than --inner-radius")
    if args.pitch <= 0:
        sys.exit("error: --pitch must be positive")
    if not (0.0 < args.crest_frac + args.root_frac <= 1.0):
        sys.exit("error: --crest-frac + --root-frac must be in (0, 1]")
    if args.length <= 0:
        sys.exit("error: --length must be positive")
    if args.sides < 6:
        sys.exit("error: --sides must be >= 6")
    if args.thread_steps < 4:
        sys.exit("error: --thread-steps must be >= 4")
    if args.taper < 0:
        sys.exit("error: --taper must be >= 0")
    if args.tolerance < 0:
        sys.exit("error: --tolerance must be >= 0")

    common = dict(
        r_inner        = args.inner_radius,
        r_outer        = args.outer_radius,
        pitch          = args.pitch,
        crest_frac     = args.crest_frac,
        root_frac      = args.root_frac,
        length         = args.length,
        taper_mm       = args.taper,
        n_sides        = args.sides,
        n_thread_steps = args.thread_steps,
    )

    print(f"Generating screw -> {args.screw_out}")
    screw_tris = build_screw(tolerance=args.tolerance, **common)
    write_stl(args.screw_out, screw_tris)
    print(f"  {len(screw_tris):,} triangles")

    print(f"Generating nut   -> {args.nut_out}")
    nut_tris = build_nut(tolerance=args.tolerance, box_margin=args.box_margin, **common)
    write_stl(args.nut_out, nut_tris)
    print(f"  {len(nut_tris):,} triangles")

    print("Done.")


if __name__ == "__main__":
    main()

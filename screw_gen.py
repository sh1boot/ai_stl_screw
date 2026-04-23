#!/usr/bin/env python3
"""Generate a matched screw and threaded-hole (nut) STL pair.

Thread profile is flat-flanked (trapezoidal) with configurable crest and root
flat widths.  The geometry is approximated by a polygonal revolution swept
along the helix axis.

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


def ring_verts(z, radius, n_sides, angle_offset=0.0):
    """Return a list of n_sides (x,y,z) vertices around a circle at height z."""
    verts = []
    for i in range(n_sides):
        a = angle_offset + 2.0 * math.pi * i / n_sides
        verts.append((radius * math.cos(a), radius * math.sin(a), z))
    return verts


def annular_tris(inner, outer, reverse=False):
    """
    Annular face between two co-planar rings of the same length.

    Both rings must be ordered consistently (e.g. both CCW from above).
    reverse=True flips winding (normal points the other way).
    """
    n = len(inner)
    tris = []
    for i in range(n):
        j = (i + 1) % n
        # quad: inner[i], inner[j], outer[j], outer[i]
        t1 = (inner[i], outer[i], outer[j])
        t2 = (inner[i], outer[j], inner[j])
        if reverse:
            tris.append((t1[0], t1[2], t1[1]))
            tris.append((t2[0], t2[2], t2[1]))
        else:
            tris.extend([t1, t2])
    return tris


# ---------------------------------------------------------------------------
# Thread profile
# ---------------------------------------------------------------------------

def thread_radius(z_frac, r_inner, r_outer, crest_frac, root_frac):
    """
    Radial position of the thread surface at fractional position z_frac in [0,1).

    Profile (symmetric, flat-flanked):
      root flat  : r_inner for root_frac/2 at each end of the pitch
      crest flat : r_outer for crest_frac centred in the pitch
      linear flanks in between
    """
    root_half  = root_frac  / 2.0
    crest_half = crest_frac / 2.0
    flank      = 0.5 - root_half - crest_half

    if flank < -1e-9:
        raise ValueError("crest_frac + root_frac must be <= 1.0")
    flank = max(flank, 0.0)

    p = z_frac % 1.0
    if p > 0.5:
        p = 1.0 - p   # mirror: profile is symmetric around 0.5

    if p <= root_half:
        return r_inner
    elif flank > 0 and p <= root_half + flank:
        t = (p - root_half) / flank
        return r_inner + (r_outer - r_inner) * t
    else:
        return r_outer


# ---------------------------------------------------------------------------
# Screw geometry
# ---------------------------------------------------------------------------

def build_screw(r_inner, r_outer, pitch, crest_frac, root_frac,
                length, taper_mm, tolerance, n_sides, n_thread_steps):
    """
    Return triangle list for a solid threaded screw shaft.

    Oriented along +Z: tip at z=0, head end at z=length.
    taper_mm: how much smaller the tip radius is vs the base (both ri and ro
              are reduced proportionally toward z=0).
    tolerance: subtracted from both radii (shrinks the screw for clearance).
    """
    ri_base = r_inner - tolerance
    ro_base = r_outer - tolerance

    total_steps = max(1, int(math.ceil(length / pitch))) * n_thread_steps

    def radii_at_z(z):
        t = (1.0 - z / length) if length > 0 else 0.0  # 1 at tip, 0 at base
        taper = taper_mm * t
        ri = max(0.0, ri_base - taper)
        ro = max(ri,  ro_base - taper)
        return ri, ro

    # Pre-compute all rings.  rings_shaft[s] = inner shaft ring at step s.
    # rings_thread[s] = outer thread-surface ring at step s.
    rings_shaft  = []
    rings_thread = []

    for s in range(total_steps + 1):
        z = length * s / total_steps
        ri, ro = radii_at_z(z)
        ang = 2.0 * math.pi * z / pitch   # helix angle
        z_frac = (z / pitch) % 1.0
        tr = thread_radius(z_frac, ri, ro, crest_frac, root_frac)
        rings_shaft.append( ring_verts(z, ri, n_sides, ang))
        rings_thread.append(ring_verts(z, tr, n_sides, ang))

    tris = []

    for s in range(total_steps):
        rs0, rs1 = rings_shaft[s],  rings_shaft[s + 1]
        rt0, rt1 = rings_thread[s], rings_thread[s + 1]

        for i in range(n_sides):
            j = (i + 1) % n_sides

            # Inner shaft lateral quad
            tris += quad_tris(rs0[i], rs0[j], rs1[i], rs1[j])

            # Outer thread tip lateral quad
            tris += quad_tris(rt0[i], rt0[j], rt1[i], rt1[j])

            # Bottom annular face of this step (faces downward = normal -Z-ish)
            # Only needed for s==0 (the actual bottom of the screw).
            # At intermediate steps this face is shared with the step above;
            # we rely on it being redundant and just omit duplicates by adding
            # only the bottom at s=0 and the top at s=total_steps-1.

        # Bottom annular ring (at s=0 only, facing down)
        if s == 0:
            tris += annular_tris(rs0, rt0, reverse=True)

    # Top annular ring (at last step, facing up)
    rs_top = rings_shaft[-1]
    rt_top = rings_thread[-1]
    tris += annular_tris(rs_top, rt_top, reverse=False)

    # Bottom cap: shaft interior disk facing down
    c_bot = (0.0, 0.0, 0.0)
    tris += fan_cap(c_bot, rings_shaft[0], reverse=True)

    # Top cap: shaft interior disk facing up
    c_top = (0.0, 0.0, length)
    tris += fan_cap(c_top, rings_shaft[-1], reverse=False)

    return tris


# ---------------------------------------------------------------------------
# Nut (threaded hole in a cuboid) geometry
# ---------------------------------------------------------------------------

def build_nut(r_inner, r_outer, pitch, crest_frac, root_frac,
              length, taper_mm, tolerance, n_sides, n_thread_steps, box_margin):
    """
    Return triangle list for a rectangular block with a matching threaded bore.

    The bore is the negative of the screw (tolerance added, taper reversed so
    the larger end of the hole matches the larger end of the screw).
    box_margin: extra wall thickness around the outer thread radius.
    """
    ri_base = r_inner + tolerance
    ro_base = r_outer + tolerance

    total_steps = max(1, int(math.ceil(length / pitch))) * n_thread_steps

    def radii_at_z(z):
        # Taper in nut: hole is wider at the tip end (z=0) to accept the tapered screw
        t = (1.0 - z / length) if length > 0 else 0.0
        taper = taper_mm * t
        ri = ri_base + taper
        ro = max(ri, ro_base + taper)
        return ri, ro

    rings_shaft  = []
    rings_thread = []

    for s in range(total_steps + 1):
        z = length * s / total_steps
        ri, ro = radii_at_z(z)
        ang = 2.0 * math.pi * z / pitch
        z_frac = (z / pitch) % 1.0
        tr = thread_radius(z_frac, ri, ro, crest_frac, root_frac)
        rings_shaft.append( ring_verts(z, ri, n_sides, ang))
        rings_thread.append(ring_verts(z, tr, n_sides, ang))

    tris = []

    # Bore surfaces: normals face into the bore (toward the axis), so windings
    # are reversed compared to the equivalent screw surfaces.
    for s in range(total_steps):
        rs0, rs1 = rings_shaft[s],  rings_shaft[s + 1]
        rt0, rt1 = rings_thread[s], rings_thread[s + 1]

        for i in range(n_sides):
            j = (i + 1) % n_sides
            # Inner bore wall (reversed winding = normal points away from axis)
            tris += quad_tris(rs0[j], rs0[i], rs1[j], rs1[i])
            # Outer thread tip (reversed)
            tris += quad_tris(rt0[j], rt0[i], rt1[j], rt1[i])

        # Bottom annular of bore (s=0): faces upward into the bore
        if s == 0:
            tris += annular_tris(rings_shaft[0], rings_thread[0], reverse=False)

    # Top annular of bore: faces downward into the bore
    tris += annular_tris(rings_shaft[-1], rings_thread[-1], reverse=True)

    # ---- Outer cuboid ----
    box = ro_base + box_margin   # half-size of the square cross-section
    z0, z1 = 0.0, length

    # Corner (x,y) pairs, CCW when viewed from above
    corners_xy = [(-box, -box), (box, -box), (box, box), (-box, box)]

    def cv(ci, z):
        x, y = corners_xy[ci]
        return (x, y, z)

    # Four side faces, normals pointing outward
    for i in range(4):
        j = (i + 1) % 4
        v00 = cv(i, z0); v01 = cv(j, z0)
        v10 = cv(i, z1); v11 = cv(j, z1)
        # outer face: winding so normal points away from centre
        tris += quad_tris(v01, v00, v11, v10)

    # Top and bottom faces of the cuboid (square with circular hole)
    box_corners_top = [cv(i, z1) for i in range(4)]
    box_corners_bot = [cv(i, z0) for i in range(4)]
    tris += _square_with_hole(box_corners_top, rings_thread[-1], reverse=False)
    tris += _square_with_hole(box_corners_bot, rings_thread[0],  reverse=True)

    return tris


def _square_with_hole(box_ring, hole_ring, reverse):
    """
    Triangulate a square face with a circular hole using radial sector decomposition.

    box_ring : 4 vertices of the square, CCW when viewed from the face normal direction.
    hole_ring: n vertices of the hole boundary, CCW when viewed from the SAME direction.
    reverse  : flip all triangle windings (for the bottom face).

    Each hole segment defines a radial sector; we project the segment endpoints
    outward to the box boundary and collect any box corners inside the sector,
    then fan-triangulate the sector polygon from the inner hole vertex.  This
    avoids any fan diagonal crossing the hole.
    """
    n   = len(hole_ring)
    z   = box_ring[0][2]
    box = max(abs(box_ring[0][0]), abs(box_ring[0][1]))  # half-size

    def _proj(px, py):
        """Radially project (px,py) onto the ±box square boundary."""
        m = max(abs(px), abs(py))
        if m < 1e-12:
            return (box, 0.0, z)
        s = box / m
        return (px * s, py * s, z)

    # Pre-compute angle of each box corner and normalise to [0, 2π)
    corner_angles = []
    for c in box_ring:
        a = math.atan2(c[1], c[0]) % (2 * math.pi)
        corner_angles.append(a)

    def _in_sector(a, a0, a1):
        """True if angle a lies in the CCW arc from a0 to a1 (exclusive endpoints)."""
        a = a % (2 * math.pi)
        if a1 > a0:
            return a0 < a < a1
        # wrap-around sector
        return a > a0 or a < a1

    def _add(tris, t):
        if reverse:
            tris.append((t[0], t[2], t[1]))
        else:
            tris.append(t)

    tris = []

    for j in range(n):
        h0 = hole_ring[j]
        h1 = hole_ring[(j + 1) % n]
        b0 = _proj(h0[0], h0[1])
        b1 = _proj(h1[0], h1[1])

        a0 = math.atan2(h0[1], h0[0]) % (2 * math.pi)
        a1 = math.atan2(h1[1], h1[0]) % (2 * math.pi)

        # Collect box corners strictly inside this angular sector (CCW from a0 to a1)
        mid_corners = []
        for ci, ca in enumerate(corner_angles):
            if _in_sector(ca, a0, a1):
                mid_corners.append(box_ring[ci])

        # Sector polygon: [h0, b0, <box corners>, b1, h1]
        # Fan from h0 — always valid because h0 can see all outward points in its sector.
        poly = [h0, b0] + mid_corners + [b1, h1]
        p0 = poly[0]
        for i in range(1, len(poly) - 1):
            _add(tris, (p0, poly[i], poly[i + 1]))
        # Closing edge: h1 -> h0 is the hole boundary, no extra triangle needed
        # unless b0 == b1 (degenerate sector), handled implicitly.

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
                   help="Polygon approximation sides (≥6)")
    p.add_argument("--thread-steps", type=int,   default=32,   metavar="N",
                   help="Axial sample steps per thread pitch (≥4)")
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
        r_inner       = args.inner_radius,
        r_outer       = args.outer_radius,
        pitch         = args.pitch,
        crest_frac    = args.crest_frac,
        root_frac     = args.root_frac,
        length        = args.length,
        taper_mm      = args.taper,
        n_sides       = args.sides,
        n_thread_steps= args.thread_steps,
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

#!/usr/bin/env python3

"""
Historical script entry point for the maintained multi_group_ace package.

Running this file delegates to multi_group_ace.cli. The helper implementations
below are archived multigroup teaching/reference code; import the package API
for maintained conversion behavior.

Supported CLI workflows:

* JSON material definitions -> simple multigroup ACE and helper files.
* One legacy ASCII multigroup ACE table -> OpenMC/OpenSn MGXS HDF5.
* Continuous-energy neutron ACE file(s) -> GNDS-named OpenMC neutron HDF5
  files and cross_sections.xml, using the optional OpenMC Python API.
* Recursive multigroup verification conversion, optionally with Gmsh meshes.

Examples:

    python simple_ace_mg_full.py --config "json examples/ex4.json"
    python simple_ace_mg_full.py --ace examples/mg08/99902.02m --h5 mgxs.h5
    python simple_ace_mg_full.py --ace He3.ace U235.ace --outdir openmc_data
    python simple_ace_mg_full.py --verification-root examples --outdir build/verification --make-mesh

ACE routing uses table contents, not filesystem extensions. The bundled
multigroup data cannot recover continuous-energy detail. Neutron conversion
processes all tables, merges compatible temperatures from one evaluation,
and rejects thermal/other particle tables, mixed types, and duplicate rounded
temperatures. --metastable-scheme defaults to mcnp; nndc is also available.
--overwrite is required to replace existing continuous-energy HDF5/XML outputs.
A custom --h5 requires one nuclide, a .h5 suffix, and no explicit --outdir.

Install the optional pinned OpenMC API with:

    python -m pip install -r requirements-ce.txt

No OpenMC transport executable is required for conversion. See README.md and
the technical reference for current interfaces and verification limits.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
from dataclasses import dataclass
from datetime import date
from typing import List, Tuple, Optional

# Physical constants for temperature conversion
KB_EV_PER_K = 8.617333262e-5  # Boltzmann constant in eV/K
MEV_TO_EV = 1.0e6             # 1 MeV = 1e6 eV


# -----------------------------------------------------------------------------
# Utility helpers
# -----------------------------------------------------------------------------

def get_today() -> str:
    """Return YYYY-MM-DD string (matches Perl get_date)."""
    d = date.today()
    return f"{d.year:04d}-{d.month:02d}-{d.day:02d}"


def flatten_matrix(mat: List[List[float]]) -> List[float]:
    """Flatten list-of-lists matrix (row-major)."""
    flat: List[float] = []
    for row in mat:
        flat.extend(row)
    return flat


def expand_to_length(vals: List[float], length: int, label: str) -> List[float]:
    """
    Expand a 0-, 1-, or length-element list to exactly length elements:
    - 0 → all zeros
    - 1 → replicate same value
    - length → unchanged
    Anything else → error.
    """
    if not vals:
        return [0.0] * length
    if len(vals) == 1:
        return [float(vals[0])] * length
    if len(vals) != length:
        raise ValueError(
            f"{label}: expected 0, 1, or {length} values, got {len(vals)}"
        )
    return [float(x) for x in vals]


def ensure_descending(energies: List[float]) -> List[float]:
    """Return energies in strictly descending order; reverse if needed."""
    if len(energies) < 2:
        raise ValueError("Energy grid must have at least 2 values (NG+1).")
    if energies[0] < energies[-1]:
        energies = list(reversed(energies))
    for i in range(len(energies) - 1):
        if energies[i] <= energies[i + 1]:
            raise ValueError(
                f"Energy grid must be strictly decreasing; "
                f"got {energies[i]} <= {energies[i+1]}"
            )
    return [float(e) for e in energies]


# -----------------------------------------------------------------------------
# Equiprobable µ-bin generator (ported from Perl get_equiprob_bins)
# -----------------------------------------------------------------------------

def get_equiprob_bins(nbins: int, N: int, A: List[float]) -> List[float]:
    """
    Return boundaries µ_k in [-1,1] such that each interval has equal probability
    for the Legendre-expansion angular distribution with coefficients A[0..N].
    """
    nx = 1000
    dx = 2.0 / nx

    x_vals: List[float] = []
    cdf: List[float] = []

    for k in range(nx + 1):
        x = -1.0 + k * dx
        x_vals.append(x)

        val = 1.0 + x
        if N > 0:
            P: List[float] = [0.0] * (N + 2)
            P[0] = 1.0
            P[1] = x
            for j in range(2, N + 2):
                P[j] = ((2.0 * j - 1.0) * x * P[j - 1] - (j - 1.0) * P[j - 2]) / j
            for j in range(1, N + 1):
                val += A[j] * (P[j + 1] - P[j - 1])
        val *= 0.5
        cdf.append(val)

    prob = 1.0 / nbins
    mu: List[float] = [-1.0] * (nbins + 1)

    for k in range(1, nbins + 1):
        pp = prob * k
        for j in range(nbins, -1, -1):
            if cdf[j] < pp:
                xj = x_vals[j]
                xj1 = x_vals[j + 1]
                cj = cdf[j]
                cj1 = cdf[j + 1]
                mu[k] = xj + (pp - cj) * (xj1 - xj) / (cj1 - cj)
                break

    return mu


# -----------------------------------------------------------------------------
# Data container
# -----------------------------------------------------------------------------

@dataclass
class MGData:
    zaid: str
    za: int
    awr: float
    tmp: float         # kT in MeV
    comment: str
    NG: int
    erg: List[float]   # length NG+1, descending
    sigt: List[float]  # total
    sigc: List[float]  # capture
    sigf: List[float]  # fission
    nu: List[float]
    chi: List[float]
    sigs: List[float]  # scatter matrix, row-major, row=to, col=from
    sig1: List[float]  # P1 scatter, row-major, row=to, col=from
    nbins: int         # angular bins if anisotropic
    comment_short: str
    date: str


# -----------------------------------------------------------------------------
# JSON loader
# -----------------------------------------------------------------------------

def load_from_json(path: str) -> Tuple[MGData, str]:
    with open(path, "r") as f:
        cfg = json.load(f)

    zaid = cfg.get("zaid")
    if not zaid:
        raise ValueError("JSON must contain 'zaid' (e.g., '99902.01m').")
    if not zaid.endswith("m"):
        raise ValueError("zaid must end in 'm' (e.g., '99902.01m').")

    za_str = zaid.split(".")[0]
    try:
        za = int(za_str)
    except ValueError:
        raise ValueError(f"Could not parse ZA from zaid '{zaid}'.")

    awr = float(cfg.get("awr", 1.0))
    tmp = float(cfg.get("tmp", 2.5301e-8))
    if tmp > 1.0:
        tmp *= 8.617385e-11  # MeV/K

    comment = str(cfg.get("comment", "multigroup ACE file"))
    comment_short = comment[:70]

    e_grid = (
        cfg.get("e_bounds")
        or cfg.get("erg")
        or cfg.get("energies")
        or cfg.get("energy_grid")
    )
    if not e_grid:
        raise ValueError("JSON must contain 'e_bounds' (or equivalent) array.")

    erg = ensure_descending([float(x) for x in e_grid])

    NG = int(cfg.get("groups", cfg.get("NG", len(erg) - 1)))
    if len(erg) != NG + 1:
        raise ValueError(
            f"Energy grid length must be NG+1; got len(e_grid)={len(erg)} vs NG={NG}."
        )

    micro_xs = bool(cfg.get("micro_xs", False))
    atom_density = cfg.get("atom_density")
    if micro_xs and atom_density is None:
        raise ValueError(
            "JSON sets 'micro_xs: true' but no 'atom_density' provided."
        )
    if atom_density is not None:
        atom_density = float(atom_density)

    def maybe_macro(vals: List[float], label: str, length: int) -> List[float]:
        base = expand_to_length(vals, length, label)
        if micro_xs and atom_density is not None:
            return [atom_density * x for x in base]
        return base

    sigt = maybe_macro(
        cfg.get("sigma_t")
        or cfg.get("Sigma_t")
        or cfg.get("sigt", []),
        "sigma_t",
        NG,
    )

    sigc = maybe_macro(
        cfg.get("sigma_c")
        or cfg.get("Sigma_c")
        or cfg.get("Sigma_a")
        or cfg.get("sigc", []),
        "sigma_c",
        NG,
    )

    sigf = maybe_macro(
        cfg.get("sigma_f")
        or cfg.get("Sigma_f")
        or cfg.get("sigf", []),
        "sigma_f",
        NG,
    )

    nu = expand_to_length(cfg.get("nu", []), NG, "nu")
    chi = expand_to_length(cfg.get("chi", []), NG, "chi")

    raw_sigs = (
        cfg.get("sigma_s")
        or cfg.get("Sigma_s")
        or cfg.get("sigs")
        or []
    )

    NS = NG * NG
    if raw_sigs and isinstance(raw_sigs[0], list):
        if len(raw_sigs) != NG or any(len(row) != NG for row in raw_sigs):
            raise ValueError("sigma_s matrix must be NG x NG.")
        sigs_flat = flatten_matrix(raw_sigs)
    else:
        sigs_flat = expand_to_length(raw_sigs, NS, "sigma_s")

    sigs = maybe_macro(sigs_flat, "sigma_s", NS)

    raw_sig1 = (
        cfg.get("sigma_s_p1")
        or cfg.get("Sigma_s_p1")
        or cfg.get("sig1")
        or []
    )
    if raw_sig1 and isinstance(raw_sig1[0], list):
        if len(raw_sig1) != NG or any(len(row) != NG for row in raw_sig1):
            raise ValueError("sigma_s_p1 matrix must be NG x NG.")
        sig1_flat = flatten_matrix(raw_sig1)
    else:
        sig1_flat = expand_to_length(raw_sig1, NS, "sigma_s_p1")

    sig1 = maybe_macro(sig1_flat, "sigma_s_p1", NS)

    FISSION = any(sf > 0.0 for sf in sigf)
    if FISSION:
        s_chi = sum(chi)
        if s_chi == 0.0:
            chi = [0.0] * NG
            chi[0] = 1.0
        else:
            chi = [x / s_chi for x in chi]
    else:
        chi = [0.0] * NG
        nu = [0.0] * NG
        sigf = [0.0] * NG

    for arr_name, arr in [
        ("erg", erg),
        ("sigt", sigt),
        ("sigc", sigc),
        ("sigf", sigf),
        ("nu", nu),
        ("chi", chi),
        ("sigs", sigs),
        ("sig1", sig1),
    ]:
        for v in arr:
            if v < 0.0:
                raise ValueError(f"Negative value not allowed in {arr_name}.")

    if all(abs(x) < 1e-12 for x in sigt):
        sigt = []
        for g in range(NG):
            c = sigc[g]
            f = sigf[g]
            s = 0.0
            for j in range(NG):
                # Internal matrix convention: row=to, col=from.
                s += sigs[j * NG + g]
            sigt.append(c + f + s)

    nbins = int(cfg.get("nbins", 1000))

    if "file" in cfg:
        file_out = str(cfg["file"])
    else:
        file_out = f"mg{NG:02d}xsecs"

    mg = MGData(
        zaid=zaid,
        za=za,
        awr=awr,
        tmp=tmp,
        comment=comment,
        NG=NG,
        erg=erg,
        sigt=sigt,
        sigc=sigc,
        sigf=sigf,
        nu=nu,
        chi=chi,
        sigs=sigs,
        sig1=sig1,
        nbins=nbins,
        comment_short=comment_short,
        date=get_today(),
    )

    return mg, file_out


# -----------------------------------------------------------------------------
# Core ACE-building logic
# -----------------------------------------------------------------------------

def build_ace(mg: MGData) -> Tuple[List[int], List[int], List[float]]:
    NG = mg.NG
    erg = mg.erg
    sigt = mg.sigt
    sigc = mg.sigc
    sigf = mg.sigf
    nu = mg.nu
    chi = mg.chi
    sigs = mg.sigs
    sig1 = mg.sig1
    NBINS = mg.nbins

    ISOTROPIC = not any(abs(x) > 0.0 for x in sig1)

    nxs = [0] * 17
    jxs = [0] * 33

    nxs[2] = mg.za
    nxs[3] = 0 if ISOTROPIC else NBINS + 1
    nxs[4] = 0
    nxs[5] = NG
    nxs[6] = NG - 1
    nxs[7] = NG - 1
    nxs[8] = 0
    nxs[9] = 0
    nxs[10] = 1
    nxs[11] = 0
    nxs[12] = 1
    nxs[13] = nxs[14] = nxs[15] = nxs[16] = 0

    jxs[1] = 1
    jxs[2] = jxs[1] + 2 * NG
    jxs[3] = jxs[2] + NG
    jxs[4] = jxs[3] + NG
    jxs[5] = jxs[4] + NG
    jxs[6] = jxs[5] + NG

    for k in range(7, 13):
        jxs[k] = 0

    jxs[13] = jxs[6] + NG * NG + 1
    jxs[14] = jxs[15] = 0
    jxs[16] = jxs[13] + NG * NG + 1

    if ISOTROPIC:
        jxs[17] = jxs[16] + 1
    else:
        jxs[17] = jxs[16] + NG * NG + 1

    for k in range(18, 33):
        jxs[k] = 0

    xss: List[float] = [0.0]

    def put(idx: int, val: float) -> None:
        if idx <= 0:
            raise IndexError("xss is 1-based; idx must be >=1")
        if idx >= len(xss):
            xss.extend([0.0] * (idx - len(xss) + 1))
        xss[idx] = float(val)

    for k in range(NG):
        mid = 0.5 * (erg[k] + erg[k + 1])
        width = erg[k] - erg[k + 1]
        put(jxs[1] + k, mid)
        put(jxs[1] + NG + k, width)

        put(jxs[2] + k, sigt[k])
        put(jxs[3] + k, sigf[k])
        put(jxs[4] + k, nu[k])
        put(jxs[5] + k, chi[k])
        put(jxs[6] + k, sigc[k])

    # Internal matrix convention is row=to, col=from.
    # ACE scatter block is row=from, col=to, so transpose at the boundary.
    put(jxs[13], jxs[13] + 1)
    for i in range(NG):
        for j in range(NG):
            idx = jxs[13] + 1 + i * NG + j
            put(idx, sigs[j * NG + i])

    if ISOTROPIC:
        put(jxs[16], 0.0)
        put(jxs[17], 0.0)
        nxs[1] = len(xss) - 1
    else:
        NEXT16 = jxs[16] + 1
        NEXT17 = jxs[17] + 1
        put(jxs[16], NEXT16)
        put(jxs[17], NEXT17)

        LOFF = 1
        for i in range(NG):
            for j in range(NG):
                k = j * NG + i
                if sigs[k] == 0.0:
                    put(NEXT16, -1.0)
                    NEXT16 += 1
                    continue
                if sig1[k] == 0.0:
                    put(NEXT16, 0.0)
                    NEXT16 += 1
                    continue

                put(NEXT16, LOFF)
                NEXT16 += 1

                LOFF += NBINS + 1

                A0 = 1.0
                A1 = sig1[k] / sigs[k]
                A = [A0, A1]
                mus = get_equiprob_bins(NBINS, 1, A)
                for m in range(NBINS + 1):
                    put(NEXT17, mus[m])
                    NEXT17 += 1

        nxs[1] = len(xss) - 1

    return nxs, jxs, xss


# -----------------------------------------------------------------------------
# ACE reader
# -----------------------------------------------------------------------------

def read_mg_ace(ace_file: str) -> Tuple[MGData, List[int], List[int], List[float]]:
    with open(ace_file, "r") as f:
        header1 = f.readline()
        header2 = f.readline()

        for _ in range(4):
            f.readline()

        nxs_vals: List[int] = []
        while len(nxs_vals) < 16:
            line = f.readline()
            if not line:
                raise ValueError("Unexpected EOF while reading NXS.")
            parts = line.split()
            if not parts:
                continue
            nxs_vals.extend(int(x) for x in parts)
        nxs = [0] + nxs_vals[:16]

        jxs_vals: List[int] = []
        while len(jxs_vals) < 32:
            line = f.readline()
            if not line:
                raise ValueError("Unexpected EOF while reading JXS.")
            parts = line.split()
            if not parts:
                continue
            jxs_vals.extend(int(x) for x in parts)
        jxs = [0] + jxs_vals[:32]

        N = nxs[1]
        xss_vals: List[float] = []
        while len(xss_vals) < N:
            line = f.readline()
            if not line:
                break
            parts = line.split()
            if not parts:
                continue
            xss_vals.extend(float(x) for x in parts)

        if len(xss_vals) < N:
            raise ValueError(f"ACE file ended early: expected {N} XSS values, "
                             f"got {len(xss_vals)}.")

        xss = [0.0] + xss_vals[:N]

    parts1 = header1.split()
    if len(parts1) < 4:
        raise ValueError("ACE header line 1 is malformed.")
    zaid = parts1[0]
    awr = float(parts1[1])
    tmp_mev = float(parts1[2])
    date_str = parts1[3]

    comment_line = header2.rstrip("\n")
    comment_short = comment_line[:70].rstrip()
    comment = comment_short

    NG = nxs[5]
    if NG <= 0:
        raise ValueError(f"Non-positive group count NG={NG} found in ACE file.")

    mids = [xss[jxs[1] + k] for k in range(NG)]
    widths = [xss[jxs[1] + NG + k] for k in range(NG)]

    ehi: List[float] = []
    elow: List[float] = []
    for k in range(NG):
        mid = mids[k]
        w = widths[k]
        ehi.append(mid + 0.5 * w)
        elow.append(mid - 0.5 * w)
    erg: List[float] = [ehi[0]] + elow

    sigt = [xss[jxs[2] + k] for k in range(NG)]
    sigf = [xss[jxs[3] + k] for k in range(NG)]
    nu = [xss[jxs[4] + k] for k in range(NG)]
    chi = [xss[jxs[5] + k] for k in range(NG)]
    sigc = [xss[jxs[6] + k] for k in range(NG)]

    p0_start = int(xss[jxs[13]])
    sigs: List[float] = []
    # Convert ACE layout (row=from, col=to) to internal (row=to, col=from).
    for to_g in range(NG):
        for from_g in range(NG):
            sigs.append(xss[p0_start + from_g * NG + to_g])

    NBINS = nxs[3] - 1 if nxs[3] > 0 else 0
    sig1 = [0.0] * (NG * NG)
    nbins = NBINS

    za = nxs[2]

    mg = MGData(
        zaid=zaid,
        za=za,
        awr=awr,
        tmp=tmp_mev,
        comment=comment,
        NG=NG,
        erg=erg,
        sigt=sigt,
        sigc=sigc,
        sigf=sigf,
        nu=nu,
        chi=chi,
        sigs=sigs,
        sig1=sig1,
        nbins=nbins,
        comment_short=comment_short,
        date=date_str,
    )

    return mg, nxs, jxs, xss


# -----------------------------------------------------------------------------
# HDF5 mgxs.h5 writer (OpenMC-style, isotropic)
# -----------------------------------------------------------------------------

def _set_fixed_ascii_attr(attrs, key, value):
    # OpenSn expects fixed-length C-string attributes in some readers.
    encoded = value.encode("ascii")
    attrs.create(key, encoded, dtype=f"S{len(encoded)}")


def write_mgxs_hdf5_isotropic(
    mg: MGData,
    h5_path: str,
    library_name: Optional[str] = None,
) -> None:
    try:
        import h5py  # type: ignore
    except ImportError as exc:
        raise SystemExit(
            "h5py is required to write HDF5 mgxs files. "
            "Install it with 'pip install h5py'."
        ) from exc

    NG = mg.NG
    if NG <= 0:
        raise ValueError("MGData has non-positive NG.")

    e_desc_mev = mg.erg
    e_asc_eV = [e * MEV_TO_EV for e in reversed(e_desc_mev)]

    kT_eV = mg.tmp * MEV_TO_EV
    if kT_eV > 0.0:
        T_K = kT_eV / KB_EV_PER_K
    else:
        T_K = 294.0
    T_int = int(round(T_K))
    temp_group_name = f"{T_int}K"

    fissionable = any(sf > 0.0 for sf in mg.sigf)

    import h5py  # local
    with h5py.File(h5_path, "w") as f:
        _set_fixed_ascii_attr(f.attrs, "filetype", "mgxs")
        f.attrs["version"] = [1, 0]
        f.attrs["energy_groups"] = NG
        f.attrs["delayed_groups"] = 0
        f.create_dataset("group structure", data=e_asc_eV)

        name = library_name or mg.zaid
        g_lib = f.create_group(name)
        g_lib.attrs["atomic_weight_ratio"] = mg.awr
        g_lib.attrs["fissionable"] = bool(fissionable)
        _set_fixed_ascii_attr(g_lib.attrs, "representation", "isotropic")
        _set_fixed_ascii_attr(g_lib.attrs, "scatter_format", "legendre")
        g_lib.attrs["order"] = 0
        _set_fixed_ascii_attr(g_lib.attrs, "scatter_shape", "[G][G'][Order]")

        g_kTs = g_lib.create_group("kTs")
        g_kTs.create_dataset(temp_group_name, data=kT_eV)

        g_T = g_lib.create_group(temp_group_name)

        total = list(mg.sigt)
        absorption = [c + f_ for c, f_ in zip(mg.sigc, mg.sigf)]
        fission_xs = list(mg.sigf)
        chi = list(mg.chi)
        nu_fission = [nu_ * f_ for nu_, f_ in zip(mg.nu, mg.sigf)]

        g_T.create_dataset("total", data=total)
        g_T.create_dataset("absorption", data=absorption)
        if fissionable:
            g_T.create_dataset("fission", data=fission_xs)
            g_T.create_dataset("chi", data=chi)
            g_T.create_dataset("nu-fission", data=nu_fission)

        g_scat = g_T.create_group("scatter_data")

        g_min = [1] * NG
        g_max = [NG] * NG
        g_scat.create_dataset("g_min", data=g_min)
        g_scat.create_dataset("g_max", data=g_max)

        scat_matrix = []
        for i in range(NG):
            for j in range(NG):
                scat_matrix.append(mg.sigs[i * NG + j])
        g_scat.create_dataset("scatter_matrix", data=scat_matrix)

        multiplicity = [1.0] * len(scat_matrix)
        g_scat.create_dataset("multiplicity_matrix", data=multiplicity)

def write_problem_mgxs_library(
    mg_list: List[MGData],
    h5_path: str,
) -> None:
    """
    Given a list of MGData objects (one per ACE/nuclide) that all share the
    same multigroup structure and temperature, write a single mgxs.h5 file
    containing all of them as separate nuclide groups.

    This is intended for "one HDF5 library per problem directory" where
    that problem has multiple ACE files (multiple isotopes).
    """
    if not mg_list:
        return

    try:
        import h5py  # type: ignore
    except ImportError as exc:
        raise SystemExit(
            "h5py is required to write HDF5 mgxs files. "
            "Install it with 'pip install h5py'."
        ) from exc

    # Use the first isotope as the reference for NG, energy grid, and kT
    ref = mg_list[0]
    NG = ref.NG
    erg_ref = ref.erg
    kT_ref = ref.tmp

    # Sanity: require all isotopes to share NG and energy grid
    from math import isclose

    for mg in mg_list[1:]:
        if mg.NG != NG:
            raise ValueError(
                f"All ACEs combined into one mgxs.h5 must have the same NG; "
                f"found {NG} and {mg.NG}."
            )
        if len(mg.erg) != len(erg_ref):
            raise ValueError("Energy grid lengths differ between ACE files.")
        for a, b in zip(mg.erg, erg_ref):
            if not isclose(a, b, rel_tol=1e-8, abs_tol=1e-12):
                raise ValueError("Energy grids differ between ACE files.")

        if not isclose(mg.tmp, kT_ref, rel_tol=1e-6, abs_tol=1e-12):
            raise ValueError(
                "All ACEs combined into one mgxs.h5 must have the same "
                "temperature (kT in MeV)."
            )

    # Common energy grid: MeV (desc) -> eV (asc)
    e_desc_mev = erg_ref
    e_asc_eV = [e * MEV_TO_EV for e in reversed(e_desc_mev)]

    # Temperature: use the common kT
    kT_eV = kT_ref * MEV_TO_EV
    if kT_eV > 0.0:
        T_K = kT_eV / KB_EV_PER_K
    else:
        T_K = 294.0
    T_int = int(round(T_K))
    temp_group_name = f"{T_int}K"

    with h5py.File(h5_path, "w") as f:
        # File-level attrs (shared)
        _set_fixed_ascii_attr(f.attrs, "filetype", "mgxs")
        f.attrs["version"] = [1, 0]
        f.attrs["energy_groups"] = NG
        f.attrs["delayed_groups"] = 0

        f.create_dataset("group structure", data=e_asc_eV)

        # One nuclide group per MGData
        for mg in mg_list:
            fissionable = any(sf > 0.0 for sf in mg.sigf)
            name = mg.zaid

            g_lib = f.create_group(name)
            g_lib.attrs["atomic_weight_ratio"] = mg.awr
            g_lib.attrs["fissionable"] = bool(fissionable)
            _set_fixed_ascii_attr(g_lib.attrs, "representation", "isotropic")
            _set_fixed_ascii_attr(g_lib.attrs, "scatter_format", "legendre")
            g_lib.attrs["order"] = 0
            _set_fixed_ascii_attr(g_lib.attrs, "scatter_shape", "[G][G'][Order]")

            # Temperature info per nuclide (same kT for all)
            g_kTs = g_lib.create_group("kTs")
            g_kTs.create_dataset(temp_group_name, data=kT_eV)

            g_T = g_lib.create_group(temp_group_name)

            total = list(mg.sigt)
            absorption = [c + f_ for c, f_ in zip(mg.sigc, mg.sigf)]
            fission_xs = list(mg.sigf)
            chi = list(mg.chi)
            nu_fission = [nu_ * f_ for nu_, f_ in zip(mg.nu, mg.sigf)]

            g_T.create_dataset("total", data=total)
            g_T.create_dataset("absorption", data=absorption)
            if fissionable:
                g_T.create_dataset("fission", data=fission_xs)
                g_T.create_dataset("chi", data=chi)
                g_T.create_dataset("nu-fission", data=nu_fission)

            # P0 scatter matrix, [order, g, g'] with order=0
            import numpy as np  # type: ignore
            scat = np.array(mg.sigs, dtype=float).reshape((NG, NG))
            g_scat = g_T.create_group("scatter_data")
            g_scat.create_dataset("scatter", data=scat[None, ...])  # [1, NG, NG]

    print(
        f"[HDF5] Wrote combined mgxs library with {len(mg_list)} nuclides to '{h5_path}'."
    )


def convert_single_ace_to_hdf5(
    ace_file: str,
    h5_out: Optional[str] = None,
    library_name: Optional[str] = None,
) -> None:
    mg, nxs, jxs, xss = read_mg_ace(ace_file)
    if h5_out is None:
        base, _ = os.path.splitext(ace_file)
        h5_out = base + ".h5"

    write_mgxs_hdf5_isotropic(mg, h5_out, library_name=library_name)
    print(f"Converted ACE → mgxs HDF5:\n  ACE: {ace_file}\n  H5 : {h5_out}\n")

# ---------------------------------------------------------------------------
# MCNP geometry -> Gmsh mesh (IN / SL / CY / SP / ISLC-style problems)
# ---------------------------------------------------------------------------

@dataclass
class Plane:
    axis: str   # 'x', 'y', or 'z'
    value: float

@dataclass
class CylinderSurf:
    axis: str   # axis cylinder is parallel to: 'x', 'y', or 'z'
    x0: float
    y0: float
    z0: float
    dx: float
    dy: float
    dz: float
    R: float

@dataclass
class SphereSurf:
    x0: float
    y0: float
    z0: float
    R: float


def _detect_geometry_kind_from_name(path: str) -> str:
    """
    Infer geometry type from problem name suffix, e.g. -IN, -SL, -CY, -SP, -ISLC.
    """
    stem = os.path.splitext(os.path.basename(path))[0].upper()
    for suffix in ("-IN", "-SL", "-CY", "-SP", "-ISLC"):
        if stem.endswith(suffix):
            return suffix[1:]  # IN, SL, CY, SP, ISLC
    return "AUTO"


def _parse_mcnp_surfaces(mcnp_input: str
                         ) -> Tuple[List[Plane], List[Plane], List[Plane],
                                    List[CylinderSurf], List[SphereSurf]]:
    """
    Minimal MCNP surface parser for the benchmark problems. [web:143]

    Supports:
      - PX / PY / PZ planes
      - CZ / CX / CY infinite cylinders
      - RCC finite right circular cylinders
      - SO sphere at origin
      - S / SPH general sphere
    """
    px_list: List[Plane] = []
    py_list: List[Plane] = []
    pz_list: List[Plane] = []
    cyl_list: List[CylinderSurf] = []
    sph_list: List[SphereSurf] = []

    with open(mcnp_input, "r") as f:
        for raw in f:
            line = raw.rstrip("\n")
            # Strip MCNP end-of-line comments. [web:143]
            if "$" in line:
                line = line.split("$", 1)[0]
            if not line.strip():
                continue
            if line[0].lower() == "c":
                continue  # comment line

            tokens = line.split()
            if len(tokens) < 3:
                continue
            # Surface cards start with integer ID
            if not tokens[0].lstrip("*+-").isdigit():
                continue

            mnemonic = tokens[1].lower()

            # Planes
            if mnemonic in ("px", "py", "pz"):
                try:
                    val = float(tokens[2])
                except ValueError:
                    continue
                axis = mnemonic[1]
                plane = Plane(axis=axis, value=val)
                if axis == "x":
                    px_list.append(plane)
                elif axis == "y":
                    py_list.append(plane)
                else:
                    pz_list.append(plane)
                continue

            # Sphere at origin: SO R
            if mnemonic == "so":
                try:
                    R = float(tokens[2])
                except ValueError:
                    continue
                if R > 0.0:
                    sph_list.append(SphereSurf(0.0, 0.0, 0.0, R))
                continue

            # General sphere: S/SPH x0 y0 z0 R
            if mnemonic in ("s", "sph") and len(tokens) >= 6:
                try:
                    x0 = float(tokens[2])
                    y0 = float(tokens[3])
                    z0 = float(tokens[4])
                    R = float(tokens[5])
                except ValueError:
                    continue
                if R > 0.0:
                    sph_list.append(SphereSurf(x0, y0, z0, R))
                continue

            # Infinite cylinders: CZ/CX/CY x0 y0 R
            if mnemonic in ("cz", "cx", "cy") and len(tokens) >= 5:
                try:
                    x0 = float(tokens[2])
                    y0 = float(tokens[3])
                    R = float(tokens[4])
                except ValueError:
                    continue
                if R <= 0.0:
                    continue
                axis = mnemonic[1]
                dx = 1.0 if axis == "x" else 0.0
                dy = 1.0 if axis == "y" else 0.0
                dz = 1.0 if axis == "z" else 0.0
                cyl_list.append(
                    CylinderSurf(axis=axis, x0=x0, y0=y0, z0=0.0,
                                 dx=dx, dy=dy, dz=dz, R=R)
                )
                continue

            # Finite cylinder: RCC x0 y0 z0 dx dy dz R
            if mnemonic == "rcc" and len(tokens) >= 9:
                try:
                    x0 = float(tokens[2])
                    y0 = float(tokens[3])
                    z0 = float(tokens[4])
                    dx = float(tokens[5])
                    dy = float(tokens[6])
                    dz = float(tokens[7])
                    R = float(tokens[8])
                except ValueError:
                    continue
                if R <= 0.0:
                    continue
                axis = "z"
                if abs(dx) >= abs(dy) and abs(dx) >= abs(dz):
                    axis = "x"
                elif abs(dy) >= abs(dx) and abs(dy) >= abs(dz):
                    axis = "y"
                cyl_list.append(
                    CylinderSurf(axis=axis, x0=x0, y0=y0, z0=z0,
                                 dx=dx, dy=dy, dz=dz, R=R)
                )
                continue

    return px_list, py_list, pz_list, cyl_list, sph_list


def _build_infinite_box(gmsh, length: float = 10.0) -> List[Tuple[int, int]]:
    half = 0.5 * length
    tag = gmsh.model.occ.addBox(-half, -half, -half, length, length, length)
    gmsh.model.occ.synchronize()
    return [(3, tag)]


def _build_slab_box(gmsh,
                    px_list: List[Plane],
                    py_list: List[Plane],
                    pz_list: List[Plane]) -> List[Tuple[int, int]]:
    def bounds(planes: List[Plane]) -> Tuple[float, float]:
        vals = sorted(p.value for p in planes)
        if len(vals) < 2:
            raise ValueError
        return vals[0], vals[-1]

    L = 1.0
    try:
        x_min, x_max = bounds(px_list)
        y_min, y_max = -0.5 * L, 0.5 * L
        z_min, z_max = -0.5 * L, 0.5 * L
    except Exception:
        try:
            y_min, y_max = bounds(py_list)
            x_min, x_max = -0.5 * L, 0.5 * L
            z_min, z_max = -0.5 * L, 0.5 * L
        except Exception:
            try:
                z_min, z_max = bounds(pz_list)
                x_min, x_max = -0.5 * L, 0.5 * L
                y_min, y_max = -0.5 * L, 0.5 * L
            except Exception:
                x_min, x_max = -0.5, 0.5
                y_min, y_max = -0.5, 0.5
                z_min, z_max = -0.5, 0.5

    dx = x_max - x_min
    dy = y_max - y_min
    dz = z_max - z_min
    tag = gmsh.model.occ.addBox(x_min, y_min, z_min, dx, dy, dz)
    gmsh.model.occ.synchronize()
    return [(3, tag)]


def _build_sphere_regions(gmsh, spheres: List[SphereSurf]
                          ) -> List[Tuple[int, int]]:
    if not spheres:
        tag = gmsh.model.occ.addSphere(0.0, 0.0, 0.0, 1.0)
        gmsh.model.occ.synchronize()
        return [(3, tag)]

    R_list = sorted({s.R for s in spheres if s.R > 0.0})
    if not R_list:
        tag = gmsh.model.occ.addSphere(0.0, 0.0, 0.0, 1.0)
        gmsh.model.occ.synchronize()
        return [(3, tag)]

    sphere_tags: List[int] = []
    for R in R_list:
        tag = gmsh.model.occ.addSphere(0.0, 0.0, 0.0, R)
        sphere_tags.append(tag)

    gmsh.model.occ.synchronize()

    regions: List[Tuple[int, int]] = [(3, sphere_tags[0])]
    for i in range(1, len(sphere_tags)):
        outer = sphere_tags[i]
        inner = sphere_tags[i - 1]
        cut_vols, _ = gmsh.model.occ.cut(
            [(3, outer)], [(3, inner)],
            removeObject=True, removeTool=False,
        )
        for dim, tag in cut_vols:
            if dim == 3:
                regions.append((dim, tag))

    gmsh.model.occ.synchronize()
    return regions


def _build_cylinder_regions(gmsh,
                            cyl_list: List[CylinderSurf],
                            pz_list: List[Plane]) -> List[Tuple[int, int]]:
    if not cyl_list:
        return _build_slab_box(gmsh, [], [], [])

    c0 = cyl_list[0]
    axis = c0.axis

    if axis == "z" and pz_list:
        vals = sorted(p.value for p in pz_list)
        z_min, z_max = vals[0], vals[-1]
        base = (c0.x0, c0.y0, z_min)
        height = (0.0, 0.0, z_max - z_min)
    else:
        if abs(c0.dx) + abs(c0.dy) + abs(c0.dz) > 0.0:
            base = (c0.x0, c0.y0, c0.z0)
            height = (c0.dx, c0.dy, c0.dz)
        else:
            base = (c0.x0, c0.y0, -1.0)
            height = (0.0, 0.0, 2.0)

    R_list = sorted({c.R for c in cyl_list if c.R > 0.0})
    cyl_tags: List[int] = []
    for R in R_list:
        tag = gmsh.model.occ.addCylinder(
            base[0], base[1], base[2],
            height[0], height[1], height[2],
            R,
        )
        cyl_tags.append(tag)

    gmsh.model.occ.synchronize()

    regions: List[Tuple[int, int]] = [(3, cyl_tags[0])]
    for i in range(1, len(cyl_tags)):
        outer = cyl_tags[i]
        inner = cyl_tags[i - 1]
        cut_vols, _ = gmsh.model.occ.cut(
            [(3, outer)], [(3, inner)],
            removeObject=True, removeTool=False,
        )
        for dim, tag in cut_vols:
            if dim == 3:
                regions.append((dim, tag))

    gmsh.model.occ.synchronize()
    return regions


def mcnp_to_gmsh_mesh(mcnp_input: str,
                      mesh_out: str = "geom.msh") -> None:
    """
    Build a simple 3D volume mesh for an MCNP benchmark input using the
    Gmsh Python API (OpenCASCADE). [web:122]
    """
    try:
        import gmsh  # type: ignore
    except ImportError as exc:
        raise SystemExit(
            "The 'gmsh' Python package is required for MCNP->Gmsh conversion.\n"
            "Install it with e.g. 'python -m pip install gmsh'."
        ) from exc

    geom_kind = _detect_geometry_kind_from_name(mcnp_input)
    px_list, py_list, pz_list, cyl_list, sph_list = _parse_mcnp_surfaces(
        mcnp_input
    )

    gmsh.initialize()
    gmsh.model.add("mcnp_geometry")

    # First, prefer explicit surfaces; only fall back to "infinite box"
    # if no recognizable surfaces are present. [web:143]
    if sph_list:
        regions = _build_sphere_regions(gmsh, sph_list)
    elif cyl_list:
        regions = _build_cylinder_regions(gmsh, cyl_list, pz_list)
    elif px_list or py_list or pz_list or geom_kind in ("SL", "ISLC"):
        # Slab / layer problems, or anything with planes only
        regions = _build_slab_box(gmsh, px_list, py_list, pz_list)
    elif geom_kind == "CY":
        regions = _build_cylinder_regions(gmsh, cyl_list, pz_list)
    elif geom_kind == "SP":
        regions = _build_sphere_regions(gmsh, sph_list)
    elif geom_kind == "IN":
        # Truly infinite problem with no useful surfaces: just give a box
        regions = _build_infinite_box(gmsh)
    else:
        # Last resort
        regions = _build_slab_box(gmsh, px_list, py_list, pz_list)

    for i, (dim, tag) in enumerate(regions, start=1):
        if dim != 3:
            continue
        gmsh.model.addPhysicalGroup(3, [tag], i)
        gmsh.model.setPhysicalName(3, i, f"region_{i}")

    gmsh.model.mesh.generate(3)
    gmsh.write(mesh_out)
    gmsh.finalize()

    print(
        f"[MCNP->Gmsh] Geometry '{geom_kind}' → {len(regions)} region(s), "
        f"mesh written to '{mesh_out}'."
    )


def make_mesh_for_problem(
    src_problem_dir: str,
    out_problem_dir: str,
    gmsh_exe: str,
    msh_name: str = "geom.msh",
) -> None:
    """
    For a given problem directory, try to generate a 3D mesh:

      1) If a .geo or .step/.stp file exists, call the Gmsh CLI on it.
      2) Otherwise, look for an MCNP input (.i/.inp/.mcnp) and use the
         Python gmsh API + MCNP surface cards to mesh the geometry.

    The mesh is stored as geom.msh in out_problem_dir.
    """
    # 1) Prefer explicit Gmsh geometry files
    geo_files = [f for f in os.listdir(src_problem_dir)
                 if f.lower().endswith(".geo")]
    step_files = [f for f in os.listdir(src_problem_dir)
                  if f.lower().endswith(".step") or f.lower().endswith(".stp")]

    geom_src = None
    geom_type = None

    if geo_files:
        geom_src = os.path.join(src_problem_dir, geo_files[0])
        geom_type = "geo"
    elif step_files:
        geom_src = os.path.join(src_problem_dir, step_files[0])
        geom_type = "step"

    os.makedirs(out_problem_dir, exist_ok=True)
    msh_path = os.path.join(out_problem_dir, msh_name)

    if geom_src is not None:
        # Use external Gmsh executable on .geo/.step
        print(f" Geometry: using {geom_type.upper()} file {geom_src}")
        print(f" Geometry: Mesh → {msh_path}")
        try:
            subprocess.check_call(
                [gmsh_exe, geom_src, "-3", "-o", msh_path],
                cwd=src_problem_dir,
            )
        except FileNotFoundError:
            raise SystemExit(
                f"Gmsh executable '{gmsh_exe}' not found. "
                "Install Gmsh or pass --gmsh with the full path."
            )
        print(f" Geometry: created mesh {msh_path}\n")
        return

    # 2) Fallback: MCNP input → Gmsh via Python API
    mcnp_candidates = [
        f for f in os.listdir(src_problem_dir)
        if f.lower().endswith((".i", ".inp", ".mcnp"))
    ]
    if not mcnp_candidates:
        print(
            f" Geometry: no .geo/.step or MCNP input (.i/.inp/.mcnp) found in "
            f"{src_problem_dir}, skipping mesh.\n"
        )
        return

    mcnp_input = os.path.join(src_problem_dir, mcnp_candidates[0])
    print(f" Geometry: using MCNP input {mcnp_input}")
    print(f" Geometry: Mesh → {msh_path}")

    # Use internal MCNP->Gmsh converter
    mcnp_to_gmsh_mesh(mcnp_input, msh_path)
    print(f" Geometry: created mesh {msh_path}\n")


def scan_verification_and_convert(
    root_dir: str,
    out_dir: str,
    library_name: Optional[str] = None,
    make_mesh: bool = False,
    gmsh_exe: str = "gmsh",
) -> None:
    """
    Walk a verification directory tree (e.g., mg01–mg75 under root_dir),
    read every file that looks like a multigroup ACE table, and for each
    *problem directory* write a single mgxs.h5 that contains all isotopes
    in that directory.

    If make_mesh is True, also look for a .geo or .step/.stp geometry file
    in each problem directory and generate geom.msh with Gmsh next to the
    .h5 file in out_dir.
    """
    converted = 0
    scanned = 0
    mesh_done_dirs = set()

    # Accumulate isotopes per problem (output) directory
    per_problem_mg: dict[str, List[MGData]] = {}

    for dirpath, _, filenames in os.walk(root_dir):
        for fname in filenames:
            ace_path = os.path.join(dirpath, fname)
            scanned += 1

            try:
                mg, nxs, jxs, xss = read_mg_ace(ace_path)
            except Exception:
                continue  # not an ACE table in this format

            # Determine where this problem lives under out_dir
            rel = os.path.relpath(ace_path, root_dir)
            rel_dir = os.path.dirname(rel)
            out_problem_dir = os.path.join(out_dir, rel_dir)
            os.makedirs(out_problem_dir, exist_ok=True)

            key = os.path.normpath(out_problem_dir)
            per_problem_mg.setdefault(key, []).append(mg)
            converted += 1

            print(f"Found ACE for problem '{rel_dir or '.'}': {ace_path}")

            # Optional: one mesh per problem directory
            if make_mesh:
                key_dir = os.path.normpath(dirpath)
                if key_dir not in mesh_done_dirs:
                    make_mesh_for_problem(
                        src_problem_dir=dirpath,
                        out_problem_dir=out_problem_dir,
                        gmsh_exe=gmsh_exe,
                    )
                    mesh_done_dirs.add(key_dir)

    # Now write one combined mgxs.h5 per problem directory
    for out_problem_dir, mg_list in per_problem_mg.items():
        h5_path = os.path.join(out_problem_dir, "mgxs.h5")
        write_problem_mgxs_library(mg_list, h5_path)

    print(
        f"Scan complete. Scanned {scanned} files; "
        f"found ACE tables in {len(per_problem_mg)} problem directory(ies), "
        f"and wrote one mgxs.h5 per problem."
    )


# -----------------------------------------------------------------------------
# Summary & ACE writers
# -----------------------------------------------------------------------------

def write_summary(basefile: str, mg: MGData, nxs: List[int],
                  jxs: List[int], xss: List[float]) -> None:
    NG = mg.NG
    erg = mg.erg

    lines: List[str] = []
    lines.append(" group Ehi     Elow      total     capture    scatter"
                 "    fission        nu        chi")
    fmt = (" {g:5d}{Ehi:9.3g}{Elow:9.3g}{tot:12.4g}{cap:12.4g}"
           "{scat:12.4g}{fis:12.4g}{nu:12.4g}{chi:12.4g}")
    for k in range(NG):
        ehi = erg[k]
        elow = erg[k + 1]
        tot = xss[jxs[2] + k]
        fis = xss[jxs[3] + k]
        nuv = xss[jxs[4] + k]
        chi = xss[jxs[5] + k]
        cap = xss[jxs[6] + k]
        scat = 0.0
        for j in range(NG):
            scat += xss[jxs[13] + 1 + k * NG + j]

        lines.append(fmt.format(
            g=k + 1, Ehi=ehi, Elow=elow, tot=tot, cap=cap,
            scat=scat, fis=fis, nu=nuv, chi=chi
        ))

    if NG > 1:
        lines.append("")
        lines.append(" scattering matrix, group-I (down) --> group-J (across)")
        header = " J -->" + "".join(f"{j+1:12d}" for j in range(NG))
        lines.append(header)
        lines.append(" I --v")
        for i in range(NG):
            row = f"{i+1:10d}"
            for j in range(NG):
                val = xss[jxs[13] + 1 + i * NG + j]
                row += f"{val:12g}"
            lines.append(row)

    TOL = 1.0e-5
    for k in range(NG):
        tot = xss[jxs[2] + k]
        fis = xss[jxs[3] + k]
        cap = xss[jxs[6] + k]
        scat = 0.0
        for j in range(NG):
            scat += xss[jxs[13] + 1 + k * NG + j]
        if abs(tot - (cap + fis + scat)) > TOL:
            lines.append(
                f"\n\tGroup {k+1:6d}:\tCap+Fis+Scat differs from Total by >{TOL:g}"
            )

    text = "\n".join(lines) + "\n"
    print(text)

    with open(basefile + ".summary.txt", "w") as f:
        f.write(text)


def write_xsdir(xsdir_file: str, mg: MGData, ace_filename: str,
                nxs: List[int]) -> None:
    za = mg.za
    awr = mg.awr
    tmp = mg.tmp
    nxs1 = nxs[1]

    lines = []
    lines.append("atomic weight ratios")
    lines.append(f"{za:5d} {awr:.6f}")
    lines.append("directory")
    lines.append(
        f"{mg.zaid:>10s} {awr:<12.6g} {ace_filename:<10s} 0 1 1 "
        f"{nxs1:d} 0 0 {tmp:<10.6g}"
    )

    text = "\n".join(lines) + "\n"

    with open(xsdir_file, "w") as f:
        f.write(text)

    print("\n XSDIR Info, written to", xsdir_file, ":\n")
    print(text)


def write_ace(ace_file: str, mg: MGData,
              nxs: List[int], jxs: List[int], xss: List[float]) -> None:
    print(f"\n Creating ACE file: {ace_file}\n")

    with open(ace_file, "w") as f:
        f.write(f"{mg.zaid:>10s} {mg.awr:11.6g} {mg.tmp:11.5e} "
                f"{mg.date:>10s}\n")
        f.write(f"{mg.comment_short:<70s}{mg.date:>10s}\n")

        for _ in range(4):
            f.write(
                f"{0:7d}{0.0:11f}{0:7d}{0.0:11f}"
                f"{0:7d}{0.0:11f}{0:7d}{0.0:11f}\n"
            )

        for k in range(1, 17):
            f.write(f"{nxs[k]:9d}")
            if k % 8 == 0:
                f.write("\n")

        for k in range(1, 33):
            f.write(f"{jxs[k]:9d}")
            if k % 8 == 0:
                f.write("\n")

        n = nxs[1]
        for k in range(1, n + 1):
            f.write(f"{xss[k]:20.12g}")
            if k % 4 == 0:
                f.write("\n")

        rem = n % 4
        if rem == 1:
            f.write(f"{' ':60s}\n")
        elif rem == 2:
            f.write(f"{' ':40s}\n")
        elif rem == 3:
            f.write(f"{' ':20s}\n")


def print_header_info(mg: MGData, nxs: List[int]) -> None:
    print("\n=====> simple_ace_mg_full.py - create simple multigroup ACE file\n")
    print(f"\t zaid = {mg.zaid}")
    print(f"\t za = {mg.za}")
    print(f"\t file = (set below)")
    print(f"\t awr = {mg.awr}")
    print(f"\t tmp = {mg.tmp}")
    print(f"\t groups = {mg.NG}")
    print(f"\t xss size = {nxs[1]}")
    print(f"\t comment = {mg.comment_short}")
    print(f"\t date = {mg.date}\n")


def print_mcnp_snippet(mg: MGData, ace_filename: str, nxs: List[int]) -> None:
    NG = mg.NG
    emax = mg.erg[0]

    print("\n MCNP input snippet for this table:\n")
    print(f"mgopt f {NG}")
    print(
        f"xs1 {mg.zaid} {emax:.5g} {mg.zaid} 0 1 1 {nxs[1]} 0 0 {mg.tmp:.6g}"
    )
    print(f"m1 {mg.zaid} 1.0\n")


# ---------------------------------------------------------------------------
# MCNP geometry → Gmsh 3D mesh for IN / SL / CY / SP / ISLC test problems
# ---------------------------------------------------------------------------

@dataclass
class Plane:
    axis: str   # 'x', 'y', or 'z'
    value: float

@dataclass
class CylinderSurf:
    axis: str   # axis cylinder is parallel to: 'x', 'y', or 'z'
    x0: float
    y0: float
    z0: float
    dx: float
    dy: float
    dz: float
    R: float

@dataclass
class SphereSurf:
    x0: float
    y0: float
    z0: float
    R: float


def _detect_geometry_kind_from_name(path: str) -> str:
    stem = os.path.splitext(os.path.basename(path))[0].upper()
    for suffix in ("-IN", "-SL", "-CY", "-SP", "-ISLC"):
        if stem.endswith(suffix):
            return suffix[1:]  # IN, SL, CY, SP, ISLC
    return "AUTO"


def _parse_mcnp_surfaces(mcnp_input: str
                         ) -> Tuple[List[Plane], List[Plane], List[Plane],
                                    List[CylinderSurf], List[SphereSurf]]:
    px_list: List[Plane] = []
    py_list: List[Plane] = []
    pz_list: List[Plane] = []
    cyl_list: List[CylinderSurf] = []
    sph_list: List[SphereSurf] = []

    with open(mcnp_input, "r") as f:
        for raw in f:
            line = raw.rstrip("\n")
            if "$" in line:
                line = line.split("$", 1)[0]
            if not line.strip():
                continue
            if line[0].lower() == "c":
                continue

            tokens = line.split()
            if len(tokens) < 3:
                continue
            if not tokens[0].lstrip("+-").isdigit():
                continue

            mnemonic = tokens[1].lower()

            # Planes PX/PY/PZ
            if mnemonic in ("px", "py", "pz"):
                try:
                    val = float(tokens[2])
                except ValueError:
                    continue
                axis = mnemonic[1]
                plane = Plane(axis=axis, value=val)
                if axis == "x":
                    px_list.append(plane)
                elif axis == "y":
                    py_list.append(plane)
                else:
                    pz_list.append(plane)
                continue

            # Spheres at origin: SO R
            if mnemonic == "so":
                try:
                    R = float(tokens[2])
                except ValueError:
                    continue
                if R > 0.0:
                    sph_list.append(SphereSurf(0.0, 0.0, 0.0, R))
                continue

            # General spheres: S/SPH x0 y0 z0 R
            if mnemonic in ("s", "sph") and len(tokens) >= 6:
                try:
                    x0 = float(tokens[2])
                    y0 = float(tokens[3])
                    z0 = float(tokens[4])
                    R = float(tokens[5])
                except ValueError:
                    continue
                if R > 0.0:
                    sph_list.append(SphereSurf(x0, y0, z0, R))
                continue

            # Infinite cylinders: CZ/CX/CY
            if mnemonic in ("cz", "cx", "cy") and len(tokens) >= 5:
                try:
                    x0 = float(tokens[2])
                    y0 = float(tokens[3])
                    R = float(tokens[4])
                except ValueError:
                    continue
                if R <= 0.0:
                    continue
                axis = mnemonic[1]
                dx = 1.0 if axis == "x" else 0.0
                dy = 1.0 if axis == "y" else 0.0
                dz = 1.0 if axis == "z" else 0.0
                cyl_list.append(
                    CylinderSurf(axis=axis, x0=x0, y0=y0, z0=0.0,
                                 dx=dx, dy=dy, dz=dz, R=R)
                )
                continue

            # Finite cylinder: RCC x0 y0 z0 dx dy dz R
            if mnemonic == "rcc" and len(tokens) >= 9:
                try:
                    x0 = float(tokens[2])
                    y0 = float(tokens[3])
                    z0 = float(tokens[4])
                    dx = float(tokens[5])
                    dy = float(tokens[6])
                    dz = float(tokens[7])
                    R = float(tokens[8])
                except ValueError:
                    continue
                if R <= 0.0:
                    continue
                axis = "z"
                if abs(dx) >= abs(dy) and abs(dx) >= abs(dz):
                    axis = "x"
                elif abs(dy) >= abs(dx) and abs(dy) >= abs(dz):
                    axis = "y"
                cyl_list.append(
                    CylinderSurf(axis=axis, x0=x0, y0=y0, z0=z0,
                                 dx=dx, dy=dy, dz=dz, R=R)
                )
                continue

    return px_list, py_list, pz_list, cyl_list, sph_list


def _build_infinite_box(gmsh, length: float = 10.0) -> List[Tuple[int, int]]:
    half = 0.5 * length
    tag = gmsh.model.occ.addBox(-half, -half, -half, length, length, length)
    gmsh.model.occ.synchronize()
    return [(3, tag)]


def _build_slab_box(gmsh,
                    px_list: List[Plane],
                    py_list: List[Plane],
                    pz_list: List[Plane]) -> List[Tuple[int, int]]:
    def bounds(planes: List[Plane]) -> Tuple[float, float]:
        vals = sorted(p.value for p in planes)
        if len(vals) < 2:
            raise ValueError
        return vals[0], vals[-1]

    L_other = 1.0
    try:
        x_min, x_max = bounds(px_list)
        y_min, y_max = -0.5 * L_other, 0.5 * L_other
        z_min, z_max = -0.5 * L_other, 0.5 * L_other
    except Exception:
        try:
            y_min, y_max = bounds(py_list)
            x_min, x_max = -0.5 * L_other, 0.5 * L_other
            z_min, z_max = -0.5 * L_other, 0.5 * L_other
        except Exception:
            try:
                z_min, z_max = bounds(pz_list)
                x_min, x_max = -0.5 * L_other, 0.5 * L_other
                y_min, y_max = -0.5 * L_other, 0.5 * L_other
            except Exception:
                x_min, x_max = -0.5, 0.5
                y_min, y_max = -0.5, 0.5
                z_min, z_max = -0.5, 0.5

    dx = x_max - x_min
    dy = y_max - y_min
    dz = z_max - z_min
    tag = gmsh.model.occ.addBox(x_min, y_min, z_min, dx, dy, dz)
    gmsh.model.occ.synchronize()
    return [(3, tag)]


def _build_sphere_regions(gmsh, spheres: List[SphereSurf]
                          ) -> List[Tuple[int, int]]:
    if not spheres:
        tag = gmsh.model.occ.addSphere(0.0, 0.0, 0.0, 1.0)
        gmsh.model.occ.synchronize()
        return [(3, tag)]

    R_list = sorted({s.R for s in spheres if s.R > 0.0})
    if not R_list:
        tag = gmsh.model.occ.addSphere(0.0, 0.0, 0.0, 1.0)
        gmsh.model.occ.synchronize()
        return [(3, tag)]

    sphere_tags: List[int] = []
    for R in R_list:
        tag = gmsh.model.occ.addSphere(0.0, 0.0, 0.0, R)
        sphere_tags.append(tag)

    gmsh.model.occ.synchronize()

    regions: List[Tuple[int, int]] = [(3, sphere_tags[0])]
    for i in range(1, len(sphere_tags)):
        outer = sphere_tags[i]
        inner = sphere_tags[i - 1]
        cut_vols, _ = gmsh.model.occ.cut(
            [(3, outer)], [(3, inner)],
            removeObject=True, removeTool=False,
        )
        for dim, tag in cut_vols:
            if dim == 3:
                regions.append((dim, tag))

    gmsh.model.occ.synchronize()
    return regions


def _build_cylinder_regions(gmsh,
                            cyl_list: List[CylinderSurf],
                            pz_list: List[Plane]) -> List[Tuple[int, int]]:
    if not cyl_list:
        return _build_slab_box(gmsh, [], [], [])

    c0 = cyl_list[0]
    axis = c0.axis

    if axis == "z" and pz_list:
        vals = sorted(p.value for p in pz_list)
        z_min, z_max = vals[0], vals[-1]
        base = (c0.x0, c0.y0, z_min)
        height = (0.0, 0.0, z_max - z_min)
    else:
        if abs(c0.dx) + abs(c0.dy) + abs(c0.dz) > 0.0:
            base = (c0.x0, c0.y0, c0.z0)
            height = (c0.dx, c0.dy, c0.dz)
        else:
            base = (c0.x0, c0.y0, -1.0)
            height = (0.0, 0.0, 2.0)

    R_list = sorted({c.R for c in cyl_list if c.R > 0.0})
    cyl_tags: List[int] = []
    for R in R_list:
        tag = gmsh.model.occ.addCylinder(
            base[0], base[1], base[2],
            height[0], height[1], height[2],
            R,
        )
        cyl_tags.append(tag)

    gmsh.model.occ.synchronize()

    regions: List[Tuple[int, int]] = [(3, cyl_tags[0])]
    for i in range(1, len(cyl_tags)):
        outer = cyl_tags[i]
        inner = cyl_tags[i - 1]
        cut_vols, _ = gmsh.model.occ.cut(
            [(3, outer)], [(3, inner)],
            removeObject=True, removeTool=False,
        )
        for dim, tag in cut_vols:
            if dim == 3:
                regions.append((dim, tag))

    gmsh.model.occ.synchronize()
    return regions


def mcnp_to_gmsh_mesh(mcnp_input: str,
                      mesh_out: str = "geom.msh") -> None:
    try:
        import gmsh  # type: ignore
    except ImportError as exc:
        raise SystemExit(
            "The 'gmsh' Python package is required for MCNP->Gmsh conversion.\n"
            "Install it with e.g. 'pip install gmsh'."
        ) from exc

    geom_kind = _detect_geometry_kind_from_name(mcnp_input)
    px_list, py_list, pz_list, cyl_list, sph_list = _parse_mcnp_surfaces(
        mcnp_input
    )

    gmsh.initialize()
    gmsh.model.add("mcnp_geometry")

    if geom_kind == "IN":
        regions = _build_infinite_box(gmsh)
    elif geom_kind in ("SL", "ISLC"):
        regions = _build_slab_box(gmsh, px_list, py_list, pz_list)
    elif geom_kind == "CY":
        regions = _build_cylinder_regions(gmsh, cyl_list, pz_list)
    elif geom_kind == "SP":
        regions = _build_sphere_regions(gmsh, sph_list)
    else:
        if sph_list:
            regions = _build_sphere_regions(gmsh, sph_list)
        elif cyl_list:
            regions = _build_cylinder_regions(gmsh, cyl_list, pz_list)
        else:
            regions = _build_slab_box(gmsh, px_list, py_list, pz_list)

    for i, (dim, tag) in enumerate(regions, start=1):
        if dim != 3:
            continue
        gmsh.model.addPhysicalGroup(3, [tag], i)
        gmsh.model.setPhysicalName(3, i, f"region_{i}")

    gmsh.model.mesh.generate(3)
    gmsh.write(mesh_out)
    gmsh.finalize()

    print(f"[MCNP->Gmsh] Geometry '{geom_kind}' → {len(regions)} region(s), mesh written to '{mesh_out}'.")


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Generate multigroup ACE from JSON, or convert ACE to "
            "OpenMC-style mgxs HDF5, optionally also generating "
            "Gmsh meshes."
        )
    )

    # JSON → ACE (original workflow)
    parser.add_argument(
        "--config", "-config", "-j",
        help=(
            "JSON configuration file describing the multigroup XS "
            "(ACE-generation workflow)."
        ),
    )

    # Single ACE → HDF5
    parser.add_argument(
        "--ace",
        help="Path to a multigroup ACE file to convert to mgxs HDF5.",
    )
    parser.add_argument(
        "--h5",
        help=(
            "Output HDF5 filename for ACE→HDF5 conversion. "
            "If omitted, uses the ACE filename with .h5 extension."
        ),
    )

    # Tree scan: many ACE → many HDF5 (+ optional meshes)
    parser.add_argument(
        "--verification-root",
        help=(
            "Root directory to scan recursively for ACE files "
            "(e.g., C:\\Users\\lukes\\MCNP\\utils\\vnvstats\\verification\\keff\\problems). "
            "All valid multigroup ACE tables are converted automatically "
            "to mgxs HDF5 under --outdir."
        ),
    )
    parser.add_argument(
        "--outdir",
        default="h5_out",
        help=(
            "Output directory for HDF5 files when using --verification-root. "
            "Defaults to ./h5_out."
        ),
    )
    parser.add_argument(
        "--library-name",
        default=None,
        help=(
            "Optional library/group name to use inside mgxs.h5. "
            "If omitted, the ZAID from the ACE/JSON is used."
        ),
    )
    parser.add_argument(
        "--make-mesh",
        action="store_true",
        help=(
            "When used with --verification-root, also create a Gmsh mesh "
            "for each problem based on an existing .geo or .step/.stp "
            "geometry file, placing geom.msh next to the .h5 file in "
            "--outdir."
        ),
    )
    parser.add_argument(
        "--gmsh",
        default="gmsh",
        help="Path to the Gmsh executable (default: 'gmsh').",
    )

    args = parser.parse_args()

    # 1) Tree scan mode: many ACE → many HDF5 (+ optional meshes)
    if args.verification_root:
        scan_verification_and_convert(
            root_dir=args.verification_root,
            out_dir=args.outdir,
            library_name=args.library_name,
            make_mesh=args.make_mesh,
            gmsh_exe=args.gmsh,
        )
        return

    # 2) Single ACE → single HDF5
    if args.ace:
        convert_single_ace_to_hdf5(
            ace_file=args.ace,
            h5_out=args.h5,
            library_name=args.library_name,
        )
        return

    # 3) JSON → ACE (original/simple_ace_mg.pl behavior)
    if not args.config:
        raise SystemExit(
            "Please provide either --config (JSON→ACE), "
            "--ace (ACE→HDF5), or --verification-root (tree scan)."
        )

    mg, ace_file = load_from_json(args.config)
    nxs, jxs, xss = build_ace(mg)

    print_header_info(mg, nxs)
    write_summary(ace_file, mg, nxs, jxs, xss)

    xsdir_file = "xsdir_mcnp6.3"
    write_xsdir(xsdir_file, mg, ace_file, nxs)
    write_ace(ace_file, mg, nxs, jxs, xss)
    print_mcnp_snippet(mg, ace_file, nxs)


if __name__ == "__main__":
    # Keep the historic script entry point on the maintained package path.
    from multi_group_ace.cli import main as package_main
    package_main()

import json
from datetime import date

from .core import MGData


def get_today():
    d = date.today()
    return f"{d.year:04d}-{d.month:02d}-{d.day:02d}"


def flatten_matrix(mat):
    flat = []
    for row in mat:
        flat.extend(row)
    return flat


def expand_to_length(vals, length, label):
    if not vals:
        return [0.0] * length
    if len(vals) == 1:
        return [float(vals[0])] * length
    if len(vals) != length:
        raise ValueError(f"{label}: expected 0, 1, or {length} values, got {len(vals)}")
    return [float(x) for x in vals]


def ensure_descending(energies):
    if len(energies) < 2:
        raise ValueError("Energy grid must have at least 2 values (NG+1).")
    if energies[0] < energies[-1]:
        energies = list(reversed(energies))
    for i in range(len(energies) - 1):
        if energies[i] <= energies[i + 1]:
            raise ValueError("Energy grid must be strictly decreasing.")
    return [float(e) for e in energies]


def _scatter_to_flat(raw, NG, label):
    if not raw:
        return [0.0] * (NG * NG)
    if isinstance(raw[0], list):
        if len(raw) != NG or any(len(row) != NG for row in raw):
            raise ValueError(f"{label} matrix must be NG x NG.")
        # Matrix convention is row=to group, col=from group:
        # [ [1->1, 2->1, ...], [1->2, 2->2, ...], ... ]
        return flatten_matrix(raw)
    return expand_to_length(raw, NG * NG, label)


def load_from_json(path):
    with open(path, "r") as f:
        cfg = json.load(f)

    
    required = ["zaid", "NG", "erg", "sigt", "sigc", "sigf", "nu", "chi", "sigs"]
    for key in required:
        if key not in cfg:
            raise KeyError(f"Missing required key in JSON: {key}")

    zaid = str(cfg["zaid"])
    if "." not in zaid:
        raise ValueError("zaid must look like '99902.01m'.")
    za = int(zaid.split(".")[0])

    NG = int(cfg["NG"])
    erg = ensure_descending([float(x) for x in cfg["erg"]])
    if len(erg) != NG + 1:
        raise ValueError(f"Energy grid length must be NG+1, got {len(erg)} for NG={NG}.")

    awr = float(cfg.get("awr", 1.0))
    tmp = float(cfg.get("tmp", 2.5301e-8))
    if tmp > 1.0:
        tmp *= 8.617385e-11

    comment = str(cfg.get("comment", "multigroup ACE file"))
    comment_short = comment[:70]

    is_micro = bool(cfg.get("is_micro", False))
    atom_density = cfg.get("atom_density")
    if is_micro and atom_density is None:
        raise ValueError("JSON sets is_micro=true but atom_density is missing.")
    if atom_density is not None:
        atom_density = float(atom_density)

    def maybe_macro(raw, label, length):
        vals = expand_to_length(raw, length, label)
        if is_micro:
            return [atom_density * x for x in vals]
        return vals

    sigt = maybe_macro(cfg["sigt"], "sigt", NG)
    sigc = maybe_macro(cfg["sigc"], "sigc", NG)
    sigf = maybe_macro(cfg["sigf"], "sigf", NG)
    nu = expand_to_length(cfg["nu"], NG, "nu")
    chi = expand_to_length(cfg["chi"], NG, "chi")
    sigs = maybe_macro(_scatter_to_flat(cfg["sigs"], NG, "sigs"), "sigs", NG * NG)

    sig1_raw = cfg.get("sig1", [])
    sig1 = maybe_macro(_scatter_to_flat(sig1_raw, NG, "sig1"), "sig1", NG * NG) if sig1_raw else [0.0] * (NG * NG)

    if any(sf > 0.0 for sf in sigf):
        s_chi = sum(chi)
        if s_chi == 0.0:
            chi = [0.0] * NG
            chi[0] = 1.0
        else:
            chi = [x / s_chi for x in chi]
    else:
        sigf = [0.0] * NG
        nu = [0.0] * NG
        chi = [0.0] * NG

    nbins = int(cfg.get("nbins", 1000))
    file_out = str(cfg.get("file", f"mg{NG:02d}xsecs"))

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

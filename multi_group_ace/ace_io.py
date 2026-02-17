from .core import MGData


def read_mg_ace(ace_file):
    with open(ace_file, "r") as f:
        header1 = f.readline()
        header2 = f.readline()

        for _ in range(4):
            f.readline()

        nxs_vals = []
        while len(nxs_vals) < 16:
            nxs_vals.extend(int(x) for x in f.readline().split())
        nxs = [0] + nxs_vals[:16]

        jxs_vals = []
        while len(jxs_vals) < 32:
            jxs_vals.extend(int(x) for x in f.readline().split())
        jxs = [0] + jxs_vals[:32]

        N = nxs[1]
        xss_vals = []
        while len(xss_vals) < N:
            line = f.readline()
            if not line:
                break
            parts = line.split()
            if parts:
                xss_vals.extend(float(x) for x in parts)
        if len(xss_vals) < N:
            raise ValueError("ACE file ended before all XSS values were read.")
        xss = [0.0] + xss_vals[:N]

    parts1 = header1.split()
    zaid = parts1[0]
    awr = float(parts1[1])
    tmp_mev = float(parts1[2])
    date_str = parts1[3]

    comment_short = header2.rstrip("\n")[:70].rstrip()
    comment = comment_short

    NG = nxs[5]
    mids = [xss[jxs[1] + k] for k in range(NG)]
    widths = [xss[jxs[1] + NG + k] for k in range(NG)]

    ehi = []
    elow = []
    for k in range(NG):
        ehi.append(mids[k] + 0.5 * widths[k])
        elow.append(mids[k] - 0.5 * widths[k])
    erg = [ehi[0]] + elow

    sigt = [xss[jxs[2] + k] for k in range(NG)]
    sigf = [xss[jxs[3] + k] for k in range(NG)]
    nu = [xss[jxs[4] + k] for k in range(NG)]
    chi = [xss[jxs[5] + k] for k in range(NG)]
    sigc = [xss[jxs[6] + k] for k in range(NG)]

    p0_start = int(xss[jxs[13]])
    sigs = []
    # Convert ACE scatter layout (row=from, col=to) to internal layout
    # (row=to, col=from).
    for to_g in range(NG):
        for from_g in range(NG):
            sigs.append(xss[p0_start + from_g * NG + to_g])

    mg = MGData(
        zaid=zaid,
        za=nxs[2],
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
        sig1=[0.0] * (NG * NG),
        nbins=max(nxs[3] - 1, 0),
        comment_short=comment_short,
        date=date_str,
    )

    return mg, nxs, jxs, xss


def write_ace(ace_file, mg, nxs, jxs, xss):
    with open(ace_file, "w") as f:
        f.write(f"{mg.zaid:>10s} {mg.awr:11.6g} {mg.tmp:11.5e} {mg.date:>10s}\n")
        f.write(f"{mg.comment_short:<70s}{mg.date:>10s}\n")

        for _ in range(4):
            f.write(f"{0:7d}{0.0:11f}{0:7d}{0.0:11f}{0:7d}{0.0:11f}{0:7d}{0.0:11f}\n")

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


def write_xsdir(xsdir_file, mg, ace_filename, nxs):
    text = "\n".join([
        "atomic weight ratios",
        f"{mg.za:5d} {mg.awr:.6f}",
        "directory",
        f"{mg.zaid:>10s} {mg.awr:<12.6g} {ace_filename:<10s} 0 1 1 {nxs[1]} 0 0 {mg.tmp:<10.6g}",
    ]) + "\n"

    with open(xsdir_file, "w") as f:
        f.write(text)


def write_summary(basefile, mg, nxs, jxs, xss):
    NG = mg.NG
    lines = []
    lines.append(" group Ehi     Elow      total     capture    scatter    fission        nu        chi")

    for k in range(NG):
        scat = 0.0
        # xss stores ACE scatter layout (row=from, col=to).
        for j in range(NG):
            scat += xss[jxs[13] + 1 + k * NG + j]
        lines.append(
            f"{k+1:6d}"
            f"{mg.erg[k]:9.3g}{mg.erg[k+1]:9.3g}"
            f"{xss[jxs[2]+k]:12.4g}{xss[jxs[6]+k]:12.4g}{scat:12.4g}"
            f"{xss[jxs[3]+k]:12.4g}{xss[jxs[4]+k]:12.4g}{xss[jxs[5]+k]:12.4g}"
        )

    with open(basefile + ".summary.txt", "w") as f:
        f.write("\n".join(lines) + "\n")

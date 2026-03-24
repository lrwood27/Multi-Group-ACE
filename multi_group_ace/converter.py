def get_equiprob_bins(nbins, N, A):
    nx = 1000
    dx = 2.0 / nx

    x_vals = []
    cdf = []

    for k in range(nx + 1):
        x = -1.0 + k * dx
        x_vals.append(x)

        val = 1.0 + x
        if N > 0:
            P = [0.0] * (N + 2)
            P[0] = 1.0
            P[1] = x
            for j in range(2, N + 2):
                P[j] = ((2.0 * j - 1.0) * x * P[j - 1] - (j - 1.0) * P[j - 2]) / j
            for j in range(1, N + 1):
                val += A[j] * (P[j + 1] - P[j - 1])
        val *= 0.5
        cdf.append(val)

    prob = 1.0 / nbins
    mu = [-1.0] * (nbins + 1)

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


def build_ace(mg):
    # Convert the package's internal data model into the three core ACE arrays:
    #   NXS = global table metadata
    #   JXS = pointers into the XSS payload
    #   XSS = the floating-point data blocks themselves
    NG = mg.NG
    NBINS = mg.nbins
    ISOTROPIC = not any(abs(x) > 0.0 for x in mg.sig1)

    nxs = [0] * 17
    jxs = [0] * 33

    # Fill the control arrays first. These values tell ACE readers what kind
    # of table this is and where each data block starts.
    nxs[2] = mg.za
    nxs[3] = 0 if ISOTROPIC else NBINS + 1
    nxs[5] = NG
    nxs[6] = NG - 1
    nxs[7] = NG - 1
    nxs[10] = 1
    nxs[12] = 1

    jxs[1] = 1
    jxs[2] = jxs[1] + 2 * NG
    jxs[3] = jxs[2] + NG
    jxs[4] = jxs[3] + NG
    jxs[5] = jxs[4] + NG
    jxs[6] = jxs[5] + NG
    jxs[13] = jxs[6] + NG * NG + 1
    jxs[16] = jxs[13] + NG * NG + 1
    jxs[17] = jxs[16] + 1 if ISOTROPIC else jxs[16] + NG * NG + 1

    xss = [0.0]

    def put(idx, val):
        if idx >= len(xss):
            xss.extend([0.0] * (idx - len(xss) + 1))
        xss[idx] = float(val)

    # ACE stores the energy grid as midpoint/width pairs rather than directly
    # storing the group boundaries.
    for k in range(NG):
        mid = 0.5 * (mg.erg[k] + mg.erg[k + 1])
        width = mg.erg[k] - mg.erg[k + 1]
        put(jxs[1] + k, mid)
        put(jxs[1] + NG + k, width)

        put(jxs[2] + k, mg.sigt[k])
        put(jxs[3] + k, mg.sigf[k])
        put(jxs[4] + k, mg.nu[k])
        put(jxs[5] + k, mg.chi[k])
        put(jxs[6] + k, mg.sigc[k])

    # Internal convention:
    #   mg.sigs[row, col] == Sigma_s(from=col -> to=row)
    # ACE block convention:
    #   row = from, col = to
    # so we transpose while writing ACE.
    # The first value at JXS(13) is itself a pointer to the start of the
    # P0 scatter matrix block.
    put(jxs[13], jxs[13] + 1)
    for from_g in range(NG):
        for to_g in range(NG):
            ace_idx = jxs[13] + 1 + from_g * NG + to_g
            mg_idx = to_g * NG + from_g
            put(ace_idx, mg.sigs[mg_idx])

    if ISOTROPIC:
        # A zeroed angular section tells the ACE reader there is no higher-
        # order angular information to recover.
        put(jxs[16], 0.0)
        put(jxs[17], 0.0)
    else:
        # For anisotropic data, write:
        #   JXS(16) -> location table for each scattering pair
        #   JXS(17) -> the actual mu-bin edges
        next16 = jxs[16] + 1
        next17 = jxs[17] + 1
        put(jxs[16], next16)
        put(jxs[17], next17)

        loff = 1
        for from_g in range(NG):
            for to_g in range(NG):
                k = to_g * NG + from_g
                if mg.sigs[k] == 0.0:
                    put(next16, -1.0)
                    next16 += 1
                    continue
                if mg.sig1[k] == 0.0:
                    put(next16, 0.0)
                    next16 += 1
                    continue

                put(next16, loff)
                next16 += 1
                loff += NBINS + 1

                A = [1.0, mg.sig1[k] / mg.sigs[k]]
                mus = get_equiprob_bins(NBINS, 1, A)
                for m in range(NBINS + 1):
                    put(next17, mus[m])
                    next17 += 1

    # NXS(1) is the final XSS payload length.
    nxs[1] = len(xss) - 1
    return nxs, jxs, xss

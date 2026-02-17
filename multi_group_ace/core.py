KB_EV_PER_K = 8.617333262e-5
MEV_TO_EV = 1.0e6


class MGData:
    # Plain container for multigroup data.
    # sigs/sig1 store NGxNG matrices in row-major order with:
    # row = destination group, column = source group.
    # Example for NG=2: [1->1, 2->1, 1->2, 2->2]
    def __init__(
        self,
        zaid,
        za,
        awr,
        tmp,
        comment,
        NG,
        erg,
        sigt,
        sigc,
        sigf,
        nu,
        chi,
        sigs,
        sig1=None,
        nbins=1000,
        comment_short="",
        date="",
    ):
        self.zaid = zaid
        self.za = za
        self.awr = awr
        self.tmp = tmp
        self.comment = comment
        self.NG = NG
        self.erg = erg
        self.sigt = sigt
        self.sigc = sigc
        self.sigf = sigf
        self.nu = nu
        self.chi = chi
        self.sigs = sigs
        self.sig1 = sig1 if sig1 is not None else [0.0] * (NG * NG)
        self.nbins = nbins
        self.comment_short = comment_short
        self.date = date

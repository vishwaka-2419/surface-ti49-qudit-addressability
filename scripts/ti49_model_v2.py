import numpy as np

MU_B_MHZ_T = 13996.24555
MU_N_MHZ_T = 7.622593285

I = 3.5


def spin(j):
    m = np.arange(j, -j - 1, -1, dtype=float)

    Jz = np.diag(m)
    Jp = np.zeros((len(m), len(m)), dtype=complex)

    for c, mc in enumerate(m):
        q = np.where(np.isclose(m, mc + 1))[0]

        if len(q):
            Jp[q[0], c] = np.sqrt(
                j * (j + 1) - mc * (mc + 1)
            )

    Jm = Jp.conj().T

    return (
        (Jp + Jm) / 2,
        (Jp - Jm) / (2j),
        Jz,
    )


sx, sy, sz = spin(0.5)
ix, iy, iz = spin(I)

Sx = np.kron(sx, np.eye(8))
Sy = np.kron(sy, np.eye(8))
Sz = np.kron(sz, np.eye(8))

Ix = np.kron(np.eye(2), ix)
Iy = np.kron(np.eye(2), iy)
Iz = np.kron(np.eye(2), iz)


def hamiltonian(
    B_T,
    Az_MHz=130.0,
    kappa_MHz=-46.373841059602654,
    gN_mag=0.3146857142857143,
    ge_z=0.56,
    Btip_T=0.0679,
    phi_deg=5.0,
    Aperp_over_Az=0.0,
    ge_perp_over_ge_z=1.0,
):
    """
    16-level effective S=1/2, I=7/2 Hamiltonian.

    Important conventions:
    - Az=130 MHz is the nominal 49Ti input.
    - kappa49 is isotope-scaled, not directly measured.
    - ge_z, Btip and phi come from the 47Ti same-geometry
      experiment and are contextual/sensitivity inputs.
    - Aperp is a sensitivity variable, not a measured 49Ti value.
    - ge_perp/ge_z is a sensitivity variable because the complete
      49Ti electron g tensor is not established here.
    - fitted tip field is excluded from the nuclear Zeeman term.
    """

    phi = np.deg2rad(phi_deg)

    Aperp = Aperp_over_Az * Az_MHz
    ge_perp = ge_perp_over_ge_z * ge_z

    qpref = (
        kappa_MHz
        / (2 * I * (2 * I - 1))
    )

    HQ = qpref * (
        -0.5 * Ix @ Ix
        -0.5 * Iy @ Iy
        + Iz @ Iz
    )

    H = (
        MU_B_MHZ_T * ge_z * B_T * Sz
        + MU_B_MHZ_T * ge_z
          * Btip_T * np.cos(phi) * Sz
        + MU_B_MHZ_T * ge_perp
          * Btip_T * np.sin(phi) * Sx
        + MU_N_MHZ_T * gN_mag * B_T * Iz
        + Az_MHz * (Sz @ Iz)
        + Aperp * (
            Sx @ Ix
            + Sy @ Iy
        )
        + HQ
    )

    return H


def nuclear_transitions(
    B_T,
    **kwargs,
):
    H = hamiltonian(
        B_T,
        **kwargs,
    )

    E, V = np.linalg.eigh(H)

    D = V.conj().T @ Ix @ V
    W = np.abs(D) ** 2

    Sz_e = np.real(
        np.diag(
            V.conj().T @ Sz @ V
        )
    )

    records = []

    for i in range(15):
        for j in range(i + 1, 16):

            f = float(E[j] - E[i])

            if not (20.0 < f < 110.0):
                continue

            si = float(Sz_e[i])
            sj = float(Sz_e[j])

            if si < 0 and sj < 0:
                branch = "ground"

            elif si > 0 and sj > 0:
                branch = "excited"

            else:
                branch = "cross"

            records.append(
                {
                    "i": i,
                    "j": j,
                    "freq_MHz": f,
                    "Ix_strength":
                        float(W[i, j]),
                    "Sz_i": si,
                    "Sz_j": sj,
                    "branch": branch,
                }
            )

    return E, V, records


def strongest_nuclear_lines(
    B_T,
    **kwargs,
):
    _, _, records = nuclear_transitions(
        B_T,
        **kwargs,
    )

    out = {}

    for branch in (
        "ground",
        "excited",
    ):
        r = [
            x for x in records
            if x["branch"] == branch
        ]

        r = sorted(
            r,
            key=lambda x:
                x["Ix_strength"],
            reverse=True,
        )[:7]

        if len(r) != 7:
            raise RuntimeError(
                f"Could not identify seven "
                f"{branch} nuclear transitions "
                f"at B={B_T}"
            )

        out[branch] = sorted(
            r,
            key=lambda x:
                x["freq_MHz"],
        )

    return out


def addressability(
    B_T,
    **kwargs,
):
    lines = strongest_nuclear_lines(
        B_T,
        **kwargs,
    )

    g = lines["ground"]
    e = lines["excited"]

    fg = np.array(
        [x["freq_MHz"] for x in g]
    )

    fe = np.array(
        [x["freq_MHz"] for x in e]
    )

    D = np.abs(
        fg[:, None]
        - fe[None, :]
    )

    ig, ie = np.unravel_index(
        np.argmin(D),
        D.shape,
    )

    intra_g = np.min(
        np.diff(
            np.sort(fg)
        )
    )

    intra_e = np.min(
        np.diff(
            np.sort(fe)
        )
    )

    purity = min(
        min(
            2 * abs(x["Sz_i"]),
            2 * abs(x["Sz_j"]),
        )
        for x in g + e
    )

    return {
        "min_inter_MHz":
            float(D[ig, ie]),

        "g_freq_MHz":
            float(fg[ig]),

        "e_freq_MHz":
            float(fe[ie]),

        "g_strength":
            float(g[ig]["Ix_strength"]),

        "e_strength":
            float(e[ie]["Ix_strength"]),

        "min_intra_g_MHz":
            float(intra_g),

        "min_intra_e_MHz":
            float(intra_e),

        "min_electron_purity":
            float(purity),

        "g_i": int(g[ig]["i"]),
        "g_j": int(g[ig]["j"]),
        "e_i": int(e[ie]["i"]),
        "e_j": int(e[ie]["j"]),
    }

def branch_basis(
    B_T,
    branch,
    **kwargs,
):
    """
    Return the eight eigenstates belonging predominantly
    to one electron-spin branch, ordered by <Iz>.

    branch = 'ground'  -> <Sz> < 0
    branch = 'excited' -> <Sz> > 0
    """

    H = hamiltonian(
        B_T,
        **kwargs,
    )

    E, V = np.linalg.eigh(H)

    Sz_e = np.real(
        np.diag(
            V.conj().T @ Sz @ V
        )
    )

    Iz_e = np.real(
        np.diag(
            V.conj().T @ Iz @ V
        )
    )

    if branch == "ground":
        ids = np.where(
            Sz_e < 0
        )[0]

    elif branch == "excited":
        ids = np.where(
            Sz_e > 0
        )[0]

    else:
        raise ValueError(
            "branch must be "
            "'ground' or 'excited'"
        )

    if len(ids) != 8:
        raise RuntimeError(
            f"Expected 8 states in {branch}, "
            f"found {len(ids)}"
        )

    # Order from mI ~ +7/2 to -7/2.
    ids = ids[
        np.argsort(
            Iz_e[ids]
        )[::-1]
    ]

    Vb = V[:, ids]

    D = (
        Vb.conj().T
        @ Ix
        @ Vb
    )

    freqs = np.abs(
        np.diff(
            E[ids]
        )
    )

    links = []

    for k in range(7):

        i = ids[k]
        j = ids[k + 1]

        links.append(
            {
                "logical_i": k,
                "logical_j": k + 1,
                "eigen_i": int(i),
                "eigen_j": int(j),
                "freq_MHz":
                    float(
                        abs(
                            E[j] - E[i]
                        )
                    ),
                "Ix_complex":
                    complex(
                        D[k, k + 1]
                    ),
                "Ix_abs":
                    float(
                        abs(
                            D[k, k + 1]
                        )
                    ),
                "Sz_i":
                    float(
                        Sz_e[i]
                    ),
                "Sz_j":
                    float(
                        Sz_e[j]
                    ),
                "Iz_i":
                    float(
                        Iz_e[i]
                    ),
                "Iz_j":
                    float(
                        Iz_e[j]
                    ),
            }
        )

    return {
        "ids": ids,
        "energies_MHz": E[ids],
        "V": Vb,
        "Ix": D,
        "links": links,
        "Sz_expect": Sz_e[ids],
        "Iz_expect": Iz_e[ids],
    }


def branch_link_table(
    B_T,
    **kwargs,
):
    out = []

    for branch in [
        "ground",
        "excited",
    ]:

        b = branch_basis(
            B_T,
            branch,
            **kwargs,
        )

        for k, row in enumerate(
            b["links"]
        ):

            out.append(
                {
                    "branch": branch,
                    "link": k,
                    **row,
                }
            )

    return out
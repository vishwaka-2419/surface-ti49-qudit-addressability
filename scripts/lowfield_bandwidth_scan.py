from pathlib import Path
import sys
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from scipy.linalg import expm
import matplotlib.pyplot as plt

sys.path.insert(
    0,
    str(
        Path(__file__)
        .resolve()
        .parent
    )
)

from ti49_model_v2 import (
    branch_basis,
    spin,
    I,
)


OUT = Path(
    "analysis_v3/02_control/bandwidth_lowfield"
)

OUT.mkdir(
    parents=True,
    exist_ok=True,
)


R = pd.read_csv(
    "analysis_v3/02_control/lowfield_inputs/"
    "robust_operating_point_summary.csv"
)

Bopt = float(
    R.loc[
        R.label == "robust_optimum",
        "B_T"
    ].iloc[0]
)


W = pd.read_csv(
    "analysis_v3/02_control/lowfield_inputs/"
    "robust_optimum_worst_case_parameters.csv"
).iloc[0]


source_backed = dict(
    Az_MHz=130.0,
    kappa_MHz=-46.373841059602654,
    ge_z=0.56,
    Btip_T=0.0679,
    phi_deg=5.0,

    # 49Ti source-backed transverse
    # hyperfine: 10 MHz / 130 MHz.
    Aperp_over_Az=10.0 / 130.0,

    # complete transverse electron g tensor
    # remains a contextual assumption.
    ge_perp_over_ge_z=1.0,
)


spectral_limiting = dict(
    Az_MHz=float(
        W.Az_MHz
    ),

    kappa_MHz=float(
        W.kappa_MHz
    ),

    ge_z=float(
        W.ge_z
    ),

    Btip_T=float(
        W.Btip_T
    ),

    phi_deg=float(
        W.phi_deg
    ),

    Aperp_over_Az=float(
        W.Aperp_over_Az
    ),

    ge_perp_over_ge_z=float(
        W.ge_perp_over_ge_z
    ),
)


cases = [
    (
        "1p35_sourcebacked",
        1.35,
        source_backed,
    ),

    (
        "lowfield_sourcebacked",
        Bopt,
        source_backed,
    ),

    (
        "lowfield_spectral_limiting_case",
        Bopt,
        spectral_limiting,
    ),
]


fR_grid_kHz = np.arange(5.0, 23.0001, 0.25)


theta = np.pi / 2


_, Iy_ideal, _ = spin(I)

Ulogical = expm(
    -1j
    * theta
    * Iy_ideal
)


def avg_gate_fidelity(
    U,
    V,
):
    d = U.shape[0]

    return float(
        (
            np.trace(
                U.conj().T
                @ U
            ).real
            +
            abs(
                np.trace(
                    V.conj().T
                    @ U
                )
            ) ** 2
        )
        /
        (
            d
            * (
                d + 1
            )
        )
    )


def unitary_similarity(
    U,
    V,
):
    d = U.shape[0]

    return float(
        (
            d
            +
            abs(
                np.trace(
                    V.conj().T
                    @ U
                )
            ) ** 2
        )
        /
        (
            d
            * (
                d + 1
            )
        )
    )


def simulate_case(
    B_T,
    pars,
    fR_kHz,
):

    branches = {
        name:
            branch_basis(
                B_T,
                name,
                **pars,
            )
        for name in [
            "ground",
            "excited",
        ]
    }


    # Build all 14 physical tones.
    tones = []

    for source_branch in [
        "ground",
        "excited",
    ]:

        Bsrc = branches[
            source_branch
        ]

        for k, link in enumerate(
            Bsrc["links"]
        ):

            # +pi/2 produces the
            # Iy phase convention.
            phase = (
                np.angle(
                    link[
                        "Ix_complex"
                    ]
                )
                + np.pi / 2
            )

            tones.append(
                {
                    "source_branch":
                        source_branch,

                    "source_link":
                        k,

                    "freq_MHz":
                        link[
                            "freq_MHz"
                        ],

                    "phase":
                        phase,
                }
            )


    fR_MHz = (
        fR_kHz
        / 1000.0
    )

    tgate_us = (
        theta
        / (
            2
            * np.pi
            * fR_MHz
        )
    )


    Uout = {}


    for branch_name in [
        "ground",
        "excited",
    ]:

        branch = branches[
            branch_name
        ]

        link_freq = np.array(
            [
                x["freq_MHz"]
                for x in branch[
                    "links"
                ]
            ]
        )

        link_d = np.array(
            [
                x["Ix_complex"]
                for x in branch[
                    "links"
                ]
            ],
            dtype=complex,
        )


        def H_int(t):

            H = np.zeros(
                (
                    8,
                    8,
                ),
                dtype=complex,
            )

            for l in range(7):

                fl = link_freq[l]
                dl = link_d[l]

                for tone in tones:

                    det = (
                        fl
                        - tone[
                            "freq_MHz"
                        ]
                    )

                    c = (
                        fR_MHz
                        * dl
                        * np.exp(
                            1j
                            * 2
                            * np.pi
                            * det
                            * t
                            - 1j
                            * tone[
                                "phase"
                            ]
                        )
                    )

                    H[
                        l,
                        l + 1
                    ] += c

                    H[
                        l + 1,
                        l
                    ] += (
                        c.conjugate()
                    )

            return H


        def rhs(
            t,
            y,
        ):

            U = y.reshape(
                8,
                8,
            )

            dU = (
                -1j
                * 2
                * np.pi
                * H_int(t)
                @ U
            )

            return dU.ravel()


        U0 = np.eye(
            8,
            dtype=complex,
        )


        detunings = []

        for fl in link_freq:

            for tone in tones:

                detunings.append(
                    abs(
                        fl
                        - tone[
                            "freq_MHz"
                        ]
                    )
                )

        detunings = np.array(
            detunings
        )

        relevant = detunings[
            detunings < 15.0
        ]

        max_det = max(
            float(
                relevant.max()
            ),
            0.1,
        )


        max_step = min(
            tgate_us / 150.0,
            1.0
            / (
                15.0
                * max_det
            ),
        )


        sol = solve_ivp(
            rhs,
            (
                0.0,
                tgate_us,
            ),
            U0.ravel(),
            method="DOP853",
            rtol=2e-8,
            atol=2e-10,
            max_step=max_step,
        )


        if not sol.success:

            raise RuntimeError(
                sol.message
            )


        U = sol.y[
            :,
            -1
        ].reshape(
            8,
            8,
        )


        Flogical = (
            avg_gate_fidelity(
                U,
                Ulogical,
            )
        )


        Uout[
            branch_name
        ] = {
            "U": U,
            "Flogical":
                Flogical,
            "nfev":
                sol.nfev,
        }


    Fbetween = (
        unitary_similarity(
            Uout[
                "ground"
            ]["U"],
            Uout[
                "excited"
            ]["U"],
        )
    )


    return {
        "gate_time_us":
            tgate_us,

        "ground_logical_Favg":
            Uout[
                "ground"
            ]["Flogical"],

        "excited_logical_Favg":
            Uout[
                "excited"
            ]["Flogical"],

        "worst_branch_logical_Favg":
            min(
                Uout[
                    "ground"
                ]["Flogical"],
                Uout[
                    "excited"
                ]["Flogical"],
            ),

        "branch_to_branch_Favg":
            Fbetween,

        "ground_nfev":
            Uout[
                "ground"
            ]["nfev"],

        "excited_nfev":
            Uout[
                "excited"
            ]["nfev"],
    }


rows = []


for (
    case_name,
    B,
    pars,
) in cases:

    for fR in fR_grid_kHz:

        print(
            f"{case_name:32s} "
            f"fR={fR:7.1f} kHz"
        )

        r = simulate_case(
            B,
            pars,
            float(fR),
        )

        rows.append(
            {
                "case":
                    case_name,

                "B_T":
                    B,

                "fR_kHz":
                    fR,

                **r,
            }
        )


df = pd.DataFrame(rows)

df.to_csv(
    OUT
    / "dual_manifold_14tone_raw.csv",
    index=False,
)


fig, ax = plt.subplots(
    figsize=(7.2, 4.8)
)

for case, g in df.groupby(
    "case"
):

    ax.plot(
        g.fR_kHz,
        g.worst_branch_logical_Favg,
        marker="o",
        label=case,
    )

ax.set_xlabel(
    "Effective nuclear generator rate "
    r"$f_R$ (kHz)"
)

ax.set_ylabel(
    r"Worst-branch logical "
    r"$I_y(\pi/2)$ fidelity"
)

ax.set_ylim(
    0.0,
    1.01,
)

ax.set_title(
    "Simultaneous dual-manifold "
    "14-tone control"
)

ax.legend(
    frameon=False,
)

fig.tight_layout()

fig.savefig(
    OUT
    / "dual_manifold_logical_fidelity.png",
    dpi=250,
)

plt.close(fig)


fig, ax = plt.subplots(
    figsize=(7.2, 4.8)
)

for case, g in df.groupby(
    "case"
):

    ax.plot(
        g.fR_kHz,
        g.branch_to_branch_Favg,
        marker="o",
        label=case,
    )

ax.set_xlabel(
    "Effective nuclear generator rate "
    r"$f_R$ (kHz)"
)

ax.set_ylabel(
    "Ground/excited branch "
    "unitary similarity"
)

ax.set_ylim(
    0.0,
    1.01,
)

ax.set_title(
    "Electron-manifold consistency "
    "of nuclear rotation"
)

ax.legend(
    frameon=False,
)

fig.tight_layout()

fig.savefig(
    OUT
    / "dual_manifold_branch_consistency.png",
    dpi=250,
)

plt.close(fig)


print("\nRESULTS")
print(
    df.to_string(
        index=False
    )
)
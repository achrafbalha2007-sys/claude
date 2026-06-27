# -*- coding: utf-8 -*-
"""
TIPE - Modele SIMPLIFIE (poutre 1D) - version slides.
Probabilite de rupture d'un ecran protecteur lors d'une chute, par Monte-Carlo.

Notions utilisees (programme MP/MPI uniquement) :
  - mecanique du point : chute libre, v = sqrt(2 g h), Ek = 1/2 m v^2
  - oscillateur / ressort : le verre stocke 1/2 k x^2
  - flexion d'une poutre (RDM de base) : I = b e^3 / 12, sigma = 6 E e x / L^2
  - probabilites : tirage aleatoire + Monte-Carlo

Dependances : numpy, matplotlib
"""

import numpy as np
import matplotlib.pyplot as plt

rng = np.random.default_rng(2026)

# --- constantes (verre trempe + smartphone) -------------------------------
E, G = 73e9, 9.81          # module d'Young (Pa), pesanteur
M_PHONE = 0.18             # masse du telephone (kg)
ETA = 0.05                 # fraction de l'energie de chute stockee en flexion
B = 0.020                  # largeur de la bande de verre qui flechit (m)
SIG_MU, SIG_SD = 850e6, 120e6   # seuil de rupture : gaussienne (Pa)
KAPPA = 0.7                # surcontrainte au coin
L_C, L_B = 0.050, 0.025    # portee effective : centre / coin (m)
E0 = 0.5e-3                # epaisseur uniforme de reference (m)
N_MC = 5000


# --- profils d'epaisseur e(r), r = distance normalisee au centre (0..1) ---
def uniforme(r, e0=E0):
    return np.full_like(np.asarray(r, float), e0)

def parabolique(r, e_c=0.5e-3, e_b=1.2e-3):
    return e_c + (e_b - e_c) * np.asarray(r, float) ** 2

def exponentiel(r, e_c=0.5e-3, e_b=1.2e-3):
    return e_c * np.exp(np.log(e_b / e_c) * np.asarray(r, float))

def epaisseur_moyenne(prof, **kw):
    r = np.linspace(0, 1, 400)
    return np.trapz(prof(r, **kw) * 2 * r, r) / np.trapz(2 * r, r)


# --- coeur du modele : contrainte de flexion (poutre 1D) ------------------
def contrainte(e, L, Ek, Kc):
    """sigma = 6 E e x_max / L^2, avec x_max issu de l'energie."""
    I = B * e ** 3 / 12.0           # moment quadratique (RDM)
    k = 48.0 * E * I / L ** 3       # raideur poutre sur 2 appuis : k ~ e^3
    x = np.sqrt(2.0 * ETA * Ek / k) # 1/2 k x^2 = eta Ek  ->  fleche max
    return 6.0 * E * e * x / L ** 2 * Kc


def gauss_tronquee(mu, s, lo, hi, n):
    out = np.empty(n); k = 0
    while k < n:
        x = rng.normal(mu, s, n - k); x = x[(x >= lo) & (x <= hi)]
        out[k:k + x.size] = x; k += x.size
    return out


def monte_carlo(prof, N=N_MC, **kw):
    theta = gauss_tronquee(65, 15, 10, 90, N)        # angle de chute (deg)
    h = gauss_tronquee(1.0, 0.25, 0.30, 2.0, N)      # hauteur (m)
    v = np.sqrt(2 * G * h); Ek = 0.5 * M_PHONE * v ** 2
    r = np.clip(theta / 90 + rng.normal(0, 0.05, N), 0, 1)   # point d'impact
    L = L_C - (L_C - L_B) * r                        # portee selon la zone
    Kc = 1 + KAPPA * r                               # coin = surcontrainte
    e = prof(r, **kw)
    sigma = contrainte(e, L, Ek, Kc)
    sigma_c = np.maximum(rng.normal(SIG_MU, SIG_SD, N), 1e6)  # seuil aleatoire
    return (sigma > sigma_c).mean(), sigma


def optimise():
    emu = epaisseur_moyenne(uniforme)
    ec_g = np.linspace(0.2e-3, 0.5e-3, 7)
    eb_g = np.linspace(0.5e-3, 2.0e-3, 16)
    P = np.full((ec_g.size, eb_g.size), np.nan); best = None
    for i, ec in enumerate(ec_g):
        for j, eb in enumerate(eb_g):
            if eb < ec or epaisseur_moyenne(parabolique, e_c=ec, e_b=eb) > 2.5 * emu:
                continue
            p, _ = monte_carlo(parabolique, N=3000, e_c=ec, e_b=eb)
            P[i, j] = p
            if best is None or p < best[0]:
                best = (p, ec, eb)
    return ec_g, eb_g, P, best


# --- figures --------------------------------------------------------------
def fig_profils_et_P():
    profs = [("Uniforme", uniforme, {"e0": E0}),
             ("Parabolique", parabolique, {"e_c": 0.5e-3, "e_b": 1.2e-3}),
             ("Exponentiel", exponentiel, {"e_c": 0.5e-3, "e_b": 1.2e-3})]
    r = np.linspace(0, 1, 100); fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    noms, Pr = [], []
    for nom, f, kw in profs:
        ax[0].plot(r, f(r, **kw) * 1e3, lw=2, label=nom)
        p, _ = monte_carlo(f, **kw); noms.append(nom); Pr.append(p * 100)
    ax[0].set_xlabel("r (centre -> bord)"); ax[0].set_ylabel("epaisseur (mm)")
    ax[0].set_title("Profils"); ax[0].grid(True); ax[0].legend()
    ax[1].bar(noms, Pr, color=["#888", "#3a7", "#37a"])
    ax[1].set_ylabel("P(rupture) (%)"); ax[1].set_title("Risque de casse")
    for i, p in enumerate(Pr): ax[1].text(i, p + 1, f"{p:.0f}%", ha="center")
    plt.tight_layout(); plt.savefig("fig_poutre_profils.png", dpi=160)

def fig_sigma_de_e():
    e = np.linspace(0.3e-3, 1.5e-3, 25)
    s = [contrainte(ei, L_C, 0.5 * M_PHONE * 2 * G * 1.0, 1.0) / 1e6 for ei in e]
    plt.figure(figsize=(6, 4)); plt.plot(e * 1e3, s, "o-")
    plt.axhline(SIG_MU / 1e6, ls="--", color="r", label="seuil ~850 MPa")
    plt.xlabel("epaisseur e (mm)"); plt.ylabel("contrainte (MPa)")
    plt.title(r"$\sigma \propto 1/\sqrt{e}$ : epaissir = moins de contrainte")
    plt.grid(True); plt.legend(); plt.tight_layout()
    plt.savefig("fig_poutre_sigma.png", dpi=160)

def fig_optim(ec_g, eb_g, P):
    plt.figure(figsize=(7, 4))
    c = plt.pcolormesh(eb_g * 1e3, ec_g * 1e3, P * 100, shading="auto", cmap="viridis")
    plt.colorbar(c, label="P(rupture) (%)")
    plt.xlabel("epaisseur de bord e_b (mm)"); plt.ylabel("epaisseur centrale e_c (mm)")
    plt.title("Carte d'optimisation"); plt.tight_layout()
    plt.savefig("fig_poutre_optim.png", dpi=160)


if __name__ == "__main__":
    v1 = np.sqrt(2 * G * 1.0); Ek1 = 0.5 * M_PHONE * v1 ** 2
    print(f"[h=1m] v={v1:.2f} m/s  Ek={Ek1:.2f} J  "
          f"sigma_centre={contrainte(E0, L_C, Ek1, 1.0)/1e6:.0f} MPa")
    for nom, f, kw in [("Uniforme 0.5mm", uniforme, {"e0": E0}),
                       ("Parabolique", parabolique, {"e_c": 0.5e-3, "e_b": 1.2e-3}),
                       ("Exponentiel", exponentiel, {"e_c": 0.5e-3, "e_b": 1.2e-3})]:
        p, s = monte_carlo(f, **kw)
        print(f"[MC] {nom:16s} P={p*100:5.1f}%  <sigma>={s.mean()/1e6:.0f} MPa")
    ec_g, eb_g, P, best = optimise()
    print(f"[OPT] P={best[0]*100:.1f}%  e_c={best[1]*1e3:.2f}mm  e_b={best[2]*1e3:.2f}mm")
    fig_profils_et_P(); fig_sigma_de_e(); fig_optim(ec_g, eb_g, P)
    print("Figures: fig_poutre_profils.png, fig_poutre_sigma.png, fig_poutre_optim.png")

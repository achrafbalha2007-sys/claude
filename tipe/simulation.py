# -*- coding: utf-8 -*-
"""
TIPE - Optimisation du profil d'epaisseur d'un ecran protecteur en verre trempe.
Simulation de Monte-Carlo de la probabilite de rupture lors d'une chute.

Auteur  : <a completer>
Filiere : MP/MPI
Dependances : numpy, matplotlib, scipy

Modele (detaille dans le dossier, section 5) :
  - Angle de chute tire selon une gaussienne tronquee centree a 65 deg.
  - Energie cinetique d'impact Ek = 1/2 m v^2, v = sqrt(2 g h).
  - Plaque mince (Kirchhoff-Love) reduite a un oscillateur 1 ddl equivalent.
  - Methode energetique : une fraction eta de Ek est stockee en flexion locale.
  - Contrainte de traction max sur une plaque circulaire encastree equivalente.
  - Rupture fragile : seuil tire selon une loi de Weibull (verre trempe).
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize

rng = np.random.default_rng(2026)

# ===========================================================================
# 1. CONSTANTES PHYSIQUES  (verre aluminosilicate trempe chimiquement)
#    Sources : cf. bibliographie du dossier (Corning Gorilla Glass, Gy 2008).
# ===========================================================================
E   = 73e9        # module d'Young (Pa)
NU  = 0.22        # coefficient de Poisson
RHO = 2450.0      # masse volumique (kg/m^3)
G   = 9.81        # acceleration de la pesanteur (m/s^2)

M_PHONE = 0.18    # masse du smartphone (kg)
E_REST  = 0.30    # coefficient de restitution effectif (-)
ETA     = 0.10    # fraction de Ek injectee en flexion locale (HYPOTHESE FORTE)

# Seuil de rupture en traction de surface (verre trempe) : loi de Weibull
SIGMA_THETA = 850e6   # contrainte caracteristique (Pa)
M_WEIBULL   = 8.0     # module de Weibull (-)

# Geometrie de reponse locale : rayon de la zone qui flechit reellement
A_CENTER = 0.030  # impact a plat -> grande zone mobilisee (m)
A_CORNER = 0.008  # impact au coin -> petite zone mobilisee (m)
KAPPA    = 1.5    # facteur de concentration de contrainte au coin (-)

# Geometrie de l'ecran
L, W, R_COIN = 0.150, 0.070, 0.012   # longueur, largeur, rayon de coin (m)
E0 = 0.5e-3                          # epaisseur uniforme de reference (m)

N_MC = 5000       # tirages de Monte-Carlo


# ===========================================================================
# 2. PROFILS D'EPAISSEUR e(r)   (r : distance normalisee au centre, 0..1)
# ===========================================================================
def profil_uniforme(r, e0=E0):
    return np.full_like(np.asarray(r, dtype=float), e0)

def profil_parabolique(r, e_c=0.5e-3, e_b=1.2e-3):
    return e_c + (e_b - e_c) * np.asarray(r, dtype=float) ** 2

def profil_exponentiel(r, e_c=0.5e-3, e_b=1.2e-3):
    alpha = np.log(e_b / e_c)
    return e_c * np.exp(alpha * np.asarray(r, dtype=float))


def epaisseur_moyenne(profil, **kw):
    """Moyenne de e ponderee par l'aire (poids 2r dr) -> proportionnelle a la masse."""
    r = np.linspace(0.0, 1.0, 400)
    return np.trapz(profil(r, **kw) * 2 * r, r) / np.trapz(2 * r, r)


# ===========================================================================
# 3. MECANIQUE DU CHOC
# ===========================================================================
def rigidite_flexion(e):
    """Rigidite de flexion D de la plaque (N.m)."""
    return E * e ** 3 / (12.0 * (1.0 - NU ** 2))

def contrainte_max(e_loc, a, Ek, Kc):
    """
    Contrainte de traction max (Pa).
    Methode energetique : 1/2 k w0^2 = eta Ek, avec k = 16 pi D / a^2.
    Plaque circulaire encastree, charge centrale :
        sigma = 2 E e w0 / ((1-nu^2) a^2), majoree par Kc (coin).
    """
    D = rigidite_flexion(e_loc)
    k = 16.0 * np.pi * D / a ** 2
    w0 = np.sqrt(2.0 * ETA * Ek / k)
    sigma = 2.0 * E * e_loc * w0 / ((1.0 - NU ** 2) * a ** 2)
    return sigma * Kc

def temps_contact_hertz(v, R_c=5e-3):
    """Theorie de Hertz : duree de contact, force et indentation max."""
    Estar = E / (1.0 - NU ** 2)
    kh = (4.0 / 3.0) * Estar * np.sqrt(R_c)
    delta = (5.0 * M_PHONE * v ** 2 / (4.0 * kh)) ** 0.4
    Fmax = kh * delta ** 1.5
    tc = 2.87 * (M_PHONE ** 2 / (R_c * Estar ** 2 * v)) ** 0.2
    return tc, Fmax, delta


def gaussienne_tronquee(mu, sigma, lo, hi, n):
    """Tirage gaussien tronque par rejet."""
    out = np.empty(n)
    k = 0
    while k < n:
        x = rng.normal(mu, sigma, n - k)
        x = x[(x >= lo) & (x <= hi)]
        out[k:k + x.size] = x
        k += x.size
    return out


# ===========================================================================
# 4. SIMULATION DE MONTE-CARLO
# ===========================================================================
def monte_carlo(profil, N=N_MC, **kw):
    theta = gaussienne_tronquee(65.0, 15.0, 10.0, 90.0, N)        # deg
    h = gaussienne_tronquee(1.0, 0.25, 0.30, 2.0, N)              # m
    v = np.sqrt(2.0 * G * h)
    Ek = 0.5 * M_PHONE * v ** 2

    r = np.clip(theta / 90.0 + rng.normal(0, 0.05, N), 0.0, 1.0)  # position d'impact
    a = A_CENTER - (A_CENTER - A_CORNER) * r
    Kc = 1.0 + KAPPA * r
    e_loc = profil(r, **kw)

    sigma = contrainte_max(e_loc, a, Ek, Kc)
    # seuil de Weibull tire par inversion : sigma_c = sigma_theta (-ln U)^(1/m)
    U = rng.random(N)
    sigma_c = SIGMA_THETA * (-np.log(1.0 - U)) ** (1.0 / M_WEIBULL)

    rupture = sigma > sigma_c
    return rupture.mean(), sigma, r


# ===========================================================================
# 5. OPTIMISATION DU PROFIL PARABOLIQUE
# ===========================================================================
def optimise_parabolique():
    e_moy_unif = epaisseur_moyenne(profil_uniforme)
    ec_grid = np.linspace(0.2e-3, 0.5e-3, 7)     # contrainte tactilite : ec <= 0.5 mm
    eb_grid = np.linspace(0.5e-3, 2.0e-3, 16)
    P = np.full((ec_grid.size, eb_grid.size), np.nan)
    best = None
    for i, ec in enumerate(ec_grid):
        for j, eb in enumerate(eb_grid):
            if eb < ec:
                continue
            emoy = epaisseur_moyenne(profil_parabolique, e_c=ec, e_b=eb)
            if emoy > 2.5 * e_moy_unif:           # contrainte de masse
                continue
            p, _, _ = monte_carlo(profil_parabolique, N=3000, e_c=ec, e_b=eb)
            P[i, j] = p
            if best is None or p < best[0]:
                best = (p, ec, eb, emoy)
    return ec_grid, eb_grid, P, best


# ===========================================================================
# 6. PRODUCTION DES FIGURES
# ===========================================================================
def figure_carte_pression():
    nx, ny = 200, 120
    xs = np.linspace(-1, 1, nx)
    ys = np.linspace(-1, 1, ny)
    X, Y = np.meshgrid(xs, ys)
    Rn = np.clip(np.sqrt(X ** 2 + (Y * W / L) ** 2), 0, 1)
    a = A_CENTER - (A_CENTER - A_CORNER) * Rn
    Kc = 1.0 + KAPPA * Rn
    e_loc = profil_uniforme(Rn)
    Ek = 0.5 * M_PHONE * (2 * G * 1.0)        # h = 1 m
    sig = contrainte_max(e_loc, a, Ek, Kc) / 1e6
    plt.figure(figsize=(7, 4))
    c = plt.pcolormesh(X * L / 2 * 1e3, Y * L / 2 * 1e3, sig, shading="auto", cmap="inferno")
    plt.colorbar(c, label="Contrainte de traction (MPa)")
    plt.title("Carte de contrainte (profil uniforme, h = 1 m)")
    plt.xlabel("x (mm)"); plt.ylabel("y (mm)"); plt.axis("equal")
    plt.tight_layout(); plt.savefig("fig_carte_pression.png", dpi=160)

def figure_P_de_e():
    e_vals = np.linspace(0.3e-3, 1.5e-3, 25)
    Ps = [monte_carlo(profil_uniforme, e0=e)[0] for e in e_vals]
    plt.figure(figsize=(6, 4))
    plt.plot(e_vals * 1e3, np.array(Ps) * 100, "o-")
    plt.axvline(0.5, ls="--", color="r", label="limite tactilite (0,5 mm)")
    plt.xlabel("Epaisseur uniforme e (mm)"); plt.ylabel("P(rupture) (%)")
    plt.title("Probabilite de rupture vs epaisseur (profil uniforme)")
    plt.grid(True); plt.legend(); plt.tight_layout()
    plt.savefig("fig_P_de_e.png", dpi=160)

def figure_comparaison_profils():
    profs = [("Uniforme", profil_uniforme, {"e0": E0}),
             ("Parabolique", profil_parabolique, {"e_c": 0.5e-3, "e_b": 1.2e-3}),
             ("Exponentiel", profil_exponentiel, {"e_c": 0.5e-3, "e_b": 1.2e-3})]
    r = np.linspace(0, 1, 100)
    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    noms, Pr = [], []
    for nom, f, kw in profs:
        ax[0].plot(r, f(r, **kw) * 1e3, label=nom)
        p, _, _ = monte_carlo(f, **kw)
        noms.append(nom); Pr.append(p * 100)
    ax[0].set_xlabel("r (distance normalisee au centre)")
    ax[0].set_ylabel("epaisseur (mm)"); ax[0].set_title("Profils testes")
    ax[0].grid(True); ax[0].legend()
    ax[1].bar(noms, Pr, color=["#888", "#3a7", "#37a"])
    ax[1].set_ylabel("P(rupture) (%)"); ax[1].set_title("Comparaison des profils")
    for i, p in enumerate(Pr):
        ax[1].text(i, p + 1, f"{p:.0f}%", ha="center")
    plt.tight_layout(); plt.savefig("fig_comparaison_profils.png", dpi=160)

def figure_optimisation(ec_grid, eb_grid, P):
    plt.figure(figsize=(7, 4))
    c = plt.pcolormesh(eb_grid * 1e3, ec_grid * 1e3, P * 100, shading="auto", cmap="viridis")
    plt.colorbar(c, label="P(rupture) (%)")
    plt.xlabel("epaisseur de bord e_b (mm)")
    plt.ylabel("epaisseur centrale e_c (mm)")
    plt.title("Carte d'optimisation du profil parabolique")
    plt.tight_layout(); plt.savefig("fig_optimisation.png", dpi=160)


# ===========================================================================
# 7. PROGRAMME PRINCIPAL
# ===========================================================================
if __name__ == "__main__":
    v1 = np.sqrt(2 * G * 1.0)
    tc, Fmax, delta = temps_contact_hertz(v1)
    D0 = rigidite_flexion(E0)
    f11 = (np.pi / 2) * np.sqrt(D0 / (RHO * E0)) * (1 / L ** 2 + 1 / W ** 2)
    print(f"[Choc h=1m] v={v1:.2f} m/s  Ek={0.5*M_PHONE*v1**2:.2f} J")
    print(f"[Hertz] t_contact={tc*1e3:.3f} ms  F_max={Fmax/1e3:.1f} kN  delta={delta*1e3:.3f} mm")
    print(f"[Plaque] f_11={f11:.0f} Hz  T_11={1e3/f11:.2f} ms  ->  t_c/T_11={tc*f11:.3f}")

    for nom, f, kw in [("Uniforme 0.5mm", profil_uniforme, {"e0": E0}),
                       ("Parabolique", profil_parabolique, {"e_c": 0.5e-3, "e_b": 1.2e-3}),
                       ("Exponentiel", profil_exponentiel, {"e_c": 0.5e-3, "e_b": 1.2e-3})]:
        p, sig, _ = monte_carlo(f, **kw)
        print(f"[MC] {nom:16s}: P(rupture)={p*100:5.1f}%  <sigma>={sig.mean()/1e6:.0f} MPa "
              f" e_moy={epaisseur_moyenne(f, **kw)*1e3:.3f} mm")

    ec_grid, eb_grid, P, best = optimise_parabolique()
    print(f"[OPT] Optimum parabolique : P={best[0]*100:.1f}%  "
          f"e_c={best[1]*1e3:.2f} mm  e_b={best[2]*1e3:.2f} mm  e_moy={best[3]*1e3:.3f} mm")

    figure_carte_pression()
    figure_P_de_e()
    figure_comparaison_profils()
    figure_optimisation(ec_grid, eb_grid, P)
    print("Figures enregistrees : fig_carte_pression.png, fig_P_de_e.png, "
          "fig_comparaison_profils.png, fig_optimisation.png")

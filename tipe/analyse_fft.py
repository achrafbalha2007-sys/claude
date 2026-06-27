# -*- coding: utf-8 -*-
"""
TIPE - Validation acoustique de la rigidite locale d'un ecran protecteur.

Principe : on excite la plaque par un impact leger en differents points d'une
grille (3x3). Le son rayonne a une frequence d'autant plus BASSE que la zone
est SOUPLE (faible epaisseur / faible precontrainte) -> f ~ sqrt(k/m).
On compare la carte de frequence propre locale a la carte de contrainte
theorique : une zone basse frequence (souple) correspond a une zone de forte
contrainte au choc.

Dependances : numpy, scipy, matplotlib
Entree      : fichiers WAV mono enregistres au smartphone (un par point).

REMARQUE (a souligner devant le jury) : l'acoustique renseigne sur la rigidite
RELATIVE, pas sur le seuil de rupture ABSOLU. C'est une validation qualitative
du classement des zones, pas une mesure de la contrainte a rupture.
"""

import sys
import glob
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import wavfile
from scipy.signal import find_peaks


def frequence_dominante(chemin_wav, f_min=200.0, f_max=12000.0):
    """Renvoie la frequence du pic spectral dominant (Hz) d'un .wav, par FFT."""
    fs, x = wavfile.read(chemin_wav)
    if x.ndim > 1:                      # stereo -> mono
        x = x.mean(axis=1)
    x = x.astype(float)
    x -= x.mean()

    # On se cale sur la fenetre qui suit l'attaque de l'impact (transitoire)
    i0 = int(np.argmax(np.abs(x)))
    seg = x[i0:i0 + int(0.20 * fs)]     # 200 ms apres l'impact
    if seg.size < 1024:
        seg = x
    seg = seg * np.hanning(seg.size)    # fenetrage pour limiter les fuites

    spectre = np.abs(np.fft.rfft(seg))
    freqs = np.fft.rfftfreq(seg.size, 1.0 / fs)

    masque = (freqs >= f_min) & (freqs <= f_max)
    f_band, s_band = freqs[masque], spectre[masque]

    pics, _ = find_peaks(s_band, height=0.2 * s_band.max())
    if pics.size == 0:
        f_dom = f_band[np.argmax(s_band)]
    else:
        f_dom = f_band[pics[np.argmax(s_band[pics])]]

    # incertitude type = resolution spectrale df = fs / N
    df = fs / seg.size
    return f_dom, df, (f_band, s_band)


def carte_frequences(fichiers_grille):
    """
    fichiers_grille : dict {(i, j): 'chemin.wav'} pour une grille 3x3.
    Position (0,0) = coin, (1,1) = centre, etc.
    Renvoie une matrice 3x3 de frequences dominantes (Hz).
    """
    F = np.full((3, 3), np.nan)
    U = np.full((3, 3), np.nan)
    for (i, j), chemin in fichiers_grille.items():
        f_dom, df, _ = frequence_dominante(chemin)
        F[i, j] = f_dom
        U[i, j] = df
        print(f"point ({i},{j}) {chemin:30s} -> f = {f_dom:7.1f} Hz  (u = {df:.1f} Hz)")
    return F, U


def correlation_rigidite_contrainte(F, sigma_theo):
    """
    Coefficient de correlation entre la frequence locale et la contrainte
    theorique. On attend une correlation NEGATIVE : basse frequence (souple)
    <-> forte contrainte au choc.
    """
    f = F.flatten()
    s = sigma_theo.flatten()
    ok = ~np.isnan(f) & ~np.isnan(s)
    r = np.corrcoef(f[ok], s[ok])[0, 1]
    return r


def figure_cartes(F, sigma_theo):
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    c0 = ax[0].imshow(F, cmap="viridis_r")
    ax[0].set_title("Frequence propre locale mesuree (Hz)")
    fig.colorbar(c0, ax=ax[0])
    c1 = ax[1].imshow(sigma_theo, cmap="inferno")
    ax[1].set_title("Contrainte theorique au choc (MPa)")
    fig.colorbar(c1, ax=ax[1])
    for a in ax:
        a.set_xticks([0, 1, 2]); a.set_yticks([0, 1, 2])
        a.set_xticklabels(["gauche", "milieu", "droite"])
        a.set_yticklabels(["haut", "milieu", "bas"])
    plt.tight_layout(); plt.savefig("fig_carte_acoustique.png", dpi=160)


if __name__ == "__main__":
    fichiers = dict(sorted({
        (i, j): f"impact_{i}{j}.wav"
        for i in range(3) for j in range(3)
        if glob.glob(f"impact_{i}{j}.wav")
    }.items()))

    if not fichiers:
        print("Aucun fichier impact_ij.wav trouve. Placez les enregistrements")
        print("d'impact dans le dossier (impact_00.wav ... impact_22.wav).")
        sys.exit(0)

    F, U = carte_frequences(fichiers)

    # Carte de contrainte theorique correspondante (MPa) issue de simulation.py
    # (a remplacer par la sortie reelle du modele ; ici structure attendue :
    #  forte au centre/bords minces, plus faible la ou c'est rigide).
    sigma_theo = np.array([[820, 560, 800],
                           [610, 350, 600],
                           [810, 570, 815]], dtype=float)

    r = correlation_rigidite_contrainte(F, sigma_theo)
    print(f"\nCorrelation f_locale / sigma_theorique : r = {r:+.2f}")
    print("Attendu : r < 0 (zone souple = basse frequence = forte contrainte).")
    figure_cartes(F, sigma_theo)
    print("Figure enregistree : fig_carte_acoustique.png")

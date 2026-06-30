
=====================================================================
%  Demarche theorique - sequence de slides (Beamer)
%  Compilable seul : pdflatex demarche_theorique.tex
% =====================================================================
\documentclass[11pt,aspectratio=169]{beamer}
\usetheme{Madrid}\usecolortheme{whale}
\setbeamertemplate{navigation symbols}{}
\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
\usepackage[french]{babel}
\usepackage{amsmath,amssymb}
\usepackage{tikz}
\usetikzlibrary{arrows.meta,positioning,decorations.pathmorphing,calc}
\definecolor{verre}{RGB}{30,120,180}
\definecolor{danger}{RGB}{200,40,40}
\definecolor{ok}{RGB}{30,150,90}
\definecolor{gristik}{RGB}{120,120,120}

\begin{document}

\section{Démarche théorique}

% ---------- 1. Modelisation ----------
\begin{frame}{Modélisation : de l'écran à une poutre}
  \begin{columns}[c]
  \column{0.46\textwidth}
    \centering
    \begin{tikzpicture}[scale=0.72]
      \draw[verre,thick](0,3.2) rectangle (4,5.6);
      \fill[verre!20](0,4.1) rectangle (4,4.6);
      
ode[font=\tiny] at (2,5.9){vue de dessus};
      
ode[font=\tiny,right] at (4.05,4.35){bande $b$};
      \draw[-{Latex}](2,4.0)--(2,3.0);
      \fill[gristik](0.4,2.5)--(0.18,2.2)--(0.62,2.2)--cycle;
      \fill[gristik](3.6,2.5)--(3.38,2.2)--(3.82,2.2)--cycle;
      \draw[verre,very thick](0.4,2.6)--(3.6,2.6);
      \draw[-{Latex},danger,thick](2,3.0)--(2,2.7);
      \draw[{Latex}-{Latex}](0.4,2.05)--(3.6,2.05) node[midway,below,font=\tiny]{portée $L$};
      
ode[font=\tiny] at (2,1.5){vue de côté (poutre)};
    \end{tikzpicture}
  \column{0.54\textwidth}
    \footnotesize
    On étudie un écran rectangulaire ($L\approx 150$~mm, $\ell\approx 70$~mm),
    d'épaisseur $e\approx 0{,}5$~mm, donc \textbf{mince} ($e\ll L,\ell$).
    \medskip\\
    On isole une bande de largeur $b$ et on la modélise par une \textbf{poutre
    en flexion} sur deux appuis.
    \medskip\\
    \textbf{Hypothèses :}
    \begin{itemize}\itemsep1pt
      \item élasticité linéaire isotrope, petites déformations ;
      \item matériau \textbf{fragile} (rupture sans plasticité) ;
      \item poutre d'Euler--Bernoulli (cisaillement négligé, $e\ll L$) ;
      \item choc traité par un \textbf{bilan d'énergie}.
    \end{itemize}
  \end{columns}
\end{frame}

% ---------- 2. Cinematique ----------
\begin{frame}{Hypothèse cinématique : les sections restent planes}
  \begin{columns}[c]
  \column{0.46\textwidth}
    \centering
    \begin{tikzpicture}[scale=1]
      \coordinate (O) at (0,3.1);
      
ode[font=\tiny,above] at (O){$O$ (centre de courbure)};
      \draw[gray] (O)--($(O)+(-112:3.1)$);
      \draw[gray] (O)--($(O)+(-68:3.1)$);
      \draw[verre,very thick] ($(O)+(-112:2.3)$) arc (-112:-68:2.3);
      \draw[danger,thick] ($(O)+(-112:2.9)$) arc (-112:-68:2.9);
      \draw[{Latex}-{Latex}] ($(O)+(-90:2.3)$)--($(O)+(-90:2.9)$)
            node[midway,right,font=\tiny]{$z$};
      
ode[verre,font=\tiny,left] at ($(O)+(-112:2.3)$){fibre neutre};
      
ode[danger,font=\tiny,left] at ($(O)+(-112:2.9)$){fibre tendue};
      \draw ($(O)+(-100:0.9)$) arc (-100:-80:0.9);
      
ode[font=\tiny] at ($(O)+(-90:1.25)$){$\mathrm{d}\theta$};
      
ode[font=\tiny,left] at ($(O)+(-90:1.7)$){$R$};
    \end{tikzpicture}
  \column{0.54\textwidth}
    \footnotesize
    On suppose qu'une section droite reste \textbf{plane et perpendiculaire} à la
    fibre moyenne après flexion.
    \medskip\\
    Une fibre à la distance $z$ de la \textbf{fibre neutre} a une longueur
    $(R+z)\,\mathrm{d}\theta$, alors que la fibre neutre garde $R\,\mathrm{d}\theta$.
    Sa déformation vaut donc :
    \[ \varepsilon(z)=\frac{(R+z)\mathrm{d}\theta-R\,\mathrm{d}\theta}{R\,\mathrm{d}\theta}
       =\frac{z}{R}=z\,\kappa \]
    où $\kappa=1/R$ est la \textbf{courbure}. La déformation est
    \textbf{proportionnelle à $z$}.
  \end{columns}
\end{frame}

% ---------- 3. Hooke -> contrainte locale ----------
\begin{frame}{Loi de Hooke et profil de contrainte}
  \begin{columns}[c]
  \column{0.52\textwidth}
    \centering
    \begin{tikzpicture}[scale=1]
      \draw[verre,thick,fill=verre!10](0,-1.5) rectangle (0.7,1.5);
      \draw[{Latex}-{Latex}](-0.35,-1.5)--(-0.35,1.5) node[midway,left,font=\tiny]{$e$};
      \draw[dashed](0,0)--(3.3,0);
      
ode[font=\tiny,right] at (2.1,0){fibre neutre, $\sigma=0$};
      \foreach \z/\l in {0.5/0.45,1.0/0.9,1.5/1.35}{
        \draw[-{Latex},ok](0.7+\l,\z)--(0.7,\z);}
      \foreach \z/\l in {0.5/0.45,1.0/0.9,1.5/1.35}{
        \draw[-{Latex},danger](0.7,-\z)--(0.7+\l,-\z);}
      
ode[ok,font=\tiny,right] at (1.9,1.5){comprimé};
      
ode[danger,font=\tiny,right] at (1.9,-1.5){tendu : rupture};
    \end{tikzpicture}
  \column{0.48\textwidth}
    \footnotesize
    La loi de Hooke relie contrainte et déformation :
    \[ \sigma(z)=E\,\varepsilon(z)=E\,z\,\kappa . \]
    La contrainte est donc \textbf{linéaire} dans l'épaisseur : nulle au milieu,
    \textbf{maximale en surface} ($z=\pm e/2$) :
    \[ \boxed{\;\sigma_{\max}=E\,\frac{e}{2}\,\kappa\;} \]
    La face tendue est le siège de la rupture.
  \end{columns}
\end{frame}

% ---------- 4. Moment flechissant ----------
\begin{frame}{Du champ de contrainte au moment fléchissant}
  \footnotesize
  Le \textbf{moment fléchissant} est le moment résultant des contraintes sur la
  section :
  \[ M=\int_{S}\sigma(z)\,z\,\mathrm{d}A
      =E\,\kappa\int_{S} z^{2}\,\mathrm{d}A=E\,I\,\kappa,
      \qquad I=\int_{S}z^{2}\,\mathrm{d}A=\frac{b\,e^{3}}{12}. \]
  On obtient ainsi deux écritures de la contrainte de surface :
  \[ \sigma_{\max}=E\,\frac{e}{2}\,\kappa \;=\; \frac{M\,(e/2)}{I}
      \;=\;\frac{6M}{b\,e^{2}} . \]
  \begin{block}{Point de vigilance}
  Ces deux formes sont \textbf{la même relation} (car $M=EI\kappa$). Elles ne
  donnent \emph{pas} la dépendance de $\sigma$ en $e$ : ni $M$ ni $\kappa$ ne
  sont des données du problème, ils dépendent du chargement, donc de $e$.
  \end{block}
\end{frame}

% ---------- 5. Ce que la chute impose ----------
\begin{frame}{Ce que la chute impose réellement}
  \begin{columns}[c]
  \column{0.4\textwidth}
    \centering
    \begin{tikzpicture}[scale=0.9]
      \draw[-{Latex}](0,3.4)--(0,-0.1);
      \draw[verre,fill=verre!12,thick,rounded corners=2pt](1.1,2.7) rectangle (1.7,3.5);
      \draw[-{Latex[length=2mm]},danger,thick](1.4,2.5)--(1.4,1.7) node[midway,right,font=\tiny]{$v$};
      \draw[thick](0.5,0)--(2.6,0);
      \foreach \x in {0.6,0.9,1.2,1.5,1.8,2.1,2.4}{\draw[gristik](\x,0)--(\x-0.2,-0.2);}
      \draw[{Latex}-{Latex}](0.15,0)--(0.15,3.1) node[midway,left,font=\tiny]{$h$};
    \end{tikzpicture}
  \column{0.6\textwidth}
    \footnotesize
    \begin{block}{Imposé par la chute}
    L'\textbf{énergie} $E_c=\tfrac12 mv^{2}=mgh$ et la \textbf{quantité de
    mouvement} $p=mv$ (fixées par $h$).
    \end{block}
    \textbf{Non imposés} : le moment $M$, la courbure $\kappa$, la force de
    contact $F$. Ce sont des \emph{conséquences} qui dépendent de la raideur de
    l'écran, donc de $e$.
    \medskip\\
    \textcolor{danger}{La seule référence fixe pour la suite est donc
    l'\textbf{énergie}.}
  \end{columns}
\end{frame}

% ---------- 6. Bilan d'energie ----------
\begin{frame}{Bilan d'énergie à l'impact}
  \begin{columns}[c]
  \column{0.42\textwidth}
    \centering
    \begin{tikzpicture}[scale=0.9]
      \draw[fill=danger!25](0,0) rectangle (1,3);
      
ode[font=\tiny] at (0.5,-0.35){$E_c$};
      \draw[-{Latex},thick](1.3,1.5)--(2.5,1.5) node[midway,above,font=\tiny]{$\eta$};
      \draw[fill=verre!30](2.8,0) rectangle (3.8,0.9);
      
ode[font=\tiny] at (3.3,-0.4){énergie};
      
ode[font=\tiny] at (3.3,-0.65){de flexion};
      
ode[gristik,font=\tiny] at (5.2,1.9){reste : rebond,};
      
ode[gristik,font=\tiny] at (5.2,1.6){son, téléphone};
    \end{tikzpicture}
  \column{0.58\textwidth}
    \footnotesize
    Une fraction $\eta$ de l'énergie de chute est stockée en \textbf{énergie
    élastique de flexion} (le reste part en rebond, son, déformation du
    téléphone) :
    \[ \tfrac12\,k\,x_{\max}^{2}=\eta\,E_c
       \;\Rightarrow\; x_{\max}=\sqrt{\frac{2\eta E_c}{k}} . \]
    $x_{\max}$ est la fléche maximale, $k$ la raideur de la poutre.
    \medskip\\
    \textcolor{danger}{$\eta$ est une hypothèse simplificatrice (on évite le
    calcul détaillé du contact).}
  \end{columns}
\end{frame}

% ---------- 7. Raideur et courbure ----------
\begin{frame}{Raideur et courbure de la poutre}
  \begin{columns}[c]
  \column{0.44\textwidth}
    \centering
    \begin{tikzpicture}[scale=0.95]
      \fill[gristik](0,0)--(-0.2,-0.3)--(0.2,-0.3)--cycle;
      \fill[gristik](3.4,0)--(3.2,-0.3)--(3.6,-0.3)--cycle;
      \draw[gray,dashed](0,0)--(3.4,0);
      \draw[verre,very thick](0,0) .. controls (1.7,-1.0) .. (3.4,0);
      \draw[-{Latex},danger,thick](1.7,0.9)--(1.7,-0.55) node[midway,right,font=\tiny]{$F$};
      \draw[{Latex}-{Latex}](3.65,0)--(3.65,-0.78) node[midway,right,font=\tiny]{$x_{\max}$};
      \draw[{Latex}-{Latex}](0,-1.05)--(3.4,-1.05) node[midway,below,font=\tiny]{$L$};
    \end{tikzpicture}
  \column{0.56\textwidth}
    \footnotesize
    Poutre sur deux appuis, charge centrale (résultat de RDM) :
    \[ x_{\max}=\frac{F L^{3}}{48\,E\,I}
       \;\Rightarrow\; k=\frac{F}{x_{\max}}=\frac{48\,E\,I}{L^{3}}
       =\frac{4Eb}{L^{3}}\,e^{3}\;\propto e^{3}. \]
    La courbure au centre se déduit de la fléche :
    \[ \kappa_{\max}=\frac{12\,x_{\max}}{L^{2}}\;\propto x_{\max}. \]
    \centering\textcolor{verre}{La raideur croît comme $e^{3}$ : doubler $e$ la
    multiplie par $8$.}
  \end{columns}
\end{frame}

% ---------- 8. Loi finale ----------
\begin{frame}{Loi finale : $\sigma_{\max}\propto 1/\sqrt{e}$}
  \begin{columns}[c]
  \column{0.55\textwidth}
    \footnotesize
    On combine les étapes précédentes :
    \[ \sigma_{\max}=E\,\frac{e}{2}\,\kappa_{\max}=\frac{6E\,e}{L^{2}}\,x_{\max},
       \qquad x_{\max}=\sqrt{\frac{2\eta E_c}{k}},\;\; k\propto e^{3}. \]
    Donc $x_{\max}\propto e^{-3/2}$ et :
    \[ \boxed{\;\sigma_{\max}\propto e\cdot e^{-3/2}=e^{-1/2}=\frac{1}{\sqrt{e}}\;} \]
    À énergie d'impact fixée, \textbf{épaissir réduit la contrainte}.
  \column{0.45\textwidth}
    \centering
    \begin{tikzpicture}[scale=1]
      \draw[-{Latex}](0,0)--(4.2,0) node[right,font=\tiny]{$e$};
      \draw[-{Latex}](0,0)--(0,3.2) node[above,font=\tiny]{$\sigma_{\max}$};
      \draw[verre,very thick,domain=0.4:4,samples=60,smooth]
        plot (\x,{1.7/sqrt(\x)});
      \draw[danger,dashed](0,1.25)--(4,1.25) node[right,font=\tiny,danger]{seuil};
      
ode[verre,font=\tiny] at (2.7,1.9){$1/\sqrt{e}$};
    \end{tikzpicture}
  \end{columns}
\end{frame}

% ---------- 9. Critere de rupture ----------
\begin{frame}{Critère de rupture du verre trempé}
  \begin{columns}[c]
  \column{0.46\textwidth}
    \centering
    \begin{tikzpicture}[scale=0.95]
      \draw[verre,thick,fill=verre!8](0,0) rectangle (3,1.6);
      \fill[ok!30](0,1.3) rectangle (3,1.6);
      \fill[ok!30](0,0) rectangle (3,0.3);
      
ode[ok,font=\tiny] at (1.5,1.45){compression de surface};
      
ode[font=\tiny] at (1.5,0.8){coeur};
      
ode[ok,font=\tiny] at (1.5,0.15){compression de surface};
      \draw[-{Latex},danger](3.3,0.8)--(3.9,0.8) node[right,font=\tiny]{traction appliquée};
    \end{tikzpicture}
  \column{0.54\textwidth}
    \footnotesize
    Le verre est \textbf{fragile} : il casse dès que la traction de surface
    dépasse un seuil.
    \medskip\\
    La \textbf{trempe} crée une \textbf{précontrainte de compression} en surface :
    la traction appliquée doit d'abord l'annuler $\Rightarrow$ seuil effectif
    élevé (plusieurs centaines de MPa).
    \medskip\\
    Rupture de nature \textbf{statistique} (micro-défauts), décrite par la loi de
    Weibull :
    \[ P_r(\sigma)=1-\exp\!\Big[-\big(\sigma/\sigma_\theta\big)^{m_W}\Big]. \]
  \end{columns}
\end{frame}

% ---------- 10. Localisation aux bords ----------
\begin{frame}{Pourquoi la rupture s'amorce aux bords}
  \begin{columns}[c]
  \column{0.46\textwidth}
    \centering
    \begin{tikzpicture}[scale=0.9]
      \fill[gristik](0,2.3)--(-0.22,2.0)--(0.22,2.0)--cycle;
      \fill[gristik](4,2.3)--(3.78,2.0)--(4.22,2.0)--cycle;
      \draw[verre,very thick](0,2.3) .. controls (2,1.85) .. (4,2.3);
      \draw[{Latex}-{Latex}](0,1.7)--(4,1.7) node[midway,below,font=\tiny]{$L$ grand (centre)};
      \fill[gristik](0,0.6)--(-0.18,0.35)--(0.18,0.35)--cycle;
      \fill[gristik](1.6,0.6)--(1.42,0.35)--(1.78,0.35)--cycle;
      \draw[danger,very thick](0,0.6) .. controls (0.8,0.2) .. (1.6,0.6);
      \draw[{Latex}-{Latex}](0,0.05)--(1.6,0.05) node[midway,below,font=\tiny]{$L$ petit (coin)};
    \end{tikzpicture}
  \column{0.54\textwidth}
    \footnotesize
    À énergie fixée, on montre que la contrainte dépend aussi de la portée :
    \[ \sigma_{\max}\propto \frac{1}{\sqrt{L}} . \]
    Au \textbf{coin}, la portée mobilisée $L$ est petite $\Rightarrow$ contrainte
    plus forte.
    \medskip\\
    S'ajoutent : une \textbf{concentration de contrainte} (facteur $K_c>1$) et une
    \textbf{trempe moins efficace} sur les tranches.
    \medskip\\
    \textcolor{danger}{C'est donc aux coins et aux bords que la fissure s'amorce.}
  \end{columns}
\end{frame}

% ---------- 11. Regard critique ----------
\begin{frame}{Regard critique sur le modèle}
  \footnotesize
  \begin{itemize}\itemsep4pt
    \item \textbf{Réduction à une poutre (1D)} : on néglige la flexion dans la
          seconde direction. On surestime la contrainte $\Rightarrow$ hypothèse
          \textbf{conservative}.
    \item \textbf{Fraction d'énergie $\eta$ constante} : hypothèse forte ; en
          réalité $\eta$ dépend du contact et de la hauteur.
    \item \textbf{Modèle énergétique vs impulsionnel} : selon l'idéalisation du
          choc, on obtient $\sigma_{\max}\propto e^{-1/2}$ ou $e^{-1}$. Dans les
          deux cas $\sigma_{\max}$ \emph{décroît} avec $e$ : la conclusion est
          \textbf{robuste}, seul l'exposant change.
    \item \textbf{Rupture statistique} : le seuil n'est pas unique ; on raisonne
          en \textbf{probabilité} (loi de Weibull), d'où la simulation de
          Monte-Carlo.
  \end{itemize}
\end{frame}

\end{document}

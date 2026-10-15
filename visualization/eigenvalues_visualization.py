r"""
# Stability Visualization

## Overview

This module provides an interactive visualization of the stability properties of the D1Q3 lattice Boltzmann scheme.

For given relaxation parameters \(t_v\) and \(t_w\), and a discrete velocity \(c\), the module computes and displays:

* the three eigenvalues of the amplification matrix in Fourier space;
* the unit circle in the complex plane for comparison;
* the three conditions obtained from the Cohn–Schur stability criterion;
* the current position of the relaxation parameters in the \((s_v,s_w)\)-parameter plane.

The visualization is interactive and allows the user to modify the relaxation parameters and the velocity through mouse and keyboard interactions.

## Main parameters

The relaxation parameters are expressed through

$$
s_v = 1+t_v,
\qquad
s_w = 1+t_w.
$$

The parameter \(c\) denotes the discrete velocity.

The parameter plane is displayed in terms of \((s_v,s_w)\), while the stability conditions are computed using \(t_v\), \(t_w\), \(c\), and the Fourier wave number \(\xi\).

## Constraint modes

The relaxation parameters can be modified according to four different constraint modes:

* **Free**: \(t_v\) and \(t_w\) are independent.
* **SRT**: \(t_w=t_v\).
* **dHG**: \(t_w=-t_v\).
* **limit**: \(t_w\) is set to the stability boundary

  $$
  t_w = \operatorname{borne\_stab}(t_v).
  $$

When switching modes, the current point is immediately projected onto the corresponding constraint in order to maintain visual consistency.

## Interactive controls

### Parameter plane

The point representing \((t_v,t_w)\) can be moved by clicking and dragging in the parameter-plane plot.

The active constraint mode determines how the coordinates of the point are modified.

### Keyboard

The relaxation parameters can also be modified using the arrow keys:

* **Left / Right**: decrease/increase \(t_v\);
* **Down / Up**: decrease/increase \(t_w\).

The increment is determined by `nbdigits_relax`.

### Speed selector

The value of \(c\) can be modified using the slider labelled \(c\).

The value is rounded according to `nbdigits_c`.

### Constraint mode selector

The radio buttons allow the user to switch between the four constraint modes:

`Free`, `SRT`, `dHG`, and `limit`.

## Main functions

### `borne_stab(t)`

Computes the stability boundary relating \(t_v\) and \(t_w\).

### `eigenvalues(tv, tw, c, xi)`

Computes the three eigenvalues of the D1Q3 amplification matrix for the Fourier wave numbers contained in `xi`.

The returned array has shape `(3, len(xi))`.

### `phi3(tv, tw, c, xi)`

Constructs the polynomial representation used in the Cohn–Schur stability analysis.

### `PtoPstar(P)`

Computes the reverse-conjugate polynomial associated with the polynomial representation `P`.

### `CohnSchur(P)`

Applies one step of the Cohn–Schur algorithm to the polynomial `P`.

### `derive(P)`

Computes the derivative of the polynomial represented by `P`.

### `condI(P)`

Computes the condition associated with the current step of the Cohn–Schur criterion.

### `lst_cond(tv, tw, c, xi)`

Computes the three stability conditions used in the visualization.

### `appliquer_contrainte(modecontrainte, tv, tw)`

Applies the selected constraint mode to the relaxation parameters and returns the corresponding values of \(t_v\) and \(t_w\).

## `ParametricInteractive`

The `ParametricInteractive` class manages the graphical interface and the interaction with the user.

### Main attributes

* `mode_contrainte`: currently selected constraint mode;
* `c`: discrete velocity;
* `tv`, `tw`: relaxation parameters;
* `xx`: array of Fourier wave numbers;
* `dragging`: indicates whether the parameter point is currently being dragged;
* `fig`: Matplotlib figure;
* `ax_curve`: eigenvalue plot;
* `ax_cond`: stability-condition plot;
* `ax_params`: relaxation-parameter plane;
* `radio_mode`: constraint-mode selector;
* `choixc`: velocity slider.

### Main methods

* `changement_mode(label)`: changes the active constraint mode and immediately applies the new constraint to the current parameters.
* `change_speed(val)`: updates the value of \(c\).
* `update_curve()`: recomputes and redraws the eigenvalues and stability conditions.
* `on_press(event)`: handles mouse-button press events in the parameter plane.
* `on_motion(event)`: updates the parameters while the mouse button is held down.
* `on_release(event)`: handles mouse-button release events.
* `on_key(event)`: updates the relaxation parameters in response to keyboard arrow events.

## Running the visualization

The application is launched by creating an instance of `ParametricInteractive`:

```python
app = ParametricInteractive()
plt.tight_layout()
plt.show()
```

The resulting Matplotlib window provides the complete interactive visualization.

"""
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.widgets import RadioButtons, Slider

mpl.rcParams['savefig.dpi'] = 150
plt.rcParams['figure.dpi'] = 80

UPS_primary = (0.38823529411764707, 0.0, 0.23529411764705882)
UPS_bleu = (0.0, 0.3058823529411765, 0.49019607843137253)
UPS_orange = (0.9333333333333333, 0.20392156862745098, 0.13725490196078433)
UPS_vert = (0.0, 0.5019607843137255, 0.47843137254901963)

nbdigits_relax, nbdigits_c = 2, 2
fmt_relax = f"{nbdigits_relax+2}.{nbdigits_relax}f"
fmt_c = f"%+{nbdigits_c+2}.{nbdigits_c}f"


def borne_stab(t):
    """Stability boundary"""
    a = t**2-6*t+1.
    b = 2*(1-t**2)
    return a/(b+np.sqrt(b*b+a*a))


def eigenvalues(tv, tw, c, xi):
    """compute the three eigenvalues of the D1Q3 operator in the Fourier space"""
    sv, sw = 1+tv, 1+tw
    d = c*c
    G0 = np.array(
        [
            [1-sw*d, sw*(1-d), sw*(1-d)],
            [
                (sw*d+sv*c),
                (2-sw*(1-d)-sv*(1-c)),
                (-sw*(1-d)+sv*(1+c))
            ],
            [
                (sw*d-sv*c),
                (-sw*(1-d)+sv*(1-c)),
                (2-sw*(1-d)-sv*(1+c))
            ]
        ]
    )
    vp_out = np.zeros((3, xi.size), dtype=complex)
    for k, xik in enumerate(xi):
        em = .5*np.exp(-1j*xik)
        ep = .5*np.exp(1j*xik)
        G = np.array([[1, 0, 0], [0, em, 0], [0, 0, ep]]) @ G0
        vp_out[:, k] = np.linalg.eig(G)[0]
    return vp_out


def phi3(tv, tw, c, xi):
    """phi3"""
    C = np.cos(xi)
    S = np.sin(xi)
    zzz = np.zeros(xi.shape)
    return [
        np.array([
            -tv*tw + zzz,
            tv*tw - (1-C)*tv*(1+tw)*c**2-C*(tv+tw),
            -1+c**2*(1+tw)*(1-C)+C*(tv+tw),
            1 + zzz
        ]), np.array([
            zzz,
            c*tw*(1+tv)*S,
            c*(1+tv)*S,
            zzz
        ])
    ]


def PtoPstar(P):
    """compute the reverse conjugate polynomial"""
    return [
        P[0][::-1], -P[1][::-1]
    ]


def CohnSchur(P):
    """Cohn-Schur algorithm"""
    Pstar = PtoPstar(P)
    c0r, c0i = P[0][0], P[1][0]
    cstar0r, cstar0i = Pstar[0][0], Pstar[1][0]
    return [
        cstar0r * P[0][1:] - cstar0i * P[1][1:]
        - c0r * Pstar[0][1:] + c0i * Pstar[1][1:],
        cstar0r * P[1][1:] + cstar0i * P[0][1:]
        - c0r * Pstar[1][1:] - c0i * Pstar[0][1:]
    ]


def derive(P):
    """derivation"""
    degree = P[0].shape[0] - 1
    ddd = np.arange(1, degree+1)[:, None]
    return [
        ddd * P[0][1:], ddd * P[1][1:]
    ]


def condI(P):
    """condition for P in I"""
    return P[0][-1]**2 + P[1][-1]**2 - P[0][0]**2 - P[1][0]**2


def lst_cond(tv, tw, c, xi):
    """conditions for the stability"""
    lst = []
    P = phi3(tv, tw, c, xi)
    cond = condI(P)
    lst.append(cond)
    for _ in range(2):
        Pd = derive(P)
        P = CohnSchur(P)
        idx_E = cond == 0
        P[0][:, idx_E] = Pd[0][:, idx_E]
        P[1][:, idx_E] = Pd[1][:, idx_E]
        cond = condI(P)
        lst.append(cond)
    return lst


def appliquer_contrainte(modecontrainte, tv, tw):
    """Compute tv and tw according to the selected mode and the click"""
    tv = min(1, max(-1, round(tv, nbdigits_relax)))
    tw = min(1, max(-1, round(tw, nbdigits_relax)))
    if modecontrainte == "SRT":
        return tv, tv
    elif modecontrainte == "dHG":
        return tv, -tv
    elif modecontrainte == "limit":
        return tv, borne_stab(tv)
    else:  # "Free" mode
        return tv, tw


class ParametricInteractive:
    """A class for the interactive plot"""
    def __init__(self):
        # Initial parameters
        self.mode_contrainte = "Free"  # Options : "Free", "SRT", "dHG", "limit"
        self.c = 0.0
        self.tv = 0.0
        self.tw = 0.0
        self.xx = np.linspace(-np.pi, np.pi, 1025)
        self.dragging = False  # click state

        # Axes creation
        self.fig = plt.figure(figsize=(15, 6))
        gs = self.fig.add_gridspec(3, 5, width_ratios=[1, 1, 1, 1, 1], height_ratios=[4, 2, 1])
        self.ax_curve = self.fig.add_subplot(gs[:, 0:2])
        self.ax_cond = self.fig.add_subplot(gs[:, 2:4])
        try:
            self.ax_cond.set_yscale('asinh', linear_width=1e-3, base=10)
        except ValueError:
            self.ax_cond.set_yscale('symlog', linthresh=0.001, linscale=0.2)
        self.ax_params = self.fig.add_subplot(gs[0, 4])
        self.ax_radio = self.fig.add_subplot(gs[1, 4])
        self.ax_radio.axis('off')
        self.ax_choixc = self.fig.add_subplot(gs[2, 4])
        self.ax_choixc.axis('off')

        # eigenvalues
        zzz = np.zeros(self.xx.shape)
        self.curve_0 = self.ax_curve.scatter(zzz, zzz, s=5, c=self.xx, cmap='hsv')
        self.curve_1 = self.ax_curve.scatter(zzz, zzz, s=5, c=self.xx, cmap='hsv')
        self.curve_2 = self.ax_curve.scatter(zzz, zzz, s=5, c=self.xx, cmap='hsv')
        self.ax_curve.plot(np.cos(self.xx), np.sin(self.xx), linewidth=0.5, color='black')
        self.ax_curve.set_xlim(-1, 1)
        self.ax_curve.set_ylim(-1, 1)
        self.ax_curve.set_xlabel("Real part")
        self.ax_curve.set_ylabel("Imaginary part")
        self.ax_curve.set_aspect('equal')
        self.ax_curve.set_title("Eigenvalues of the amplification matrix", fontsize=10, fontweight='bold')

        # stability condition
        self.ax_cond.axhline(0, linestyle='dotted', linewidth=0.5, color='black')
        self.cond_0, = self.ax_cond.plot(self.xx, np.nan*self.xx, linewidth=2, color='navy', label='cond_3')
        self.cond_1, = self.ax_cond.plot(self.xx, np.nan*self.xx, linewidth=2, color='orange', label='cond_2')
        self.cond_2, = self.ax_cond.plot(self.xx, np.nan*self.xx, linewidth=2, color='green', label='cond_1')
        self.ax_cond.set_xlim(self.xx[0], self.xx[-1])
        self.ax_cond.set_ylim(-1, 100)
        self.ax_cond.set_xlabel("Wave number")
        self.ax_cond.set_ylabel("The three conditions")
        self.ax_cond.set_title("The three conditions")
        self.ax_cond.legend()

        # parameters area
        self.ax_params.set_xlim(0, 2)
        self.ax_params.set_ylim(0, 2)
        self.ax_params.set_aspect('equal')
        titre = f"$s_v={self.tv+1:{fmt_relax}}$ $s_w={self.tw+1:{fmt_relax}}$"
        self.ax_params.set_title(titre, fontsize=10, fontweight='bold')
        self.ax_params.set_xlabel(r"$s_v$")
        self.ax_params.set_ylabel(r"$s_w$")
        self.point, = self.ax_params.plot(self.tv+1, self.tw+1, 'o', color=UPS_bleu)
        vecttv = np.linspace(-1, 1, 125)
        self.ax_params.plot(
            vecttv+1, vecttv+1, color='black', linewidth=1
        )
        self.ax_params.plot(
            vecttv+1, -vecttv+1, color='black', linewidth=1
        )
        self.ax_params.plot(
            vecttv+1, borne_stab(vecttv)+1,
            color=UPS_orange, linewidth=1
        )

        # Add the mode selector (RadioButtons)
        self.ax_radio.set_title("Constraint mode", fontsize=10, fontweight='bold')
        choix_modes = ('Free', 'SRT', 'dHG', 'limit')
        self.radio_mode = RadioButtons(self.ax_radio, choix_modes, active=0, activecolor='tomato')

        # Add the speed selector (Slider)
        self.ax_choixc.set_title("Speed selector", fontsize=10, fontweight='bold')
        self.choixc = Slider(
            ax=self.ax_choixc,
            label=r"$c$",
            valmin=-1.1,
            valmax=1.1,
            valinit=self.c,
            valfmt=fmt_c
        )

        # mouse events
        self.fig.canvas.mpl_connect('button_press_event', self.on_press)
        self.fig.canvas.mpl_connect('button_release_event', self.on_release)
        self.fig.canvas.mpl_connect('motion_notify_event', self.on_motion)
        self.fig.canvas.mpl_connect("key_press_event", self.on_key)
        self.radio_mode.on_clicked(self.changement_mode)
        self.choixc.on_changed(self.change_speed)

        # First plot
        self.update_curve()
    
    def changement_mode(self, label):
        self.mode_contrainte = label
        # When switching modes, we immediately apply
        # the new rule to the current point 
        # to ensure visual consistency.
        self.tv, self.tw = appliquer_contrainte(
                self.mode_contrainte, self.tv, self.tw
            )
        self.point.set_data([self.tv+1], [self.tw+1])
        self.update_curve()
    
    def change_speed(self, val):
        self.c = round(val, nbdigits_c)
        self.update_curve()

    def update_curve(self):
        """update the main curve according to (tv, tw)."""
        vp = eigenvalues(self.tv, self.tw, self.c, self.xx)
        self.curve_0.set_offsets(np.asarray([np.real(vp[0, :]), np.imag(vp[0, :])]).T)
        self.curve_1.set_offsets(np.asarray([np.real(vp[1, :]), np.imag(vp[1, :])]).T)
        self.curve_2.set_offsets(np.asarray([np.real(vp[2, :]), np.imag(vp[2, :])]).T)
        titre = f"$s_v={self.tv+1:{fmt_relax}}$ $s_w={self.tw+1:{fmt_relax}}$"
        self.ax_params.set_title(titre, fontsize=10, fontweight='bold')
        cond = lst_cond(self.tv, self.tw, self.c, self.xx)
        self.cond_0.set_ydata(cond[0])
        self.cond_1.set_ydata(cond[1])
        self.cond_2.set_ydata(cond[2])
        self.fig.canvas.draw_idle()

    def on_press(self, event):
        """mouse-button press event"""
        if event.inaxes == self.ax_params:
            self.dragging = True
            self.tv, self.tw = appliquer_contrainte(
                self.mode_contrainte, event.xdata-1, event.ydata-1
            )
            self.point.set_data([self.tv+1], [self.tw+1])
            self.update_curve()

    def on_motion(self, event):
        """mouse dragging"""
        if self.dragging and event.inaxes == self.ax_params:
            if event.xdata is not None and event.ydata is not None:
                self.tv, self.tw = appliquer_contrainte(
                    self.mode_contrainte, event.xdata-1, event.ydata-1
                )
                self.point.set_data([self.tv+1], [self.tw+1])
                self.update_curve()
                self.fig.canvas.draw_idle()

    def on_release(self, event):
        """mouse-button release event"""
        self.dragging = False
    
    def on_key(self, event):
        """arrow keys"""
        tv, tw = self.tv, self.tw
        if event.key == "left":
            tv = tv - 10**(-nbdigits_relax)
        elif event.key == "right":
            tv = tv + 10**(-nbdigits_relax)
        if event.key == "down":
            tw = tw - 10**(-nbdigits_relax)
        elif event.key == "up":
            tw = tw + 10**(-nbdigits_relax)
        self.tv, self.tw = appliquer_contrainte(
            self.mode_contrainte, tv, tw
        )
        self.point.set_data([self.tv+1], [self.tw+1])
        self.update_curve()
        self.fig.canvas.draw_idle()

# run application
app = ParametricInteractive()
plt.tight_layout()
plt.show()

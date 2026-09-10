# description.py
'''
Descriptions of pimc data types and formats.
'''

from dataclasses import dataclass


PARAM_NAMES = ("T", "V", "u", "t", "N", "n", "R", "L", "M")

# All estimators
EST_TYPE = (
    "estimator",
    "obdm",
    "pair",
    "pcycle",
    "super",
    "worm",
    "radial",
    "radwind",
    "radarea",
    "planedensity",
    "planewind",
    "planearea",
    "virial",
    "linedensity",
    "linepotential",
    "energy",
    "position",
    "ssf",
    "isf",
    "ssfq",
    "planeavedensity",
    "planeaveVext",
    "lineardensity",
    "number",
    "swap",
    "entpart"
)

# Cumulative estimators
CUM_EST_TYPE = ("position", "planeavedensity", "planeaveVext", "locsuper")

# Scalar estimators
SCALAR_EST_TYPE = ('estimator','energy','virial','super','swap','entpart')

# Vector estimators
VECTOR_EST_TYPE = (est for est in EST_TYPE if est not in SCALAR_EST_TYPE)

# -----------------------------------------------------------------------------
# Pure constants (no NDIM dependence)
# -----------------------------------------------------------------------------

LENGTH_T_UNIT = r"$[\mathrm{\AA}]$"
LENGTH_UNIT   = r"$\mathrm{\AA}$"

DENSITY_NAME = ("Line", "Area", "Volume")

PARAM_SHORT_NAME = {
    "T": "T",
    "V": "V",
    "M": "M",
    "u": r"$\mu$",
    "t": r"$\tau$",
    "N": "N",
    "n": r"$\rho$",
    "R": "R",
    "L": "L",
}

PARAM_FORMAT = {
    "T": r"%4.2f",
    "V": r"%3d",
    "M": r"%3d",
    "u": r"%+3.1f",
    "t": r"%5.3f",
    "N": "%3d",
    "n": r"%f",
    "R": r"% 4.1f",
    "L": r"%3d",
}

# -----------------------------------------------------------------------------
# NDIM-specific description object
# -----------------------------------------------------------------------------
@dataclass
class Desc:
    NDIM: int

    # "static" info (same for all NDIM)
    lengthTUnit: str = LENGTH_T_UNIT
    lengthUnit: str = LENGTH_UNIT
    densityName: tuple = DENSITY_NAME
    paramNames: tuple = PARAM_NAMES
    paramShortName: dict = None
    paramFormat: dict = None
    estType: tuple = EST_TYPE
    cumEstType: tuple = CUM_EST_TYPE

    # NDIM-dependent maps (filled in __post_init__)
    paramUnit: dict = None
    paramLongName: dict = None
    estimatorLongName: dict = None
    estimatorShortName: dict = None
    estimatorXLongName: dict = None

    def __post_init__(self):
        NDIM = self.NDIM
        if NDIM not in (1, 2, 3):
            raise ValueError("NDIM must be 1, 2, or 3")

        # bind shared dicts (you can copy() if you prefer them isolated too)
        self.paramShortName = PARAM_SHORT_NAME
        self.paramFormat = PARAM_FORMAT

        self.paramUnit = {
            "T": "K",
            "M": "",
            "V": rf"$\mathrm{{\AA^{{{NDIM}}}}}$",
            "u": "K",
            "t": r"$K^{-1}$",
            "N": "",
            "n": rf"$\mathrm{{\AA}}^{{-{NDIM}}}$",
            "R": self.lengthUnit,
            "L": self.lengthUnit,
        }

        self.paramLongName = {
            "T": "Temperature  [K]",
            "V": rf"Volume  $[\mathrm{{\AA^{{{NDIM}}}}}]$",
            "u": "Chemical Potential  [K]",
            "t": "Imaginary Time Step  [1/K]",
            "N": "Number of Particles",
            "n": rf"{self.densityName[NDIM-1]} Density  $[\mathrm{{\AA}}^{{-{NDIM}}}]$",
            "R": f"Pore Radius  {self.lengthTUnit} ",
            "L": f"Length {self.lengthTUnit}",
            "W": "Virial Window [1/K]",
            "M": "Update Length",
            "D": r"CoM Delta [$\AA$]",
        }

        self.estimatorLongName = {
            "K": "Kinetic Energy [K]",
            "V": "Potential Energy [K]",
            "V_op": "Potential Energy [K]",
            "Vcv": "Potential Energy [K]",
            "E": "Energy [K]",
            "Ecv": "Energy [K]",
            "Eth": "Energy [K]",
            "E_mu": r"$E - \mu N$",
            "K/N": "Kinetic Energy per Particle [K]",
            "V/N": "Potential Energy per Particle [K]",
            "E/N": "Energy per Particle [K]",
            "Ecv/N": "Energy per Particle [K]",
            "Cv": "Specific Heat [K]",
            "P": rf"Pressure $[\mathrm{{K}}\mathrm{{\AA}}^{{-{NDIM}}}]$",
            "N": "Number of Particles",
            "N^2": r"(Number of Particles)$^2$",
            "density": rf"{self.densityName[NDIM-1]} Density $[\mathrm{{\AA}}^{{-{NDIM}}}]$",
            "diagonal": "Diagonal Fraction",
            "kappa": r"$\rho^2 \kappa [units]$",
            "pair": "Pair Correlation Function [units]",
            "radial": r"Radial Density $[\mathrm{\AA}^{-3}]$",
            "number": "Number Distribution",
            "planedensity": r"Density $[\mathrm{\AA}^{-3}]$",
            "obdm": "One Body Density Matrix",
            "rho_s/rho": "Superfluid Fraction",
            "Area_rho_s": "Area Superfluid Fraction",
            "rho_s/rho|Z": r"$\rho_s/\rho_0$",
            "radwind": r"Radial Winding Superfliud Density $[\mathrm{\AA}^{-3}]$",
            "radarea": r"Radial Area Superfliud Density $[\mathrm{\AA}^{-3}]$",
            "x^2": r"$ \langle x^2 \rangle [\AA^2]$",
        }

        self.estimatorShortName = {
            "K": "K [K]",
            "V": "V [K]",
            "E": "E [K]",
            "Eth": "Eth [K]",
            "Ecv": "Ecv [K]",
            "E_mu": r"$E - \mu N$",
            "K/N": "K/N [K]",
            "V/N": "V/N [K]",
            "E/N": "E/N [K]",
            "N": "N",
            "N^2": r"N$^2$",
            "P": rf"$P [\mathrm{{K}}\mathrm{{\AA}}^{{-{NDIM}}}]$",
            "density": rf"$\rho [\mathrm{{\AA}}^{{-{NDIM}}}]$",
            "diagonal": "D",
            "kappa": r"$\rho^2 \kappa [units]$",
            "pair": "g(r) [units]",
            "radial": r"Radial Density $[\mathrm{\AA}^{-3}]$",
            "number": "Number Distribution",
            "obdm": "One Body Density Matrix",
            "rho_s/rho": r"$\rho_s/\rho$",
            "rho_s/rho|Z": r"$\rho_s/\rho_0$",
            "Area_rho_s": "Area Superfluid Fraction",
        }

        self.estimatorXLongName = {
            "number": "Number of Particles",
            "pair": f"r  {self.lengthTUnit}",
            "obdm": f"r  {self.lengthTUnit}",
            "planedensity": "Gridbox Number",
            "radial": f"r  {self.lengthTUnit}",
            "radwind": f"r  {self.lengthTUnit}",
            "radarea": f"r  {self.lengthTUnit}",
        }

        # Note 2025-10-06: ssf and ssfq are duplicated for backwards compatibility
        self.axisLabel = {
            'obdm':['r [Å]','n(r)'], 
            'pair':['r [Å]','g(r)'], 
            'radial':['r [Å]','ρ(r)'], 
            'number':['N','P(N)'],
            'radwind':['r [Å]','ρₛ(r)'],
            'radarea':['r [Å]','ρₛ(r)'],
            'planewind':['n','ρₛ(x,y)'],
            'planearea':['n','ρₛ(x,y)'],
            'planedensity':['n','ρ(x,y)'], 
            'linedensity':['r [Å]','ρ1d(r)'],
            'linepotential':['r [Å]','V1d(r)'], 
            'ssf':['q_index', 'S(q)'],
            'ssfq':['q_index', 'S(q)'], 
            'planeavedensity':['n','ρ(x,y)'],
            'lineardensity':['r [Å]','ρ(r)'], 
            'isf':['τ [1/K]','F(q,τ)'],
            'planeaveVext':['n','Vext(x,y)'],
            'pcycle':['m','P(m)'],
            'position':['n','ρ'],
            'estimator':'',
            'energy':'',
            'super':'',
            'virial':'',
            'swap':'',
            'entpart':''
        }

# --------------------------------------------------------------------------------
def get_desc(NDIM=3):
    """Return a fresh Desc instance each call (no caching)."""
    return Desc(NDIM)


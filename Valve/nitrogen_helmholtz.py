#################################################################################
# IDAES PSE Examples Repository
# https://github.com/virajdesai0309/IDAES-PSE-Examples
#
# This file adds nitrogen to the IDAES general Helmholtz equation of state
# framework, which ships without nitrogen.
#
# The equation of state is the reference equation used by CoolProp (and so by
# DWSIM's CoolProp property package):
#
# Span, R., Lemmon, E.W., Jacobsen, R.T, Wagner, W., Yokozeki, A. (2000).
#     "A Reference Equation of State for the Thermodynamic Properties of
#     Nitrogen for Temperatures from 63.151 to 1000 K and Pressures to
#     2200 MPa." Journal of Physical and Chemical Reference Data, 29, 1361.
#
# The coefficients are read from CoolProp's fluid database instead of being
# typed by hand. Run this script once to (re)generate the files in
# helmholtz_data/; it requires the CoolProp package.
#################################################################################
"""Generate IDAES Helmholtz parameter and expression files for nitrogen."""

__author__ = "Viraj Desai"

import json
import os

import CoolProp.CoolProp as CP
from idaes.models.properties.general_helmholtz.helmholtz_parameters import (
    WriteParameters,
)

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "helmholtz_data")


def numbered(values, start=1):
    """IDAES parameter files index coefficients from 1 (or another start)."""
    return {i: v for i, v in enumerate(values, start=start)}


def nitrogen_parameters():
    """Build the IDAES Helmholtz parameter dictionary from CoolProp's EOS data."""
    fluid = json.loads(CP.get_fluid_param_string("Nitrogen", "JSON"))[0]
    eos = fluid["EOS"][0]
    reducing = eos["STATES"]["reducing"]

    mw = eos["molar_mass"] * 1000  # g/mol
    Tc = reducing["T"]  # K
    rhoc = reducing["rhomolar"] * eos["molar_mass"]  # kg/m3
    Tt = eos["Ttriple"]  # K

    def sat(prop, q):
        return CP.PropsSI(prop, "T", Tt, "Q", q, "Nitrogen")

    # Ideal part: ln(delta) + a1 + a2*tau + a*ln(tau) + sum(n*tau^t)
    #             + sum(n*ln(1 - exp(-v/Tc*tau)))   (IDAES ideal type 02)
    alpha0 = {term["type"]: term for term in eos["alpha0"]}
    lead = alpha0["IdealGasHelmholtzLead"]
    power0 = alpha0["IdealGasHelmholtzPower"]
    planck = alpha0["IdealGasHelmholtzPlanckEinsteinFunctionT"]
    n0 = [lead["a1"], lead["a2"], alpha0["IdealGasHelmholtzLogTau"]["a"]]
    n0 += power0["n"] + planck["n"]
    g0 = power0["t"] + [v / planck["Tcrit"] for v in planck["v"]]
    last_power = 3 + len(power0["n"])

    # Residual part: polynomial, exponential exp(-delta^l) and Gaussian terms
    # (IDAES residual type 02). CoolProp stores polynomial terms as l = 0.
    alphar = {term["type"]: term for term in eos["alphar"]}
    power = alphar["ResidualHelmholtzPower"]
    gauss = alphar["ResidualHelmholtzGaussian"]
    terms = sorted(
        zip(power["n"], power["d"], power["t"], power["l"]), key=lambda x: x[3] > 0
    )
    n = [x[0] for x in terms] + gauss["n"]
    d = [x[1] for x in terms] + gauss["d"]
    t = [x[2] for x in terms] + gauss["t"]
    n_poly = sum(1 for x in terms if x[3] == 0)
    n_exp = len(terms)
    first_gauss = n_exp + 1

    # Approximate saturated densities (initial guesses for phase equilibrium)
    rho_l = fluid["ANCILLARIES"]["rhoL"]  # rho/rhoc = 1 + sum(n*theta^t)
    rho_v = fluid["ANCILLARIES"]["rhoV"]  # rho/rhoc = exp(Tc/T*sum(n*theta^t))

    return {
        "comp": "nitrogen",
        "basic": {
            "R": eos["gas_constant"] / mw,  # kJ/kg/K
            "MW": mw,
            "T_star": Tc,
            "rho_star": rhoc,
            "Tc": Tc,
            "rhoc": rhoc,
            "Pc": reducing["p"] / 1000,  # kPa
            "Tt": Tt,
            "Pt": sat("P", 0) / 1000,  # kPa
            "rhot_l": sat("D", 0),
            "rhot_v": sat("D", 1),
            "P_min": 1e-9,
            "P_max": 2.2e6,  # kPa
            "rho_max": 1300.0,
            "T_min": Tt,
            "T_max": 1000.0,
        },
        "eos": {
            "reference": ["Span, R., Lemmon, E.W., Jacobsen, R.T, Wagner, W., Yokozeki, A. (2000), J. Phys. Chem. Ref. Data, 29, 1361."],
            "phi_ideal_type": 2,
            "last_term_ideal": [last_power, len(n0)],
            "n0": numbered(n0),
            "g0": numbered(g0, start=4),
            "phi_residual_type": 2,
            "last_term_residual": [n_poly, n_exp, len(n)],
            "n": numbered(n),
            "d": numbered(d),
            "t": numbered(t),
            "c": {i + 1: terms[i][3] for i in range(n_poly, n_exp)},
            "a": numbered(gauss["eta"], start=first_gauss),
            "e": numbered(gauss["epsilon"], start=first_gauss),
            "b": numbered(gauss["beta"], start=first_gauss),
            "g": numbered(gauss["gamma"], start=first_gauss),
        },
        "aux": {
            "delta_l_sat_approx": {"type": 1, "c": 1.0, "n": numbered(rho_l["n"]), "t": numbered(rho_l["t"])},
            "delta_v_sat_approx": {"type": 3, "c": 1.0, "n": numbered(rho_v["n"]), "t": numbered(rho_v["t"])},
        },
    }


def main():
    os.makedirs(DATA_DIR, exist_ok=True)
    parameters = nitrogen_parameters()
    cwd = os.getcwd()
    os.chdir(DATA_DIR)  # WriteParameters writes its files to the working directory
    try:
        with open("nitrogen.json", "w") as f:
            json.dump(parameters, f, indent=4)
        WriteParameters(parameters="nitrogen.json").write()
    finally:
        os.chdir(cwd)


if __name__ == "__main__":
    main()

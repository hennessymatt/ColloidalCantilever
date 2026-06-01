"""
A class containing all of the parameter values needed to simulate
the model
"""

class Params():
    def __init__(self, RH = 50, G_0 = 1e5, t_g = 700):
        """
        The constructor.  Takes as optional inputs:

        RH = the value of the relative humidity
        G_0 = a reference value of the gel stiffness for non-dimensionalisation
        t_g = the salt-induced gelation time with 1 M salinity
        """

        #---------------------------------------------
        # Simulation end time (normalised by t_d)
        #---------------------------------------------
        self.t_end = 0.18

        #---------------------------------------------
        # constants
        #---------------------------------------------
        self.g = 9.81

        #---------------------------------------------
        # dimensional quantities for the beam
        #---------------------------------------------
        self.h_b = 70e-6                        # thickness
        self.L = 22e-3                          # length
        self.W = 8e-3                           # width
        self.E_b = 200e9                        # Young's modulus
        self.nu_b = 0.3                         # Poisson's ratio
        self.B = self.E_b * self.h_b**3 / 12    # Bending modulus

        # density (material assumed to be steel)
        self.rho_b = 8000

        # aspect ratio
        self.delta = self.h_b / self.L

        #---------------------------------------------
        # dimensional quantities for the film
        #---------------------------------------------
        self.V_0 = 30e-6 / 1000                 # Volume (30 uL)
        
        self.h_f = 289e-6                       # Film thickness
        self.V_e = 2e-8 * (1 - RH / 100) * 2    # Evap rate
        self.t_d = self.h_f / self.V_e          # Drying time

        # gel time (= 700 s for 1M salt soln, e.g. the default)
        self.t_g = t_g

        # Shear modulus
        self.G_0 = G_0          # characteristic value
        self.G_p = 6.11e8       # packing stiffness (obtained from fitting)

        # Poisson's ratio of the packing
        self.nu_p = 0.2

        # solid fractions
        self.phi_g = 0.32   # gel point
        self.phi_0 = 0.10   # initial
        self.phi_p = 0.64   # random close packing fraction

        # fluid fractions
        self.phi_f_0 = 1 - self.phi_0   # initial fluid
        self.phi_f_inf = 0.0            # air

        # Particle radius
        self.a = 6e-9

        # Permeability at the gel point
        self.k_0 = 1 / 45 * (self.a)**2 * (1 - self.phi_g)**3 / self.phi_g**2

        # Viscosities
        self.mu_f = 1e-3            # interstitial fludid (water)
        self.mu_m = 8e-3            # mixture (water + nanoparticles)

        # Surface tension
        self.gamma = 64e-3

        # Densities
        self.rho_w = 1000       # water
        self.rho_p = 2216       # SiO2
        self.rho_m = lambda phi_f: self.rho_w * phi_f + self.rho_p * (1 - phi_f) # mixture
        self.rho_0 = self.rho_m(self.phi_f_0)   # initial mixture

        # particle diffusivity
        self.D_p = 3.6e-11 # Stokes-Einstein (checked)

        #----------------------------------------------------------------
        # Fitted params depending on RH
        #----------------------------------------------------------------

        if RH <= 50:
            self.m = 46
            self.G_c = 0.24
        elif RH == 60:
            self.m = 39
            self.G_c = 0.22
        elif RH > 60: # 70 and 80
            self.m = 38
            self.G_c = 0.22
        else:
            raise Exception(f'Need values for m and G_c when RH = {RH}')
            # 1 + 1

        # Flag to include stress in colloidal gel
        self.incremental = 1

        # Pore blockage factor (0 if unblocked)
        self.blockage = 0

        # Flag to apply contact-stress model
        self.contact_stress = True

        #----------------------------------------------------------------
        # Non-dimensional quantities
        #----------------------------------------------------------------
        self.eps = self.h_f / (self.L / 2)
        self.delta = self.h_b / self.L
        self.alpha = self.h_f / self.h_b

        # Bending moment
        self.E = self.G_0 * self.h_f * self.L**2 / self.B

        # Pe numbers
        self.Pe = self.mu_f * self.V_e * (self.L / 2)**2 / self.k_0 / self.h_f / self.G_0
        self.Pe_p = self.V_e * (self.L / 2) / self.eps / self.D_p

        # film weight
        self.G = self.rho_0 * self.g * self.h_f * self.L**4 / self.B / self.h_b

        # Actual bond number
        self.Bond = self.rho_0 * self.g * (self.L / 2)**2 / self.gamma

        # Bond-ish number
        self.Bo = self.eps**4 * self.rho_0 * self.g * (self.L / 2)**2 / self.mu_m / self.V_e

        # Capillary number
        self.Ca = self.mu_m * self.V_e / self.eps**4 / self.gamma

        # Normalised shear modulus of packed region
        self.E_p = self.G_p / self.G_0

        #----------------------------------------------------------------
        # flags
        #----------------------------------------------------------------
        self.drying_gelation = True

    def to_dict(self):
        """
        Converts the object into a dict.  Excludes
        functions
        """

        new_dict = {
            i:self.__dict__[i] for i in self.__dict__ if (
                type(self.__dict__[i]) is float
                or
                type(self.__dict__[i]) is int
                or
                type(self.__dict__[i]) is bool
            )
        }

        return new_dict
    

    def convert_dict(self, d):
        """
        Converts a loaded param set from a json file into
        a Param object
        """

        for i in d:
            self.__setattr__(i, d[i])

        # update density
        self.rho_m = lambda phi_f: self.rho_w * phi_f + self.rho_p * (1 - phi_f)
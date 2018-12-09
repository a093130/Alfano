# -*- coding: utf-8 -*-
"""
Created on Sun Apr 30 12:25:44 2017

Description: Class intended to perform a Sims Flanagen trajectory optimization for a LEO to GEO
    96kW SEP spacecraft using Busek Clustered BHT-8000 Hall Effect Thrusters
    
@author: Colin Helms (chelms@socal.rr.com), Francesco Biscani (bluescarni@gmail.com)
    R0 - Copied from PyKEP example1 and modified for conda pygmo installation (no caps in module name)
"""
try:
    from pygmo.core import problem as PyGMO_problem

    """
    This, is a low-thrust optimization solved using one of the available PagMO solvers. 
    The problem is a non-linear constrained problem that uses the Sims-Flanagan transcription to model the low-thrust trajectory. 
    """
    class NLP_earth_GEO(PyGMO_problem):

        """
        This constructs a PaGMO.problem object that represents a low-thrust transfer between LEO and GEO. The decision vector
        contains [t0,T,mf,Vx,Vy,Vz,[throttles]] in the following units: [mjd2000, days, kg, m/s,m/s,m/s, [non-dimensional]]

       This class represents a box-bounded, multiobjective, mixed-integer, constrained optimisation problem defined by:
        - a global dimension, i.e., the number of dimensions of the global search space,
        - the dimension of the integral (or combinatorial) part of the problem,
        - lower and upper bounds of the global search space,
        - the number of inequality constraints (never exceeding the total number of constraints),
        - a constraint computation function,
        - a constraints tolerance vector.
        
        All dimensions are invariant in the life cycle of a problem object.

        """

        def __init__(self, mass=8000, Tmax=0.05, Isp=2020, Vinf=7.612, nseg=20):
            # First we call the constructor for the base PyGMO problem
            # (dim, integer dim, number of obj, number of con, number of inequality con, tolerance on con violation)
            super(NLP_earth_GEO, self).__init__(6 + nseg * 3, 0, 1, 8 + nseg, nseg + 1, 1e-4)

            # We then define some data members (we use the double underscore to
            # indicate they are private)
            from PyKEP import MU_EARTH
            from PyKEP.planet import jpl_lp
            from PyKEP.sims_flanagan import spacecraft, leg
            
            self.__earth = jpl_lp('earth')
            #self.__mars = jpl_lp('mars')
            self.__sc = spacecraft(mass, Tmax, Isp)
            self.__Vinf = Vinf * 1000
            self.__leg = leg()
            self.__leg.set_mu(MU_EARTH)
            self.__leg.set_spacecraft(self.__sc)
            self.__nseg = nseg
            """
             The bounds of the problem are allowed to vary over the whole range of double-precision values for continuous optimisation,
             while for combinatorial optimisation the bounds must be in the [-32767,32767] range (corresponding to the INT_MIN and INT_MAX
             constants defined in the C++ standard "climits" header file). 
             All bounds setting functions will make sure that the following conditions are respected:
                 - lower bounds are not greater than upper bounds,
                 - the bounds of the integer part of the problem are integer and they are within the allowed range.
                 
            void set_bounds(const double (&v1)[N], const double (&v2)[N])
            Set lower and upper bounds to the content of the raw arrays v1 and v2. Will fail if N is different
            from the global size of the problem or if at least one lower bound is greater than the corresponding upper bound.
                @param[in] v1 lower bounds for the problem.
                @param[in] v2 upper bounds for the problem.
           """
            self.set_bounds([2480, 2400, self.__sc.mass / 10, -self.__Vinf, -self.__Vinf, -self.__Vinf] + [-1]
                            * 3 * nseg, [2490, 2500, self.__sc.mass, self.__Vinf, self.__Vinf, self.__Vinf] + [1] * 3 * nseg)

        # This is the objective function
        def _objfun_impl(self, x):
            return (-x[2],)

        # And these are the constraints
        def _compute_constraints_impl(self, x):
            from PyKEP import epoch, AU, EARTH_VELOCITY
            from PyKEP.sims_flanagan import leg, sc_state

            start = epoch(x[0])
            end = epoch(x[0] + x[1])

            r, v = self.__earth.eph(start)
            v = [a + b for a, b in zip(v, x[3:6])]
            x0 = sc_state(r, v, self.__sc.mass)

            r, v = self.__earth.eph(end)
            xe = sc_state(r, v, x[2])

            self.__leg.set(start, x0, x[-3 * self.__nseg:], end, xe)
            v_inf_con = (x[3] * x[3] + x[4] * x[4] + x[5] * x[5] -
                         self.__Vinf * self.__Vinf) / (EARTH_VELOCITY * EARTH_VELOCITY)
            retval = list(self.__leg.mismatch_constraints() +
                          self.__leg.throttles_constraints()) + [v_inf_con]

            # We then scale all constraints to non-dimensional values
            retval[0] /= AU
            retval[1] /= AU
            retval[2] /= AU
            retval[3] /= EARTH_VELOCITY
            retval[4] /= EARTH_VELOCITY
            retval[5] /= EARTH_VELOCITY
            retval[6] /= self.__sc.mass
            return retval

        # This transforms the leg into a high fidelity one
        def high_fidelity(self, boolean):
            self.__leg.high_fidelity = boolean

        # And this helps to visualize the trajectory
        def plot(self, x):
            import matplotlib as mpl
            from mpl_toolkits.mplot3d import Axes3D
            import matplotlib.pyplot as plt
            from PyKEP import epoch, AU
            from PyKEP.sims_flanagan import sc_state
            from PyKEP.orbit_plots import plot_planet, plot_sf_leg

            # Making sure the leg corresponds to the requested chromosome
            self._compute_constraints_impl(x)
            start = epoch(x[0])
            end = epoch(x[0] + x[1])

            # Plotting commands
            fig = plt.figure()
            axis = fig.gca(projection='3d')
            # The Sun
            axis.scatter([0], [0], [0], color='y')
            # The leg
            plot_sf_leg(self.__leg, units=AU, N=10, ax=axis)
            # The planets
            plot_planet(
                self.__earth, start, units=AU, legend=True, color=(0.8, 0.8, 1), ax = axis)
            plot_planet(
                self.__mars, end, units=AU, legend=True, color=(0.8, 0.8, 1), ax = axis)
            plt.show()

    def run_PlotTraj(n_restarts=5):
        from PyGMO import algorithm, island
        
        prob = mga_lt_earth_GEO(nseg=15)
        prob.high_fidelity(True)
        algo = algorithm.scipy_slsqp(max_iter=500, acc=1e-5)
        # algo = algorithm.snopt(major_iter=1000, opt_tol=1e-6, feas_tol=1e-11)
        algo2 = algorithm.mbh(algo, n_restarts, 0.05)
        algo2.screen_output = True
        isl = island(algo2, prob, 1)
        
        print("Running Monotonic Basin Hopping .... this will take a while.")
        
        isl.evolve(1)
        isl.join()
        
        print("Is the solution found a feasible trajectory? " +
              str(prob.feasibility_x(isl.population.champion.x)))
        
        prob.plot(isl.population.champion.x)

except:
    print("Could not import pygmo - is required for mga_lt_earth_GEO.")


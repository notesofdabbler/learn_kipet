#
# https://github.com/kwmcbride/kipet_examples/blob/master/examples/example_1/Ex_1_ode_sim.py
#

# Standard library imports
import sys

# Kipet library imports
import kipet

import matplotlib.pyplot as plt

# Create the ReactionModel instance
r1 = kipet.ReactionModel('reaction-1')

# Change the desired time basis here (if different from default)
r1.unit_base.time = 's'

# Add the model parameters
k1 = r1.parameter('k1', value=2, units='1/s')
k2 = r1.parameter('k2', value=0.2, units='1/s')

# Declare the components and give the initial values
A = r1.component('A', value=1.0, units='M')
B = r1.component('B', value=0.5, units='M')
C = r1.component('C', value=0.0, units='M')
    
# Input the reactions as expressions
rA = r1.add_reaction('rA', k1*A)
rB = r1.add_reaction('rB', k2*B)

# Input the ODEs
r1.add_ode('A', -rA )
r1.add_ode('B', rA - rB )
r1.add_ode('C', rB )

# Option to check the units of your models
r1.check_model_units(display=True)

# Dosing requires the component, the time, the concentration, and volume
# The concentrations and volumes are converted ot the unit_base values
r1.add_dosing_point('A', time=3, conc=(2, 'M'), vol=(200, 'mL'))
r1.add_dosing_point('A', time=5, conc=(3, 'M'), vol=(0.33, 'L'))

# Simulations require a time span
r1.set_time(10)

# Change some of the default settings
r1.settings.collocation.ncp = 3
r1.settings.collocation.nfe = 50
r1.settings.solver.linear_solver = 'ma27'
r1.settings.simulator.tee = True

r1.simulate()
    
# Create plots
r1.report()

fig, ax = plt.subplots()
ax.plot(r1.results_dict['simulator'].Z['A'])
plt.savefig('tmp.png')
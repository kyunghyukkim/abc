__author__ = 'Keagan Moo'

import MooScripts as moo
import KimScripts as kim
import pandas as pd

param_input = moo.input_to_df("10 0 100 0.4 2 100 0.1 0.02 10", 3, 4)
particles = kim.initial_particles(param_input, 10)
result = moo.PerturbationKernel(particles, 10)
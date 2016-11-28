__author__ = 'Keagan Moo'

import MooScripts as moo
import KimScripts as kim
import pandas as pd

param_input = moo.input_to_df("10 0 100 0.4 2 100 0.1 0.02 10", 3, 4)
Thet = kim.initial_particles(param_input, 1)
moo.priorProb(Thet, [0.4, 1, 200, 0.01], [0.4, 1, 200, 0.01])

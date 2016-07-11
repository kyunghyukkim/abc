import MooScripts as moo
import KimScripts as kim

# math model is defined in c codes. 
# three species: init copy numbers = [10, 0, 100]
# 4 parameters: 
# example: species 1 synthesis rate, deg rate = 0.4, 2
# Third last param: time interval
# Second last param: Number of data points.
# the last param: not necessary at this point. 

param_input = '10 0  100  0.4    2   100   0.1  0.02   10000   240000'

# input: a threshold epsilon
epsilon = 0.1
# output: a weighted sample of particles from p(theta|x)
particles = #keagan, fill this out.
t = 0
N_sample = 1
for i in range(1:N_sample):
	t++
 	epsilon = epsilon_next(epsilon)
#	determine thhe parameters of the perturbation kernel K(.|.)
############
##  This is justa rough draft. 
##  Keagan, can you make the rest of the code?
############
#	i=1
#	repeat
#		if t=1:
			sample theta_tilde from pi(theta)
		else:
			sample theta from the rpevious population with weights
			sample theta_tilde from K(.|theta) and such that pi(theta_tilde) > 0 
		sample y from f(.|theta_tilde)
		if moo.euclidd(y,x) <= epsilon:
			theta[i][t] = theat_tilde
			y[i][t] = y
# stochastic simulation
# something like this...
tseries =  kim.ssa(param_input)


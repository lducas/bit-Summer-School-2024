##########
# Exercise 1
# 
# Generate the systematic basis of random SIS lattices
# given parameters n, m and q.
# Make sure that the q-vectors are put first !
##########

##########
# Exercise 2
# Practical approximation factor
#
# Setting a = sqrt{4/3} + epsilon
# The Theory of LLL predicts that ||b_1|| < a^(n-1)/2 * det(L)^1/n
# Equivalently,  1/n log ||b_1|| - 1/n log(det(L)) <= 2 log(a) 
# This quantity is called the Root-Hermite-Factor. 
# Plot thew root hermite factor in practice for q = 997, and
# m = n, for n in {5, 10, 15, 20, 25, 30}, and compare it to
# the theoretical upper bound
#
##########

##########
# Exercise 3
#
# implement a formula for predicting how short of a vector LLL
# finds for a SIS lattice of parameter n,m,q; using the practical
# RHF rather than the theoretical RHF bound. 
#
# Also implement a function that choose the optimal m for a given
# n and q.
##########

##########
# Exercise 4
# Compare this prediction to practice for SIS parameters: q = 997,
# n = 1, and growning m in {5, 10, 15, ...}. Also show the optimal
# m on this plot.
#
# Note that the predictions seems to fail for large m. Try to
# find a explanation.
# 
# HINT: look at the profile of the LLL basis for those large m.
###########



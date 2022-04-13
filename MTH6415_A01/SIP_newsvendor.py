"""
Newsvendor problem with discrete uniform distribution.
"""


import pulp


if __name__ == "__main__":
    # CREATE LP MINIMIZATION PROBLEM
    mdl = pulp.LpProblem('Newsvendor', pulp.LpMinimize)

    # PROBLEM CONSTANTS
    c = 0.75  # Buying cost
    q = 2.25  # Selling price
    r = 0.25  # Return price
    LB = 60   # Uniform distribution lower bound
    UB = 100  # Uniform distribution upper bound

    # CREATE PROBLEM VARIABLES
    x = []  # Decision variable for number of papers bought
    y = []  # Number of papers sold under scenario s
    w = []  # Number of papers returned under scenario s
    x = pulp.LpVariable("x", lowBound=LB, upBound=UB, cat=pulp.const.LpInteger)   # Variable papers to buy, LB <= x <= UB
    for s in range(1, (UB + 2) - LB):
        # Create a variable y >= 0 and y <= s
        y.append(pulp.LpVariable("y_" + str(s), lowBound=0, upBound=(s + (LB - 1)), cat=pulp.const.LpInteger))
        # Create a variable w >= 0 and w <= upper bound
        w.append(pulp.LpVariable("w_" + str(s), lowBound=0, upBound=UB, cat=pulp.const.LpInteger))

    # OBJECTIVE FUNCTION
    revenue = 0.0
    for s in range(0, (UB + 1) - LB):
        revenue += ((1.0 / ((UB + 1) - LB)) * ((q * y[s]) + (r * w[s])))
    mdl += (c * x) - revenue

    # CONSTRAINTS
    for s in range(0, (UB + 1) - LB):
        mdl += x - y[s] - w[s] == 0

    # DISPLAY THE PROBLEM
    print(mdl)

    # SOLVE
    status = mdl.solve()
    print(pulp.LpStatus[status])   # The solution status

    # PRINTING THE FINAL SOLUTION
    print('x = ', pulp.value(x))
    for s in range(1, (UB + 2) - LB):
        print('y_' + str(s) + ' = ', pulp.value(y[s - 1]), '\t', 'w_' + str(s) + ' = ', pulp.value(w[s - 1]))
    print('Objective function value = ', pulp.value(mdl.objective))

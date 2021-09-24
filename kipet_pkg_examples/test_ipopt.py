from pyomo.environ import *

M = ConcreteModel()
M.x = Var()
M.y = Var()
expr1 = (M.x-1)**2 + 100*(M.y-M.x**2)**2

M.o  = Objective(expr=expr1)

solver = SolverFactory('ipopt')
res = solver.solve(M, tee = True)

M.x()

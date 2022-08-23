import openmdao.api as om
from mpi4py import MPI

comm = MPI.COMM_WORLD

prob = om.Problem()
model = prob.model

model.set_input_defaults('x', 1.)

parallel = model.add_subsystem('parallel', om.ParallelGroup(), 
                               promotes_inputs=[('c1.x', 'x'), ('c2.x', 'x'), 
                                                ('c3.x', 'x'), ('c4.x', 'x')])
parallel.add_subsystem('c1', om.ExecComp(['y=-2.0*x']))
parallel.add_subsystem('c2', om.ExecComp(['y=5.0*x']))
parallel.add_subsystem('c3', om.ExecComp(['y=-3.0*x']))
parallel.add_subsystem('c4', om.ExecComp(['y=4.0*x']))

model.add_subsystem('c5', om.ExecComp(['y=3.0*x1 + 7.0*x2 - 2.0*x3 + x4']))

model.connect("parallel.c1.y", "c5.x1")
model.connect("parallel.c2.y", "c5.x2")
model.connect("parallel.c3.y", "c5.x3")
model.connect("parallel.c4.y", "c5.x4")

prob.setup(check=False, mode='fwd')
prob.run_model()

om.n2(prob)

print(prob['c5.y'])
# ----------------------------------------------------------------------------
#                   Classes for setting up sellar MDA
# ----------------------------------------------------------------------------

import openmdao.api as om
import numpy as np
import pickle

# Class for Discipline 1
class SellarDis1(om.ExplicitComponent):
    """
        Component containing Discipline 1 -- no derivatives version.
    """
    def setup(self):
        # Design Variable
        self.add_input('x', val=np.zeros(3))

        # Coupling parameter
        self.add_input('y2', val=1.0)

        # Coupling output
        self.add_output('y1', val=1.0)

        # Finite difference all partials.
        self.declare_partials('*', '*', method='fd')

    def compute(self, inputs, outputs):
        """
            Evaluates the equation
            y1 = x1**2 + x2 + x3 - 0.2*y2
        """
        x1 = inputs['x'][0]
        x2 = inputs['x'][1]
        x3 = inputs['x'][2]
        y2 = inputs['y2']

        outputs['y1'] = x1**2 + x2 + x3 - 0.2*y2

# Class for Displine 2
class SellarDis2(om.ExplicitComponent):
    """
        Component containing Discipline 2 -- no derivatives version.
    """
    def setup(self):
        # Design Variable
        self.add_input('x', val=np.zeros(3))

        # Coupling parameter
        self.add_input('y1', val=1.0)

        # Coupling output
        self.add_output('y2', val=1.0)

        # Finite difference all partials.
        self.declare_partials('*', '*', method='fd')

    def compute(self, inputs, outputs):
        """
            Evaluates the equation
            y2 = sqrt(y1) + x1 + x3
        """
        x1 = inputs['x'][0]
        x2 = inputs['x'][1]
        x3 = inputs['x'][2]
        y1 = inputs['y1']

        # Note: this may cause some issues. However, y1 is constrained to be
        # above 3.16, so lets just let it converge, and the optimizer will
        # throw it out
        if y1.real < 0.0:
            y1 *= -1

        # outputs['y2'] = y1**0.5 + x1 + x2
        outputs['y2'] = y1**0.5 + x1 + x3

# Class for MDA
class SellarMDA(om.Group):
    """
        Group containing the Sellar MDA.
    """
    def readInput(self):
        # Reading input file for the analysis
        filehandler = open("input.pickle", 'rb') 
        input = pickle.load(filehandler)
        filehandler.close()

        self.x = input["x"]

    def setup(self):
        # Adding MDA components
        cycle = self.add_subsystem('cycle', om.Group(), promotes=['*'])
        cycle.add_subsystem('d1', SellarDis1(), promotes_inputs=['x', 'y2'],promotes_outputs=['y1'])
        cycle.add_subsystem('d2', SellarDis2(), promotes_inputs=['x', 'y1'],promotes_outputs=['y2'])                         

        # Nonlinear Block Gauss Seidel is a gradient free solver
        # iprint = 0 stops printing mda sovler details
        cycle.nonlinear_solver = om.NonlinearBlockGS(iprint=0, maxiter=50, use_aitken=True)

        # Adding output components
        self.add_subsystem('obj_comp', om.ExecComp('obj = x[1]**2 + x[2] + y1 + exp(-y2)', x=np.zeros(3)), promotes=['x', 'y1', 'y2', 'obj'])
        self.add_subsystem('con_comp1', om.ExecComp('con1 = 1 - y1/8.0'), promotes=['con1', 'y1'])
        self.add_subsystem('con_comp2', om.ExecComp('con2 = y2/10.0 - 1'), promotes=['con2', 'y2'])

prob = om.Problem()
prob.model = SellarMDA()
prob.model.readInput()
prob.setup()
prob['x'] = prob.model.x
prob.run_model()

output = {}
output["y"] = np.array([prob['obj'][0], prob['con1'][0], prob['con2'][0]])
output["d1_counts"] = prob.model.cycle.d1.iter_count
output["d2_counts"] = prob.model.cycle.d2.iter_count

filehandler = open("output.pickle", "xb")
pickle.dump(output, filehandler)
filehandler.close()

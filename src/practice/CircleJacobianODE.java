package practice;

import org.apache.commons.math3.ode.ExpandableStatefulODE;
import org.apache.commons.math3.ode.JacobianMatrices;
import org.apache.commons.math3.ode.MainStateJacobianProvider;
import org.apache.commons.math3.ode.ParameterJacobianProvider;
import org.apache.commons.math3.ode.nonstiff.DormandPrince853Integrator;

import java.util.Arrays;
import java.util.Collection;

public class CircleJacobianODE extends CircleODE implements MainStateJacobianProvider, ParameterJacobianProvider {
	public static final String CENTER_X = "cx";
	public static final String CENTER_Y = "cy";
	public static final String OMEGA = "omega";

	private final double[][] savedDfDp;
	private final double[] c;
	private final double omega;

	public CircleJacobianODE(double[] c, double omega) {
		super(c, omega);
		this.c = c;
		this.omega = omega;
		this.savedDfDp = new double[2][3];
	}

	public static void main(String[] args) {
		CircleJacobianODE circle = new CircleJacobianODE(new double[] {1.0, 1.0}, 0.1);

		// here, we could select only a subset of the parameters, or use another order
		JacobianMatrices jm = new JacobianMatrices(circle, CircleJacobianODE.CENTER_X, CircleJacobianODE.CENTER_Y, CircleJacobianODE.OMEGA);
		jm.addParameterJacobianProvider(circle);

		ExpandableStatefulODE efode = new ExpandableStatefulODE(circle);
		efode.setTime(0);
		double[] y = {1.0, 0.0};
		efode.setPrimaryState(y);
		System.out.println("Initial state:\n" + Arrays.toString(y));

		// create the variational equations and append them to the main equations, as secondary equations
		jm.registerVariationalEquations(efode);

		// integrate the compound state, with both main and additional equations
		DormandPrince853Integrator integrator = new DormandPrince853Integrator(1.0e-6, 1.0e3, 1.0e-10, 1.0e-12);
		integrator.setMaxEvaluations(5000);
		integrator.integrate(efode, 20.0);

		// retrieve the Jacobian of the final state with respect to initial state
		double[][] dYdY0 = new double[2][2];
		jm.getCurrentMainSetJacobian(dYdY0);

		// retrieve the Jacobian matrices of the final state with respect to the various parameters
		double[] dYdCx = new double[2];
		double[] dYdCy = new double[2];
		double[] dYdOm = new double[2];
		jm.getCurrentParameterJacobian(CircleJacobianODE.CENTER_X, dYdCx);
		jm.getCurrentParameterJacobian(CircleJacobianODE.CENTER_Y, dYdCy);
		jm.getCurrentParameterJacobian(CircleJacobianODE.OMEGA, dYdOm);

		System.out.println("dY/dCx:\n" + Arrays.toString(dYdCx));
		System.out.println("dY/dCy:\n" + Arrays.toString(dYdCy));
		System.out.println("dY/dOMEGA:\n" + Arrays.toString(dYdOm));
		System.out.println("Jacobian matrix dY/dY0:");
		for (double[] row : dYdY0) System.out.println(Arrays.toString(row));

	}

	public Collection<String> getParametersNames() {
		return Arrays.asList(CENTER_X, CENTER_Y, OMEGA);
	}

	public boolean isSupported(String name) {
		return CENTER_X.equals(name) || CENTER_Y.equals(name) || OMEGA.equals(name);
	}

	public void computeMainStateJacobian(double t, double[] y, double[] yDot, double[][] dFdY) {

		// compute the Jacobian of the main state
		dFdY[0][0] = 0;
		dFdY[0][1] = -omega;
		dFdY[1][0] = omega;
		dFdY[1][1] = 0;

		// precompute the derivatives with respect to the parameters,
		// they will be provided back when computeParameterJacobian are called later on
		savedDfDp[0][0] = 0;
		savedDfDp[0][1] = omega;
		savedDfDp[0][2] = c[1] - y[1];
		savedDfDp[1][0] = -omega;
		savedDfDp[1][1] = 0;
		savedDfDp[1][2] = y[0] - c[0];

	}

	public void computeParameterJacobian(double t, double[] y, double[] yDot,
	                                     String paramName, double[] dFdP) {
		// we simply return the derivatives precomputed earlier
		if (CENTER_X.equals(paramName)) {
			dFdP[0] = savedDfDp[0][0];
			dFdP[1] = savedDfDp[1][0];
		} else if (CENTER_Y.equals(paramName)) {
			dFdP[0] = savedDfDp[0][1];
			dFdP[1] = savedDfDp[1][1];
		} else {
			dFdP[0] = savedDfDp[0][2];
			dFdP[1] = savedDfDp[1][2];
		}
	}

}

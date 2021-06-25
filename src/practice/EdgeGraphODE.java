package practice;

import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.FirstOrderIntegrator;
import org.apache.commons.math3.ode.nonstiff.DormandPrince853Integrator;

import java.util.Arrays;

public record EdgeGraphODE(double beta, double gamma) implements FirstOrderDifferentialEquations {
	public static final int S0 = 0;
	public static final int S1 = 1;
	public static final int I0 = 2;
	public static final int I1 = 3;
	public static final int S0I1 = 4;
	public static final int I0S1 = 5;

	public static void main(String[] args) {
		FirstOrderIntegrator integrator = new DormandPrince853Integrator(1.0e-8, 100.0,
				1.0e-10, 1.0e-10);
		FirstOrderDifferentialEquations ode = new EdgeGraphODE(0.8, 0.1);

		double[] y = new double[] {0, 1, 1, 0, 0, 1}; // Initial state
		System.out.println("Initial state:\n" + Arrays.toString(y) + "\n");

		integrator.integrate(ode, 0, y, 2, y);
		System.out.println("Final state at t = 16:\n" + Arrays.toString(y));
	}


	@Override
	public int getDimension() {
		return 6;
	}

	@Override
	public void computeDerivatives(double t, double[] y, double[] yDot) throws MaxCountExceededException, DimensionMismatchException {
		yDot[S0] = -beta * y[S0I1];
		yDot[S1] = -beta * y[I0S1];
		yDot[I0] = (gamma * y[S0I1]) - (gamma * y[I0]);
		yDot[I1] = (gamma * y[I0S1]) - (gamma * y[I1]); // I1
		yDot[S0I1] = (-beta * y[S0I1]) + (gamma * y[I0S1]) - (gamma * y[I1]);
		yDot[I0S1] = (gamma * y[S0I1]) - (gamma * y[I0]) - (beta * y[I0S1]);
	}

}

package practice;

import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.FirstOrderIntegrator;
import org.apache.commons.math3.ode.nonstiff.DormandPrince853Integrator;

import java.util.Arrays;

public record TriangleGraphODE(double beta, double gamma) implements FirstOrderDifferentialEquations {
	// Singles
	public static final int S0 = 0;
	public static final int S1 = 1;
	public static final int S2 = 2;
	public static final int I0 = 3;
	public static final int I1 = 4;
	public static final int I2 = 5;
	// Doubles
	public static final int S0_I1 = 6;
	public static final int I0_S1 = 7;
	public static final int S0_I2 = 8;
	public static final int I0_S2 = 9;
	public static final int S1_I2 = 10;
	public static final int I1_S2 = 11;
	// Triples
	public static final int S0_S1_I2 = 12;
	public static final int S0_I1_S2 = 13;
	public static final int S0_I1_I2 = 14;
	public static final int I0_S1_S2 = 15;
	public static final int I0_I1_S2 = 16;
	public static final int I0_S1_I2 = 17;

	public static void main(String[] args) {
		FirstOrderIntegrator integrator = new DormandPrince853Integrator(1.0e-8, 100.0,
				1.0e-10, 1.0e-10);
		FirstOrderDifferentialEquations ode = new TriangleGraphODE(0.8, 0.1);

		double[] y = new double[18]; // Initial state
		for (int i = 0; i < y.length; i++) y[i] = Math.random();
		System.out.println("Initial state:\n" + Arrays.toString(y) + "\n");

		integrator.integrate(ode, 0, y, 2, y);
		System.out.println("Final state at t = 16:\n" + Arrays.toString(y));
	}

	@Override
	public int getDimension() {
		return 18;
	}

	@Override
	public void computeDerivatives(double v, double[] y, double[] yDot) throws MaxCountExceededException, DimensionMismatchException {
		// Singles
		yDot[S0] = (-beta * y[S0_I1]) - (beta * y[S0_I2]);
		yDot[S1] = (-beta * y[I0_S1]) - (beta * y[S1_I2]);
		yDot[S2] = (-beta * y[I1_S2]) - (beta * y[I0_S2]);

		yDot[I0] = (beta * y[S0_I1]) + (beta * y[S0_I2]) - (gamma * y[I0]);
		yDot[I1] = (beta * y[I0_S1]) + (beta * y[S1_I2]) - (gamma * y[I1]);
		yDot[I2] = (beta * y[I1_S2]) + (beta * y[I0_S2]) - (gamma * y[I2]);

		// Doubles
		yDot[S0_I1] = (-(beta + gamma) * y[S0_I1]) + (beta * y[S0_S1_I2]) - (beta * y[S0_I1_I2]);
		yDot[S0_I2] = (-(beta + gamma) * y[S0_I2]) - (beta * y[S0_I1_I2]) + (beta * y[S0_I1_S2]);
		yDot[S1_I2] = (-(beta + gamma) * y[S1_I2]) - (beta * y[I0_S1_I2]) + (beta * y[I0_S1_S2]);

		yDot[I0_S1] = (-(beta + gamma) * y[I0_S1]) - (beta * y[I0_S1_I2]) + (beta * y[S0_S1_I2]);
		yDot[I0_S2] = (-(beta + gamma) * y[I0_S2]) - (beta * y[I0_I1_S2]) + (beta * y[S0_I1_S2]);
		yDot[I1_S2] = (-(beta + gamma) * y[I1_S2]) + (beta * y[I0_S1_S2]) - (beta * y[I0_I1_S2]);

		// Triples
		yDot[S0_S1_I2] = (-((2 * beta) + gamma) * y[S0_S1_I2]);
		yDot[S0_I1_S2] = (-((2 * beta) + gamma) * y[S0_I1_S2]);
		yDot[S0_I1_I2] = (-((2 * beta) + gamma) * y[S0_I1_I2]) + (beta * y[S0_S1_I2]) + (beta * y[S0_I1_S2]);

		yDot[I0_S1_S2] = (-((2 * beta) + gamma) * y[I0_S1_S2]);
		yDot[I0_S1_I2] = (-((2 * beta) + gamma) * y[I0_S1_I2]) + (beta * y[I0_S1_S2]) + (beta * y[S0_S1_I2]);
		yDot[I0_I1_S2] = (-((2 * beta) + gamma) * y[I0_I1_S2]) + (beta * y[I0_S1_S2]) + (beta * y[S0_I1_S2]);
	}
}

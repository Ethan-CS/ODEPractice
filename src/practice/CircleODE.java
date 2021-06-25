package practice;

import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.FirstOrderIntegrator;
import org.apache.commons.math3.ode.nonstiff.DormandPrince853Integrator;
import org.apache.commons.math3.ode.sampling.StepHandler;
import org.apache.commons.math3.ode.sampling.StepInterpolator;

import java.util.Arrays;

/**
 * This example is taken from the <a href="https://commons.apache.org/proper/commons-math/userguide/ode.html">user guide
 * for the Apache Math ODE package.</a>
 * <p>
 * The following example shows how to implement the simple two-dimensional problem using double primitives:
 * <p>
 * y'0(t) = ω × (c1 - y1(t))
 * y'1(t) = ω × (y0(t) - c0)
 * with some initial state y(t0) = (y0(t0), y1(t0)). The exact solution of this problem is that y(t) moves along a
 * circle centered at c = (c0, c1) with constant angular rate ω.
 */
public class CircleODE implements FirstOrderDifferentialEquations {

	private final double[] c;
	private final double omega;

	public CircleODE(double[] c, double omega) {
		this.c = c;
		this.omega = omega;

	}

	public static void main(String[] args) {
		FirstOrderIntegrator integrator = new DormandPrince853Integrator(1.0e-8, 100.0,
				1.0e-10, 1.0e-10);
		FirstOrderDifferentialEquations ode = new CircleODE(new double[] {1.0, 1.0}, 0.1);

		StepHandler stepHandler = new StepHandler() {
			public void init(double t0, double[] y0, double t) {
			}

			public void handleStep(StepInterpolator interpolator, boolean isLast) {
				double t = interpolator.getCurrentTime();
				double[] y = interpolator.getInterpolatedState();
				System.out.println(t + " " + y[0] + " " + y[1]);
			}
		};
		integrator.addStepHandler(stepHandler);

		double[] y = new double[] {0.0, 1.0}; // The initial state
		System.out.println("Initial state: " + Arrays.toString(y));

		integrator.integrate(ode, 0.0, y, 16.0, y); // Now y contains the final state at time t=16
		System.out.println("Final state at t = 16:" + Arrays.toString(y));
	}

	@Override
	public int getDimension() {
		return 2;
	}

	@Override
	public void computeDerivatives(double t, double[] y, double[] yDot) throws MaxCountExceededException, DimensionMismatchException {
		yDot[0] = omega * (c[1] - y[1]);
		yDot[1] = omega * (y[0] - c[0]);
	}

	@Override
	public boolean equals(Object o) {
		if (this == o) return true;
		if (o == null || getClass() != o.getClass()) return false;

		CircleODE circleODE = (CircleODE) o;

		if (Double.compare(circleODE.omega, omega) != 0) return false;
		return Arrays.equals(c, circleODE.c);
	}

	@Override
	public int hashCode() {
		int result;
		long temp;
		result = Arrays.hashCode(c);
		temp = Double.doubleToLongBits(omega);
		result = 31 * result + (int) (temp ^ (temp >>> 32));
		return result;
	}

	@Override
	public String toString() {
		return "CircleODE{" +
		       "c=" + Arrays.toString(c) +
		       ", omega=" + omega +
		       '}';
	}
}

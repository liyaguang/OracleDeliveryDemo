package edu.usc.infolab.main;

import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.charset.CharsetEncoder;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;

import ilog.concert.*;
import ilog.cplex.*;

public class SimpleLPTest {

	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub
		try {
			// LPTest1();
			deliveryMIPTest();
		} catch (IloException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	private static void LPTest2() {
		try {
			IloCplex cplex = new IloCplex();
			double[] lb = { 0.0, 0.0, 0.0 };
			double[] ub = { 40.0, Double.MAX_VALUE, Double.MAX_VALUE };
			IloNumVar[] x = cplex.numVarArray(3, lb, ub);
			double[] objvals = { 1.0, 2.0, 3.0 };
			cplex.addMaximize(cplex.scalProd(x, objvals));
			cplex.addLe(cplex.sum(cplex.prod(-1.0, x[0]),
					cplex.prod(1.0, x[1]), cplex.prod(1.0, x[2])), 20.0);
			cplex.addLe(cplex.sum(cplex.prod(1.0, x[0]),
					cplex.prod(-3.0, x[1]), cplex.prod(1.0, x[2])), 30.0);
			if (cplex.solve()) {
				cplex.output()
						.println("Solution status = " + cplex.getStatus());
				cplex.output().println(
						"Solution value = " + cplex.getObjValue());
				double[] val = cplex.getValues(x);
				int ncols = cplex.getNcols();
				for (int j = 0; j < ncols; ++j)
					cplex.output().println(
							"Column: " + j + " Value = " + val[j]);
			}
			cplex.end();
		} catch (IloException e) {
			System.err.println("Concert exception '" + e + "' caught");
		}
	}

	/**
	 * Linear Programming test
	 * 
	 * <pre>
	 * Maximize 
	 * x0 + 2x1 
	 * Subject to 
	 * x0 - x1 + 5 <= 0 
	 * x0 + x1 - 10 <= 10 
	 * Bounds 
	 * x0 >= 0 
	 * x1 >= 0 
	 * End
	 * </pre>
	 * 
	 * @throws IloException
	 */
	public static void LPTest1() throws IloException {
		IloCplex cplex = new IloCplex();
		// ArrayList<IloNumExpr> exps = new ArrayList<IloNumExpr>();
		IloNumVar[] x = cplex.numVarArray(2, 0, 1e10);
		IloLinearNumExpr exp = cplex.linearNumExpr();
		exp.addTerm(1, x[0]);
		exp.addTerm(1, x[1]);
		cplex.addLe(exp, 10, "c1");
		exp = cplex.linearNumExpr();
		exp.addTerm(1, x[0]);
		exp.addTerm(-2, x[1]);
		cplex.addGe(exp, -10, "c2");
		// cplex.addLe(cplex.sum(cplex.prod(arg0, arg1)

		// objective
		IloLinearNumExpr target = cplex.linearNumExpr();
		target.addTerm(2, x[1]);
		target.addTerm(1, x[0]);
		IloObjective obj = cplex.maximize(target);
		cplex.add(obj);
		if (cplex.solve()) {
			cplex.output().println("Solution status = " + cplex.getStatus());
			cplex.output().println("Solution value = " + cplex.getObjValue());
			double[] val = cplex.getValues(x);
			int ncols = cplex.getNcols();
			for (int j = 0; j < ncols; ++j)
				cplex.output().println("Column: " + j + " Value = " + val[j]);
		}
		cplex.end();
	}

	public static void deliveryMIPTest() throws IloException, IOException {
		double[][] T = DataRetrievalTest.getDistArray();
		int N = 6;
		int TN = 30; // Time period number
		double f = 15; // time period length
		double P = 50; // penalty
		double B = 1e10;
		double[] S = new double[N]; // service time
		double[] w = new double[N]; // weight associated with customer
		double[] U = new double[N]; // upper bound of delivery time
		double[] L = new double[N]; // lower bound of delivery time
		double U2 = 60 * f;
		double maxU = 40;
		for (int i = 0; i < N; ++i) {
			S[i] = 5;
			w[i] = 0.3;
			L[i] = 0;
			U[i] = maxU;
		}
		// variables
		IloCplex cplex = new IloCplex();
		IloIntVar[] X = cplex.intVarArray(N * N * TN, 0, 1);
		IloIntVar[] Yl = cplex.intVarArray(N, 0, 1);
		IloIntVar[] Yu = cplex.intVarArray(N, 0, 1);
		IloNumVar[] T2 = cplex.numVarArray(N, 0, U2); // leave time
		// objective
		IloLinearNumExpr target = cplex.linearNumExpr();
		// add \sum_{i,j,t}{X_{ij} * T_{ijt}}
		for (int i = 0; i < N; ++i) {
			for (int j = 0; j < N; ++j) {
				for (int t = 0; t < TN; ++t) {
					target.addTerm(T[i][j], X[(i * N + j) * N + t]);
				}
			}
		}
		// \sum_i^n{Yui + Yli} * P
		for (int i = 0; i < N; ++i) {
			target.addTerm(P, Yl[i]);
			target.addTerm(P, Yu[i]);
		}
		// \sum_i^n{w * T_i}
		for (int i = 0; i < N; ++i) {
			target.addTerm(w[i], T2[i]);
		}
		cplex.addMinimize(target);

		// constraints
		// (2) T_0 = 0
		cplex.addEq(cplex.prod(1.0, T2[0]), 0, "(2)");
		// (3) \sum_{t,j} X_{ijt} = 1, (4) sum_{t,j} X_{ijt} = sum_{t,j} X_{jit}
		for (int i = 0; i < N; ++i) {
			IloLinearNumExpr exp = cplex.linearNumExpr();
			IloLinearNumExpr exp2 = cplex.linearNumExpr();
			for (int j = 0; j < N; ++j) {
				if (j == i)
					continue;
				for (int t = 0; t < TN; ++t) {
					exp.addTerm(1, X[(i * N + j) * TN + t]);
					exp2.addTerm(1, X[(j * N + i) * TN + t]);
				}
			}
			cplex.addEq(exp, 1, String.format("(3):(i=%d)", i));
			cplex.addEq(exp2, 1, String.format("(4):(i=%d)", i));
		}
		// (5) (6) (7)
		for (int t = 0; t < TN; ++t) {
			for (int i = 0; i < N; ++i) {
				for (int j = 1; j < N; ++j) {
					// (5) T_j - T_i - B X_{ijt} >= T_{ijt} + S_j - B
					IloLinearNumExpr exp = cplex.linearNumExpr();
					exp.addTerm(1, T2[j]);
					exp.addTerm(-1, T2[i]);
					exp.addTerm(-B, X[(i * N + j) * TN + t]);
					cplex.addGe(exp, T[i][j] + S[j] - B,
							String.format("(5):(%d, %d, %d)", i, j, t));
					// (6) T_i + B * X{ijt} <= f * (t + 1) + B
					exp = cplex.linearNumExpr();
					exp.addTerm(1, T2[i]);
					exp.addTerm(B, X[(i * N + j) * TN + t]);
					cplex.addLe(exp, f * (t + 1) + B,
							String.format("(6):(%d, %d, %d)", i, j, t));

					// (7) T_i >= f * t * X{ijt}
					exp = cplex.linearNumExpr();
					exp.addTerm(1, T2[i]);
					exp.addTerm(-f * t, X[(i * N + j) * TN + t]);
					cplex.addGe(exp, 0,
							String.format("(7):(%d, %d, %d)", i, j, t));
				}
			}

		}
		for (int i = 0; i < N; ++i) {
			// (8) T_i - B * Yu_i <= U_i, Upper bound constraint
			IloLinearNumExpr exp = cplex.linearNumExpr();
			exp.addTerm(1, T2[i]);
			exp.addTerm(-B, Yu[i]);
			cplex.addLe(exp, U[i], String.format("(8):(%d)", i));

			// (9) T_i + B * Yl_i >= L_i, Lower bound constraint
			exp = cplex.linearNumExpr();
			exp.addTerm(1, T2[i]);
			exp.addTerm(B, Yl[i]);
			cplex.addGe(exp, L[i], String.format("(9):(%d)", i));
		}
		if (cplex.solve()) {

			int[] deliverySeq = new int[N];
			double[] departureTimes = new double[N];
			int[][] violation = new int[2][N];
			double objVal = cplex.getObjValue();
			// Show and interpret the result
			cplex.output().println("Solution status = " + cplex.getStatus());
			// X[i,j,t]
			double[] val = cplex.getValues(X);
			int ncols = X.length;
			for (int col = 0; col < ncols; ++col)
				if (val[col] > 0) {
					int s = col / (N * TN);
					int e = (col - s * N * TN) / TN;
					int t = col % TN;
					deliverySeq[s] = e;
					// cplex.output().println(
					// String.format("X[%d,%d,%d] = %.3f", s, e, t,
					// val[col]));
				}
			// T, departure times
			val = cplex.getValues(T2);
			for (int col = 0; col < T2.length; ++col) {
				departureTimes[col] = val[col];
				// cplex.output().println(
				// String.format("T[%d] = %.3f", col, val[col]));
			}

			val = cplex.getValues(Yl);
			for (int col = 0; col < Yl.length; ++col) {
				violation[0][col] = (int) val[col];
				// cplex.output().println(
				// String.format("Yl[%d] = %.3f", col, val[col]));
			}
			// bound violation
			val = cplex.getValues(Yu);
			for (int col = 0; col < Yu.length; ++col) {
				violation[1][col] = (int) val[col];
				// cplex.output().println(
				// String.format("Yu[%d] = %.3f", col, val[col]));
			}

			interpretResult(deliverySeq, departureTimes, violation, objVal);
		}
		cplex.end();
		// interpret result
	}

	private static void interpretResult(int[] deliverySeq,
			double[] departureTimes, int[][] violation, double objVal) {
		// TODO Auto-generated method stub
		// int N = deliverySeq.length;
		int cur = 0, N = deliverySeq.length;
		System.out.println("==========================================");
		System.out.println(String.format("Objective Val: %.3f", objVal));
		do {
			System.out.println(String.format("%d->%d: %.3f", cur,
					deliverySeq[cur], departureTimes[cur]));
			cur = deliverySeq[cur];
		} while (cur != 0);
		// time window violation summary
		System.out.println("Time window violation summary:");
		for (int i = 0; i < N; ++i) {
			if (violation[0][i] == 1) {
				System.out.println(String.format(
						"Lower bound of node %d is voilated.", i));
			}
			if (violation[1][i] == 1) {
				System.out.println(String.format(
						"Upper bound of node %d is voilated.", i));
			}
		}
	}
}

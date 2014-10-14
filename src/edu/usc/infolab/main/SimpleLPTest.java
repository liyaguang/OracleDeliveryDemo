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
//			LPTest1();
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
		int N = T.length;
		int TN = 30; // Time period number
		double f = 15; // time peroid length
		double P = 50; // penalty
		double B = 1e10;
		double[] S = new double[N]; // service time
		double[] w = new double[N]; // weight associated with customer
		double[] U = new double[N]; // upper bound of delivery time
		double[] L = new double[N]; // lower bound of delivery time
		double U2 = 60 * f;
		for (int i = 0; i < N; ++i) {
			S[i] = 5;
			w[i] = 0.3;
			L[i] = 0;
			U[i] = 100;
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
					target.addTerm(T[i][j], X[(i * N + j) * TN + t]);
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
		// (2) start time = 0
		cplex.addEq(cplex.prod(1.0, T2[0]), 0, "Start Time");
		// (3) \sum_t \sum_j X_ij = 1, (4) \sum_t \sum_j X_ji = 1
		for (int i = 0; i < N; ++i) {
			IloLinearNumExpr exp = cplex.linearNumExpr();
			IloLinearNumExpr exp2 = cplex.linearNumExpr();
			for (int j = 0; j < N; ++j) {
				for (int t = 0; t < TN; ++t) {
					exp.addTerm(1, X[(i * N + j) * TN + t]);
					exp2.addTerm(1, X[(j * N + i) * TN + t]);
				}
			}
			cplex.addEq(exp, 1);
			cplex.addEq(exp2, 1);
		}
		// (5) (6) (7) 
		for (int t = 0; t < TN; ++t) {
			for (int i = 0; i < N; ++i) {
				for (int j = 1; j < N; ++j) {
					// (5) T_j - T_i - B X_{ijt} >= T_{ijt} + S_j - B
					IloLinearNumExpr exp = cplex.linearNumExpr();
					exp.addTerm(1, T2[j]);
					exp.addTerm(-1, T2[i]);
					exp.addTerm(-B, X[(i * N + j) * N + t]);
					cplex.addGe(exp, T[i][j] + S[j] - B);
					// (6) 
					exp = cplex.linearNumExpr();
					exp.addTerm(1, T2[i]);
					exp.addTerm(B, X[(i * N + j) * N + t]);
					cplex.addLe(exp, f * (t + 1) + B);
					
					// (7) 
					exp = cplex.linearNumExpr();
					exp.addTerm(1, T2[i]);
					exp.addTerm(-f * t, X[i * N + j]);
					cplex.addGe(exp, 0);
				}
			}
			
		}
		// (8) (9)
		for (int i = 0; i < N; ++i) {
			IloLinearNumExpr exp = cplex.linearNumExpr();
			exp.addTerm(1, T2[i]);
			exp.addTerm(-B, Yu[i]);
			cplex.addLe(exp, U[i]);
			
			exp = cplex.linearNumExpr();
			exp.addTerm(1, T2[i]);
			exp.addTerm(B, Yl[i]);
			cplex.addGe(exp, L[i]);
		}
		if (cplex.solve()) {
			cplex.output().println("Solution status = " + cplex.getStatus());
			cplex.output().println("Solution value = " + cplex.getObjValue());
			double[] val = cplex.getValues(X);
			// X
			int ncols = X.length;
			for (int j = 0; j < ncols; ++j)
				if (val[j] > 0) {
					int i1 = j / (N * N);
					int i2 = (j - i1 * N * N) / N;
					int i3 = j % N;
					cplex.output().println(String.format("X[%d,%d,%d] = %.3f", 
							i1, i2, i3, val[j]));
				}
			// T
			val = cplex.getValues(T2);
			for (int j = 0; j < T2.length; ++j)
				cplex.output().println(String.format("T%d = %.3f", j, val[j]));
			// bound
			val = cplex.getValues(Yu);
			for (int j = 0; j < Yu.length; ++j)
				cplex.output().println(String.format("Yu%d = %.3f", j, val[j]));
			
			val = cplex.getValues(Yl);
			for (int j = 0; j < Yl.length; ++j)
				cplex.output().println(String.format("Yl%d = %.3f", j, val[j]));
		}
		cplex.end();
	}
}

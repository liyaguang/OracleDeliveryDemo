package edu.usc.infolab.main;

import java.util.*;
import ilog.concert.*;
import ilog.cplex.*;

public class SimpleLPTest {

	public static void main(String[] args) {
		// TODO Auto-generated method stub
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
	 * Linear Programming test Maximize x0 + 2x1 Subject to x0 - x1 + 5 <= 0 x0
	 * + x1 - 10 <= 10 Bounds x0 >= 0 x1 >= 0 End
	 * 
	 * @throws IloException
	 */
	public static void LPTest1() throws IloException {
		IloCplex cplex = new IloCplex();
		// ArrayList<IloNumExpr> exps = new ArrayList<IloNumExpr>();
		IloNumVar[] x = cplex.numVarArray(2, 0, 1e10);
		IloLinearNumExpr exp = cplex.linearNumExpr();
		exp.addTerm(1, x[0]);
		exp.addTerm(-1, x[1]);
		cplex.addLe(exp, 5, "c1");
		cplex.addLe(exp, 5, "c1");
		// cplex.addLe(cplex.sum(cplex.prod(arg0, arg1)

		// objective
		IloLinearNumExpr target = cplex.linearNumExpr();
		target.addTerm(1, x[0]);
		target.addTerm(2, x[1]);
		IloObjective obj = cplex.maximize(target);
		cplex.add(obj);
	}

}

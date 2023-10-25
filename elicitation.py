import os
import time
import random
import copy
import numpy as np
#env "CFLAGS=-I/opt/homebrew/include -L/opt/homebrew/lib" pip install pycddlib
from gurobipy import Model, LinExpr, GRB

from constraint_utils import get_constraints, copy_polytope, get_obj
import Polytope
from diversipy import polytope
import logging
import json


def sample_from_polytope(W, num_rollouts):
    constraints = W.formatted_constraints
    A = []
    b = []
    for j in range(len(constraints)):
        a = []
        items = constraints[j].split(' ')
        sign = items[-2]
        for i in items[:-2]:
            if 'w' not in i and i != '' and i != '+':
                if sign == '>=':
                  item = float(i)*(-1) if float(i) != 0 else float(i)
                else:
                  item = float(i)
                a.append(item)
        b.append(float(items[-1]))
        A.append(a)

    X = polytope.sample(
        n_points=num_rollouts,
        lower=[0]*len(A[0]),
        upper=[1]*len(A[0]),
        A1=np.array(A),
        b1=np.array(b)
        )
    return X

def PMR_LP_gurobi(s, t, W, model, gurVars):
    obj_function_lin = ''
    for i in range(len(s)):
        obj_function_lin += str(t[i] - s[i]) + ' ' + W.vars[i]
        if i != len(s) - 1:
            obj_function_lin += ' + '
    var_names, obj_function_coefficients = get_obj(obj_function_lin)
    obj = 0
    for i, coef in enumerate(obj_function_coefficients):
        obj += coef * gurVars[i]
    model.setObjective(obj, GRB.MAXIMIZE)
    model.optimize()
    if GRB.OPTIMAL != 3:
        alpha = model.objVal
    else:
        raise Exception('Unfeasible liner programming problem computing MR')
    return alpha


def MR_LP_gurobi(x, data, W, model, gurVars):
    MRs = -float('inf')
    opponentForx = None

    for y in data:
        regret = PMR_LP_gurobi(x, y, W, model, gurVars)
        if regret > MRs:
            MRs = regret
            opponentForx = y

    return MRs, opponentForx

def mMR_LP_gurobi(data, W):
    mMR = float('inf')
    currentSolution = None
    opponent = None
    formatted_constraints = copy.deepcopy(W.formatted_constraints)
    coefficients_matrix, constant_terms, constraints_sense = get_constraints(formatted_constraints)

    model = Model()
    model.Params.LogToConsole = 0
    gurVars = {}
    var_names = ['w' + str(i) for i in range(num_vars)]
    for ix, var in enumerate(var_names):
        gurVars[ix] = model.addVar(vtype="C", name=var_names[ix])

    for i, (var_names, coeffs) in enumerate(coefficients_matrix):
        consVars = [gurVars[idx] for idx, i in enumerate(coeffs)]
        expr = LinExpr(coeffs, consVars)
        model.addLConstr(expr, sense=constraints_sense[i], rhs=constant_terms[i], name='c' + str(i))

    for x in data:
        MRx, opponentForx = MR_LP_gurobi(x, data, W, model, gurVars)
        if MRx < mMR:
            mMR = MRx
            opponent = opponentForx
            currentSolution = x
    model.terminate()
    return mMR, currentSolution, opponent, W


def compare_alternatives(currentSolution, opponent, dm_ws):
    x = np.array(currentSolution)
    y = np.array(opponent)

    u_x = dm_ws.dot(x)
    u_y = dm_ws.dot(y)
    if u_x >= u_y:
        coeffs = x - y
    else:
        coeffs = y - x
    return coeffs


def get_set_alternatives_gurobi(db, W, k):
    Q_prime = {}
    MRs = {}
    formatted_constraints = copy.deepcopy(W.formatted_constraints)
    coefficients_matrix, constant_terms, constraints_sense = get_constraints(formatted_constraints)

    model = Model()
    model.Params.LogToConsole = 0
    gurVars = {}
    var_names = ['w' + str(i) for i in range(num_vars)]
    for ix, var in enumerate(var_names):
        gurVars[ix] = model.addVar(vtype="C", name=var_names[ix])

    for i, (var_names, coeffs) in enumerate(coefficients_matrix):
        consVars = [gurVars[idx] for idx, i in enumerate(coeffs)]
        expr = LinExpr(coeffs, consVars)
        model.addLConstr(expr, sense=constraints_sense[i], rhs=constant_terms[i], name='c' + str(i))

    for x in db:
        MRx, opponentForx = MR_LP_gurobi(x, db, W, model, gurVars)
        MRs[tuple(x)] = MRx
        Q_prime[tuple(x)] = opponentForx
    model.terminate()
    MRs_sorted = sorted(MRs, key=MRs.get, reverse=False)[:k]
    Q_prime_sorted = {x: Q_prime[x] for x in MRs_sorted}
    return Q_prime_sorted


def run_simulation(data, a, b, W, dm_ws, vars, eps=0.5, constant=0, sign=">="):
    W_prime = copy_polytope(W, vars)
    coeffs = compare_alternatives(a, b, dm_ws)
    constraint = W_prime.get_formatted_constraint_from_vectors(vars, coeffs, sign, constant)
    W_prime.add_formatted_constraints(np.array([constraint]))
    maxRegret, a, b, W_prime = mMR_LP_gurobi(data, W_prime)

    count = 0
    while maxRegret > eps:
        # logging.info("maxRegret: {}, count: {}".format(maxRegret, count))
        coeffs = compare_alternatives(a, b, dm_ws)
        constraint = W_prime.get_formatted_constraint_from_vectors(vars, coeffs, sign, constant)
        W_prime.add_formatted_constraints(np.array([constraint]))
        maxRegret, a, b, W_prime = mMR_LP_gurobi(data, W_prime)

        count += 1

    return count

def run_simulation2(data, a, b, W, dm_ws, vars, eps=0.5, constant=0, sign=">="):
    W_prime = copy_polytope(W, vars)
    coeffs = compare_alternatives(a, b, dm_ws) # CSS with random wights
    constraint = W_prime.get_formatted_constraint_from_vectors(vars, coeffs, sign, constant)
    W_prime.add_formatted_constraints(np.array([constraint]))
    count_queries, _, _, _, _ = run_mcts_elicitation1(data, W_prime, vars, num_rollouts=10)
    # print(maxRegret )
    # count = 0
    # maxRegret, a, b, W_prime = mMR_LP_gurobi(data, W_prime)
    # while maxRegret > eps:
    #     print("maxRegret: %s, count: %s" % (maxRegret, count))
    #     logging.info("maxRegret: %s, count: %s" % (maxRegret, count))
    #     coeffs = compare_alternatives(a, b, dm_ws)
    #     constraint = W_prime.get_formatted_constraint_from_vectors(vars, coeffs, sign, constant)
    #     W_prime.add_formatted_constraints(np.array([constraint]))
    #     maxRegret, a, b, W_prime = run_mcts_elicitation1(data, W_prime, vars, num_rollouts=10)
    #
    #     count += 1

    return count_queries

def run_mcts_elicitation(data, W, vars, num_rollouts=5, eps=0.5, constant=0, sign=">="):
    true_max_regret = float('inf')
    true_dm_weights = sample_from_polytope(W, 1)[0]
    count_queries = 0
    while true_max_regret > eps:
        Q_prime = get_set_alternatives_gurobi(data, W, k)
        scores = {}
        dm_weights = sample_from_polytope(W, num_rollouts)

        for x, y in Q_prime.items(): # k =2
            counts = []
            for i in range(num_rollouts):

                # Run rollout
                start = time.process_time()
                count = run_simulation(data, x, y, W, dm_weights[i], vars, eps)
                counts.append(count)
                logging.info("Rollout %s score: %s, time: %s" %(i, count, time.process_time() - start))

            # Compute average score for all rollouts
            scores[tuple([tuple(x), tuple(y)])] = np.mean(counts)
            logging.info("Average score for the query: %s" % np.mean(counts))

        # Get the pair with lowest average score
        a, b = min(scores, key=scores.get)
        coeffs = compare_alternatives(a, b, true_dm_weights)
        constraint = W.get_formatted_constraint_from_vectors(vars, coeffs, sign, constant)
        W.add_formatted_constraints(np.array([constraint]))
        count_queries += 1
        true_max_regret, a, b, W = mMR_LP_gurobi(data, W)
        logging.info("---------------------")
        logging.info('True max regret %s' % true_max_regret)
        logging.info("---------------------")
    # save W constraints
    logging.info("For MCTS elicitation %s queries are needed" % count_queries)
    return count_queries, true_max_regret, a, b


def run_mcts_elicitation1(data, W, vars, num_rollouts=5, eps=0.5, constant=0, sign=">="):
    true_max_regret = float('inf')
    true_dm_weights = sample_from_polytope(W, 1)[0]
    count_queries = 0
    while true_max_regret > eps:
        Q_prime = get_set_alternatives_gurobi(data, W, k)
        scores = {}
        dm_weights = sample_from_polytope(W, num_rollouts)

        for x, y in Q_prime.items():  # k =2
            counts = []
            for i in range(num_rollouts):
                # Run rollout
                start = time.process_time()
                count = run_simulation(data, x, y, W, dm_weights[i], vars, eps)
                counts.append(count)
                # logging.info("Rollout %s score: %s, time: %s" % (i, count, time.process_time() - start))

            # Compute average score for all rollouts
            scores[tuple([tuple(x), tuple(y)])] = np.mean(counts)
            # logging.info("Average score for the query: %s" % np.mean(counts))

        # Get the pair with lowest average score
        a, b = min(scores, key=scores.get)
        coeffs = compare_alternatives(a, b, true_dm_weights)
        constraint = W.get_formatted_constraint_from_vectors(vars, coeffs, sign, constant)
        W.add_formatted_constraints(np.array([constraint]))
        count_queries += 1
        true_max_regret, a, b, W = mMR_LP_gurobi(data, W)
        # logging.info("---------------------")
        logging.info('True max regret of simulation mcts_elicitation %s' % true_max_regret)
        # logging.info("---------------------")
    # save W constraints
    logging.info("For MCTS elicitation %s queries are needed" % count_queries)
    return count_queries, true_max_regret, a, b, W


def run_mcts_elicitation2(data, W, vars, num_rollouts=5, eps=0.5, constant=0, sign=">="):
    true_max_regret = float('inf')
    true_dm_weights = sample_from_polytope(W, 1)[0]
    count_queries = 0
    while true_max_regret > eps:
        Q_prime = get_set_alternatives_gurobi(data, W, k)
        scores = {}
        dm_weights = sample_from_polytope(W, num_rollouts)

        for x, y in Q_prime.items():  # k =2
            counts = []
            for i in range(num_rollouts):
                # Run rollout

                start = time.process_time()
                count = run_simulation2(data, x, y, W, dm_weights[i], vars, eps)
                counts.append(count)
                logging.info("Rollout %s score: %s, time: %s" % (i, count, time.process_time() - start))

            # Compute average score for all rollouts
            scores[tuple([tuple(x), tuple(y)])] = np.mean(counts)
            logging.info(" : %s" % np.mean(counts))

        # Get the pair with lowest average score
        a, b = min(scores, key=scores.get)
        coeffs = compare_alternatives(a, b, true_dm_weights)
        constraint = W.get_formatted_constraint_from_vectors(vars, coeffs, sign, constant)
        W.add_formatted_constraints(np.array([constraint]))
        count_queries += 1
        true_max_regret, a, b, W = mMR_LP_gurobi(data, W)
        logging.info("---------------------")
        logging.info('True max regret %s' % true_max_regret)
        logging.info("---------------------")
    # save W constraints
    logging.info("For MCTS elicitation %s queries are needed" % count_queries)
    return count_queries, true_max_regret, a, b


def run_elicitation(data, W, vars, eps=0.5, constant=0, sign=">="):
    maxRegret = float('inf')
    dm_weights = sample_from_polytope(W, 1)[0]

    count_queries = 0
    while maxRegret > eps:

        maxRegret, a, b, W = mMR_LP_gurobi(data, W)
        coeffs = compare_alternatives(a, b, dm_weights) # Query
        constraint = W.get_formatted_constraint_from_vectors(vars, coeffs, sign, constant)
        W.add_formatted_constraints(np.array([constraint]))
        count_queries += 1

    logging.info("For standard elicitation %s queries are needed" % count_queries)
    return count_queries

if __name__ == "__main__":


    # Try this first
    maxsize = 10000
    num_vars = 10 # number of criterias: for both increase to 10 to make eliciitation more difficult
    num_items = 30 # increase the items: 100?
    k = 2
    eps = 0.5
    num_rollouts = 10
    num_runs = 10
    vars = ['w' + str(i) for i in range(num_vars)]
    data = [[random.randrange(0, 10) for i in range(num_vars)] for i in range(num_items)] # create a dataset

    if not os.path.exists("logs"):
        os.makedirs("logs")
    logging.basicConfig(level=logging.DEBUG, filename="logs/logfile%sitems_%srollouts" % (num_items, num_rollouts),
                        filemode="a+",
                        format="%(asctime)-15s %(levelname)-8s %(message)s")
    runs = {}
    if not os.path.exists("runs"):
        os.makedirs("runs")
    for i in range(num_runs):
        seed_value = random.randrange(maxsize)
        random.seed(seed_value)
        W = Polytope.Polytope(vars, frac=True)
        mcts_queries, true_max_regret, x, y = run_mcts_elicitation(data, W, vars, num_rollouts=num_rollouts)
        W = Polytope.Polytope(vars, frac=True)
        queries = run_elicitation(data, W, vars)
        runs[i] = [seed_value, queries, mcts_queries]

        with open('run_info_%sitems_%srollouts.json'%(num_items, num_rollouts), 'w') as f:
            json.dump(runs, f)



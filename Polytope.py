"""
Interface for the library pycddlib to
define user-friendly constraints.
This class is used to represent the user preferences :math:`w`.
"""

import cdd
import Utils
from fractions import Fraction
import numpy as np
import time
import sys
import copy


class Polytope:
    """
    Class used to represents a convex p-dimensional polytope :math:`W` where each vector
    :math:`w = (w_1, ..., w_p)\in W` is such that
    :math:`w_i \geq 0, w_i \leq 1` and :math:`w_1 + ... + w_n = 1`.

    This class is an interface to interact with the library pycddlib defining constraints
    with a user-friendly syntax.
    
    
    :ivar vars: List of strings: id used to represent the dimensions of the polytope.
    :ivar n_vars: Integer: number of dimension of the polytope.
    :ivar formatted_constraints: List of string: formatted constraints defining the polytope.
    :ivar extreme_points: List of p-dimensional float: extreme points of the polytope considering the
        formatted constraints.

    """

    def __init__(self, vars, frac=True, epi_var=''):
        """
        Initialize polytope as the set W of vectors w such that
        :math:`w_i \geq 0, w_i \leq 1` and :math:`w_1 + ... + w_n = 1`.

        :param vars: list of string: List of id used to represent the dimensions of the polytope.
        :param frac: boolean: If True, the coefficients of the formatted constraints must be
            fractions. This increase the precision of the polytope representation and avoid numerical errors.
        :param epi_var: string: Name of variable used to represent the epigraph value

        """

        #constraints vars[i] >= 0
        formatted_constraints = []
        for j in range(len(vars)):
            if vars[j] == epi_var: continue
            l_constraints = ''
            for i in range(len(vars)):
                if i != j:
                    l_constraints += ' 0 ' + vars[i]
                else:
                    l_constraints += ' 1 ' + vars[i]
                if i < len(vars) - 1:
                    l_constraints += ' + '
                else:
                    l_constraints += ' >= 0'
            formatted_constraints.append(l_constraints)

        # constraint sum_i vars[i] = 1 expresssed as sum_i vars[i] >= 1 and sum_i vars[i] <= 1
        sum_constraint1 = ''
        sum_constraint2 = ''
        for i in range(len(vars)):
            if vars[i] == epi_var: continue
            sum_constraint1 += ' 1 ' + vars[i]
            sum_constraint2 += ' 1 ' + vars[i]
            if i < len(vars) - 1:
                sum_constraint1 += ' + '
                sum_constraint2 += ' + '
            else:
                sum_constraint1 += ' >= 1'
                sum_constraint2 += ' <= 1'
        formatted_constraints.append(sum_constraint1)
        formatted_constraints.append(sum_constraint2)



        self.vars = vars
        self.n_vars = len(vars)
        self.formatted_constraints = formatted_constraints
        self.__matrix = Polytope.__get_matrix_from_constraints(formatted_constraints, vars, frac)
        self.__cdd_polytope = Polytope.__get_cdd_polytope_from_matrix(self.__matrix, frac)
        self.extreme_points = Polytope.__get_extreme_points_from_cdd_polytope(self.__cdd_polytope)
        self.frac = frac


    def add_formatted_constraints(self, formatted_constraints):
        """
        Method to add a constraints to the polytope.

        :param formatted_constraints: list of string. Each string represents a constraint;
            constraint format:  \\number \\space \\var1 \\space + ... + \\space \\number \\space \\varp \\space {<=, >=} \\number.
            Example: '2.0 w1 + 0 w2 <= 1'

        """
        for constraint in formatted_constraints:
            self.formatted_constraints.append(constraint)
        M = Polytope.__get_matrix_from_constraints(formatted_constraints, self.vars, self.frac)
        self.__matrix = np.vstack((self.__matrix, M))
        self.__cdd_polytope = Polytope.__get_cdd_polytope_from_matrix(self.__matrix, self.frac)
        self.extreme_points = Polytope.__get_extreme_points_from_cdd_polytope(self.__cdd_polytope)

    @staticmethod
    def get_formatted_constraint_from_vectors(vars, vals, sign, constant):
        """
        Method to generate a formatted constraint.

        :param vars: list of string: Id of variables.
        :param vals: list of float: Coefficients of variables.
        :param sign: string: '<=' or '>='.
        :param constant: float: Constant value
        :return: string: Formatted constraint.

        """
        if len(vars) != len(vals):
            sys.exit('different length vars and vals')
        c = ''
        first = True
        for i in range(len(vars)):
            if first:
                first = False
            else:
                c += ' + '
            c +='%f %s' % (vals[i], vars[i])
        c+= ' %s %f' % (sign, constant)
        return c

    @staticmethod
    def __get_extreme_points_from_constraints(formatted_constraints, vars, frac=False):

        M = Polytope.__get_matrix_from_constraints(formatted_constraints, vars, frac)
        #print(M)

        start_time = time.time()
        poly = Polytope.__get_cdd_polytope_from_matrix(M, frac)
        extreme_points = Polytope.__get_UL_bounds_matrix_constraints(poly)

        tot_time = time.time() - start_time
        #print('time to compute top points of A: %.4f; n extreme points: %i' % (tot_time, len(extreme_points)))

        if len(extreme_points) == 0:
            print('Null set of extreme points')
            sys.exit()
        return extreme_points

    @staticmethod
    def __get_matrix_from_constraints(formatted_constraints, vars, frac=False):

        # For a polyhedron described as P = {x | A x <= b}, the H-representation is the matrix [b -A].

        n_vars = len(vars)

        #creating 0 matrix: [b -A];
        #rows: n formatted_constraints
        #cols: constant term + n vars = n+1
        matrix = [[0 for j in range(1 + n_vars)] for i in range(len(formatted_constraints))]


        for i in range(len(formatted_constraints)):
            constraint = formatted_constraints[i]
            if '<=' in constraint:
                split_elmnt = '<='
                sign = 1
            elif '>=' in constraint:
                split_elmnt = '>='
                sign = -1
            else:
                continue

            #print(constraint)
            if frac:
                costant_term = Fraction(constraint.split(split_elmnt)[1].strip())
            else:
                costant_term = float(constraint.split(split_elmnt)[1].strip())
            matrix[i][0] = sign*costant_term

            constraint_elms = list(constraint.split(split_elmnt)[0].strip().split('+'))
            constraint_elms = [x.strip() for x in constraint_elms]
            for element in constraint_elms:
                var_name = element.split(' ')[1].strip()
                if frac:
                    coefficient = Fraction(element.split(' ')[0].strip())
                else:
                    coefficient = float(element.split(' ')[0].strip())

                j = Utils.index_containing_substring(vars, var_name) + 1
                matrix[i][j] = -sign*coefficient

        return matrix

    @staticmethod
    def __get_cdd_polytope_from_matrix(M, frac=False):
        if frac:
            mat = cdd.Matrix(M, number_type='fraction')
        else:
            mat = cdd.Matrix(M, number_type='float')
        mat.rep_type = cdd.RepType.INEQUALITY
        poly = cdd.Polyhedron(mat)
        return poly

    @staticmethod
    def __get_cdd_polytope_from_constraints(formatted_constraints, vars, frac=False):
        M = Polytope.__get_matrix_from_constraints(formatted_constraints, vars, frac)
        return Polytope.__get_cdd_polytope_from_matrix(M, frac)

    @staticmethod
    def __get_extreme_points_from_cdd_polytope(polytope):
        # get_generators return a matrix;
        # a row represents an extreme point:
        #   the extreme point start at index 1
        ext_points = polytope.get_generators()
        ext_points_list = set([])
        for point in ext_points:
            if point[0] == 0:
                continue
            ext_point = point[1:]
            ext_points_list.add(ext_point)
        return(ext_points_list)

    @staticmethod
    def __get_UL_bounds_matrix_constraints_and_sum_constraint(l_constraints, u_constraints, vars, frac=False):

        # For a polyhedron described as P = {x | A x <= b}, the H-representation is the matrix [b -A].
        M1 = Polytope.__get_UL_bounds_matrix_constraints(vars, l_constraints, u_constraints, frac)
        M2 = Polytope.__get_sum_constraint(vars)
        M = np.vstack((M1, M2))
        return M

    @staticmethod
    def __get_UL_bounds_matrix_constraints(l_constraints, u_constraints, frac=False):
        # For a polyhedron described as P = {x | A x <= b}, the H-representation is the matrix [b -A].

        n_vars = len(l_constraints)

        #creating 0 matrix: [b -A];
        #rows: n upperbound constraint, n lowerbound constraint = 2n
        #cols: constant term + n vars = n+1

        matrix = []

        for i in range(n_vars):
            # [b -A]
            matrix_row_l = [0] * (1 + n_vars)
            matrix_row_u = [0] * (1 + n_vars)

            # setting b value; note that we need to change the sign of the lower bound coefficent
            if frac:
                matrix_row_l[0] = -Fraction(l_constraints[i])
                matrix_row_u[0] = Fraction(u_constraints[i])
            else:
                matrix_row_l[0] = -l_constraints[i]
                matrix_row_u[0] = u_constraints[i]

            # setting scalar coefficient -a_i of current var i to 1;
            # note that we need to change the sign of the lower bound coefficent
            matrix_row_l[i + 1] = -(-1)
            matrix_row_u[i + 1] = -1

            matrix.append(matrix_row_l)
            matrix.append(matrix_row_u)


        return matrix

    @staticmethod
    def __get_sum_constraint(vars):

        n_vars = len(vars)

        #sum constraints

        #setting scalar coefficient to 1 for all variables; note -A and we need to change the sign of the lower bound coefficent
        sum_constraint_l = [-(-1)] * (1 + n_vars)
        sum_constraint_u = [-1] * (1 + n_vars)


        #setting constant term b; note that we need to change the sign of the lower bound coefficent
        sum_constraint_l[0] = -1
        sum_constraint_u[0] = 1

        matrix = []

        matrix.append(sum_constraint_l)
        matrix.append(sum_constraint_u)

        return matrix

    @staticmethod
    def __add_formatted_constraints_to_matrix(M, formatted_constraints, vars, frac=False):
        Mc = copy.deepcopy(M)
        return np.vstack((Mc, Polytope.__get_matrix_from_constraints(formatted_constraints, vars, frac)))

"""
Set of methods for common operations
"""
from fractions import Fraction

def dot(a, b):
   """
   Dot product of two vectors. The input vectors must have the same size.

    :param a: list of float.
    :param b: list of float.
    :return: float.
   """
   if len(a) != len(b):
       raise Exception('The input vectors have different length')
   return sum([a_i * b_i for a_i, b_i  in zip(a, b)])

def vector_sub(a, b):
   """
   Pointwise sum of two vectors. The input vectors must have the same size.

    :param a: list of float.
    :param b: list of float.
    :return: list of float.
   """
   if len(a) != len(b):
       raise Exception('The input vectors have different length')
   return [a_i - b_i for a_i, b_i in zip(a, b)]

def vector_sum(a, b):
   """
   Pointwise difference of two vectors. The input vectors must have the same size.

    :param a: list of float.
    :param b: list of float.
    :return: list of float.
   """
   if len(a) != len(b):
       raise Exception('The input vectors have different length')
   return [a_i + b_i for a_i, b_i in zip(a, b)]

def index_containing_substring(vars_list, var):
    """
    Return the index of the first string in vars_list that contains the substring var

    :param vars_list: list of string.
    :param var: string.
    :return: integer.
    """
    for i, s in enumerate(vars_list):
        if var in s:
              return i
    return -1

def get_cplex_constraints(obj_function, formatted_constraints):
    """
    Returns objective function and constraints formatted for Cplex. For details see the Cplex documentation.

    :param obj_function: string.
        Constraint format:  \\number \\space \\var1 \\space + ... + \\space \\number \\space \\varp \\space.
        Example: '0 w1 + 1.1 w2'
    :param formatted_constraints: list of string. Each string represents a contraint;
        constraint format:  \\number \\space \\var1 \\space + ... + \\space \\number \\space \\varp \\space {<=, ==} \\number.
        Example: '2.0 w1 + 0 w2 <= 1'
    :return: var_names: list of string representing the id of the variables.

        obj_function_coefficients: coefficient of the objective function.

        coefficients_matrix: matrix of the coefficients of the formatted constraints.

        constant_terms: vector of the constant terms of the formatted constraints.

        constraints_sense: vector of the constraints sense of the formatted constraints.
    """

    obj_function_vals = obj_function.split('+')
    obj_function_coefficients = list([])
    var_names = list([])
    for val in obj_function_vals:
        val = val.strip()
        coefficient = float(val.split(' ')[0].strip())
        obj_function_coefficients.append(coefficient)
        var_names.append(val.split(' ')[1].strip())

    constant_terms = list([])
    coefficients_matrix = list([])
    constraints_sense = list([])


    for constraint in formatted_constraints:
        if '==' in constraint:
            split_elmnt = '=='
            constraints_sense.append('E')
        elif '<=' in constraint:
            split_elmnt = '<='
            constraints_sense.append('L')
        else:
            split_elmnt = '>='
            constraints_sense.append('G')

        costant_term = float(Fraction(constraint.split(split_elmnt)[1].strip()))
        constant_terms.append(costant_term)

        constraint_vals = list(constraint.split(split_elmnt)[0].strip().split('+'))
        constraint_coefficients = list([])
        constraint_vars = list([])
        for val in constraint_vals:
            val = val.strip()
            coefficient = float(Fraction(val.split(' ')[0].strip()))
            constraint_coefficients.append(coefficient)
            var = val.split(' ')[1].strip()
            constraint_vars.append(var)

        constraint_left = list([constraint_vars, constraint_coefficients])
        coefficients_matrix.append(constraint_left)

    return var_names, obj_function_coefficients, coefficients_matrix, constant_terms, constraints_sense


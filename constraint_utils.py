from fractions import Fraction
import numpy as np
import Polytope
import copy

def copy_polytope(W_orig, vars):
    W_copy = Polytope.Polytope(vars, frac=True)
    for c in W_orig.formatted_constraints:
      W_copy.add_formatted_constraints(np.array([c]))

    return W_copy


def get_constraints(formatted_constraints):
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
    return coefficients_matrix, constant_terms, constraints_sense

def get_obj(obj_function):
    obj_function_vals = obj_function.split('+')
    obj_function_coefficients = list([])
    var_names = list([])
    for val in obj_function_vals:
        val = val.strip()
        coefficient = float(val.split(' ')[0].strip())
        obj_function_coefficients.append(coefficient)
        var_names.append(val.split(' ')[1].strip())

    return var_names, obj_function_coefficients

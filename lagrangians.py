import numpy as np


class LMLagrangian():
    def __init__(self, fun, init_guess, args, constraints) -> None:
        self.args = args
        self.init_guess = init_guess
        self.dim = int(len(init_guess)/2)
        self.fun = fun
        self.n_constr = len(constraints)
        self.constraints = constraints        

        G = np.zeros((self.n_constr, self.dim*2))
        for n, constr in enumerate(constraints):
            jac = constr['args'][0]
            G[n] = jac
        self.constr_mat = np.copy(G)

    def get_guess(self):
        guess = np.zeros(self.dim*2 + self.n_constr)
        guess[:self.dim*2] = self.init_guess
        return guess

    def min_func(self, x_in, *args):
        lam = x_in[self.dim*2:]
        x = x_in[0:self.dim*2]
        f, df = self.fun(x, *self.args)
        
        constr_vals = np.zeros(self.n_constr)
        for n, constr in enumerate(self.constraints):
            constr_vals[n] = constr['fun'](x, constr['args'])
        f += lam @ constr_vals

        deriv = np.zeros_like(x_in)
        deriv[0:self.dim*2] = df - self.constr_mat.T @ lam
        deriv[self.dim*2:] = constr_vals

        return f, deriv

class CostLagrangian():
    def __init__(self, fun, init_guess, args, constraints, cost) -> None:
        self.args = args
        self.init_guess = init_guess
        self.dim = int(len(init_guess)/2)
        self.fun = fun
        self.n_constr = len(constraints)
        self.constraints = constraints
        self.cost = cost

        G = np.zeros((self.n_constr, self.dim*2))
        for n, constr in enumerate(constraints):
            jac = constr['args'][0]
            G[n] = jac
        self.constr_mat = np.copy(G)

    def get_guess(self):
        return self.init_guess

    def min_func(self, x_in, *args):
        x = x_in[0:self.dim*2]
        f, df = self.fun(x, *self.args)
        
        constr_vals = np.zeros(self.n_constr)
        for n, constr in enumerate(self.constraints):
            constr_vals[n] = constr['fun'](x, constr['args'])
        f += self.cost * constr_vals @ constr_vals

        deriv = np.zeros_like(x_in)
        deriv[0:self.dim*2] = df + self.constr_mat.T @ constr_vals*self.cost

        return f, deriv


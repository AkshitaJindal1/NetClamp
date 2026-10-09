import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve
from deap import base, creator, tools, algorithms
import random
import matplotlib.pyplot as plt
import os

# Redefine necessary functions
def opsinmodel_ODE(t, y, vp, light):
    O1, C2, O2 = y
    C1 = 1 - (O1 + O2 + C2)

    if len(vp) >= 11:
        Gd1, Gd2, Gr, e12, e21, epsilon1, epsilon2, I, lambda_, sigma, w = vp[:11]
    else:
        raise ValueError(f"Insufficient number of parameters. Expected 11, got {len(vp)}.")
        
    I = light(t)

    Phi = sigma * I * lambda_ / w
    k1 = epsilon1 * Phi
    k2 = epsilon2 * Phi

    dO1_dt = (k1 * C1) + (e21 * O2) - O1 * (e12 + Gd1)
    dC2_dt = (Gd2 * O2) - C2 * (Gr + k2)
    dO2_dt = (k2 * C2) + (e12 * O1) - O2 * (e21 + Gd2)

    return [dO1_dt, dC2_dt, dO2_dt]

def setDefaultParameters_opsinmodel_ODE():
    modelpars = {
        'Gd1': 0.1043,
        'Gd2': 0.3,
        'Gr': 5.494e-4,
        'e12': 0.079,
        'e21': 0.25,
        'epsilon1': 0.65,
        'epsilon2': 0.35,
        'I': 0,
        'lambda': 470,
        'sigma': 6e-4,
        'w': 0.77
    }

    varpars = list(modelpars.values())
    sysPars = {'tstart': 0, 'tfinal': 600}
    fun = lambda t, y, vp, light: opsinmodel_ODE(t, y, vp, light)
    return modelpars, varpars, sysPars, fun

def light_fun(t):
    return (t >= 100) * (t <= 500) * 0.77

def computeResidual(pars, fun, Data):
    try:
        varpars = pars + [470, 6e-4, 0.77]  # Append fixed parameters
        y0 = [0.55, 0.225, 0.225]
        y0 = fsolve(lambda y: fun(0, y, varpars, light_fun), y0, xtol=1e-6)

        # Ensure fsolve converged
        if not np.all(np.isfinite(y0)):
            return float('inf'),

        # ODE solver and fitness calculation
        t_start, t_end, dt = 0, 600, 10
        t_eval = np.arange(t_start, t_end + dt, dt)
        sol = solve_ivp(fun, [t_start, t_end], y0, t_eval=t_eval, args=(varpars, light_fun), rtol=1e-8, atol=1e-10)

        if sol.success:
            O1, O2 = sol.y[0], sol.y[2]
            G, I_leak = pars[6], pars[7]
            gamma, E, V = 0.05, 0, -80

            I_opsin = G * (O1 + gamma * O2) * (V - E) - I_leak
            I_interp = np.interp(Data['t'].values, sol.t, I_opsin)
            
            res = np.sum((I_interp - Data['I'].values) ** 2)

            # Print fitness (residual) value for debugging
            print(f"Fitness (Residual) for params {pars}: {res}")
            
            return res,
        else:
            return float('inf'),
    except Exception as e:
        print(f"Error during evaluation: {e}")
        return float('inf'),

# Run GA
def run_GA():
    # Load Data
    Data = pd.read_csv(os.path.join('Data_2d.csv'))
    t_exp, I_exp = Data['t'].values, Data['I'].values

    modelpars, varpars, sysPars, fun = setDefaultParameters_opsinmodel_ODE()

    creator.create("FitnessMin", base.Fitness, weights=(-1.0,))  # Minimize fitness
    creator.create("Individual", list, fitness=creator.FitnessMin)

    toolbox = base.Toolbox()
    param_bounds = [
        (0.01, 0.1), (0.01, 0.1), (1e-5, 1e-3), (0.01, 0.1),
        (0.01, 0.5), (0.1, 5), (0.1, 0.5), (1, 5), (0.0, 40)
    ]
    
    # Register individual creation functions
    for i, (low, high) in enumerate(param_bounds):
        toolbox.register(f"attr_param{i}", random.uniform, low, high)

    toolbox.register("individual", tools.initCycle, creator.Individual,
                     [toolbox.__getattribute__(f"attr_param{i}") for i in range(len(param_bounds))])
    toolbox.register("population", tools.initRepeat, list, toolbox.individual)
    toolbox.register("mate", tools.cxBlend, alpha=0.3)
    toolbox.register("mutate", tools.mutPolynomialBounded, low=[low for low, _ in param_bounds], up=[high for _, high in param_bounds], eta=1.0, indpb=0.05)
    toolbox.register("select", tools.selTournament, tournsize=3)
    toolbox.register("evaluate", lambda ind: computeResidual(ind, fun, Data))

    # Initialize population
    population = toolbox.population(n=100)
    NGEN, CXPB, MUTPB = 100, 0.7, 0.2
    best_per_gen = []

    # Evolutionary loop
    for gen in range(NGEN):
        offspring = algorithms.varAnd(population, toolbox, CXPB, MUTPB)
        fits = toolbox.map(toolbox.evaluate, offspring)

        for fit, ind in zip(fits, offspring):
            ind.fitness.values = fit

        # Check if fitness is assigned properly
        if not offspring[0].fitness.values:
            print(f"Error: No fitness values assigned for generation {gen}")
            break

        # Print fitness of best individual in this generation
        best_individual = tools.selBest(population, k=1)[0]
        print(f"Generation {gen}: Best fitness (Residual) {best_individual.fitness.values[0]}")

        population = toolbox.select(offspring, k=len(population))

        best_per_gen.append([best_individual.fitness.values[0]] + list(best_individual))

    # Save fitness results
    columns = ["Fitness", "Gd1", "Gd2", "Gr", "e12", "e21", "epsilon1", "epsilon2", "G", "I_leak"]
    best_df = pd.DataFrame(best_per_gen, columns=columns)
    best_df.to_csv('best_fitness_per_generation.csv', index_label="Generation")

    # Final optimization and plotting
    best_params = tools.selBest(population, k=1)[0]
    print(f"Best individual: {best_params}")
    print(f"Best fitness: {best_params.fitness.values[0]}")

    y0 = [0.55, 0.225, 0.225]
    y0 = fsolve(lambda y: fun(0, y, list(best_params) + [470, 6e-4, 0.77], light_fun), y0)

    # Ensure fsolve converged
    if not np.all(np.isfinite(y0)):
        print("Error: fsolve did not converge.")
        return

    t_span = [0, 600]
    t_eval = np.arange(0, 600 + 10, 10)
    sol = solve_ivp(fun, t_span, y0, t_eval=t_eval, args=(list(best_params) + [470, 6e-4, 0.77], light_fun))

    plt.plot(t_eval, sol.y[0], label='O1')
    plt.plot(t_eval, sol.y[1], label='C2')
    plt.plot(t_eval, sol.y[2], label='O2')
    plt.xlabel('Time (s)')
    plt.ylabel('Concentration')
    plt.legend()
    plt.savefig('optimised_plot.png')
    plt.show()

# Execute the GA
run_GA()

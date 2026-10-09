import os
import random
import warnings
import multiprocessing
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve, minimize
from deap import base, creator, tools, algorithms

# Local imports
from modelparams import setDefaultParameters_opsinmodel_ODE
from modeleq import opsinmodel_ODE

# ============================================================
# 1. SETUP & DATA LOADING
# ============================================================
RANDOM_SEED = 6
random.seed(RANDOM_SEED)
np.random.seed(RANDOM_SEED)
warnings.filterwarnings("ignore")

base_folder = r"C:\Users\aj609\PatchDataAnalysisToolbox\opsincode\pythonfiles\gafitting_allcurrents\New_fitting\Fitting_currents_together"
os.makedirs(base_folder, exist_ok=True)
os.chdir(base_folder)

# ------------------------------------------------------------
# 14-DIMENSIONAL PARAMETER CONFIGURATION
# 9 Shared Parameters (7 Kinetics + G + I_leak) 
# + (1 Trace-Specific Parameter * 5 Currents)
# ------------------------------------------------------------
columns = ["I1", "I2", "I3", "I4", "I5"]

shared_bounds = [
    (0.001, 1.5),   # Gd1 
    (0.001, 1.5),   # Gd2 
    (1e-6,  0.05),  # Gr  
    (0.001, 0.5),   # e12 
    (0.001, 1.0),   # e21 
    (0.01,  2.0),   # epsilon1
    (0.01,  2.0),   # epsilon2
    (0.05,  15.0),  # G (Conductance) - Now Shared
    (0.0,   20.0)   # I_leak - Now Shared
]
shared_names = ["Gd1", "Gd2", "Gr", "e12", "e21", "epsilon1", "epsilon2", "G", "I_leak"]

indiv_bounds = [
    (0.0,   50.0),  # Irradiance - Only individual parameter
]
indiv_names = ["Irradiance"]

# Build Global Parameter Lists
param_bounds = list(shared_bounds)
param_names = list(shared_names)

for col in columns:
    param_bounds.extend(indiv_bounds)
    param_names.extend([f"{n}_{col}" for n in indiv_names])

# ------------------------------------------------------------
# Load Data
# ------------------------------------------------------------
Data = pd.read_csv("voltage_trace_chr2_cell6.csv")
Data.columns = Data.columns.str.strip()

t_data_full = Data["t"].values.astype(float)
t_fit = np.linspace(t_data_full[0], t_data_full[-1], 1500)

experimental_fit, normalisation = {}, {}
for col in columns:
    y = np.interp(t_fit, t_data_full, Data[col].values)
    experimental_fit[col] = y
    normalisation[col] = max(np.std(y), 1e-12)

# Smoothed Custom Weights
weights_fit = np.ones_like(t_fit)
weights_fit[(t_fit >= 98) & (t_fit <= 140)] = 40.0   # Entire peak and initial rapid decay
weights_fit[(t_fit > 140) & (t_fit <= 500)] = 5.0    # Steady-state plateau
weights_fit[(t_fit > 500) & (t_fit <= 550)] = 15.0   # Light-off deactivation

def light_fun(t, irradiance):
    return irradiance if 100 <= t <= 500 else 0.0

# ============================================================
# 2. MODEL SIMULATION & GLOBAL FITNESS
# ============================================================
def simulate_model(pars, t_eval):
    # pars expects exactly 10 parameters to run the ODE for a specific trace:
    # [Gd1, Gd2, Gr, e12, e21, epsilon1, epsilon2, Irradiance, G, I_leak]
    model_varpars = list(pars[:8]) + [470, 6e-4, 0.77]
    y_initial = np.array([0.55, 0.225, 0.225])
    
    try:
        y0 = fsolve(lambda y: opsinmodel_ODE(0, y, model_varpars, light_fun), y_initial, xtol=1e-5)
        if not np.all(np.isfinite(y0)): return None
        
        sol = solve_ivp(
            opsinmodel_ODE, [t_eval[0], t_eval[-1]], y0, t_eval=t_eval,
            args=(model_varpars, light_fun), method="LSODA", 
            max_step=0.5, rtol=1e-5, atol=1e-7
        )
        if not sol.success: return None
        
        I_opsin = pars[8] * (sol.y[0, :] + 0.05 * sol.y[2, :]) * (-80.0 - 0) - pars[9]
        return I_opsin if np.all(np.isfinite(I_opsin)) else None
    except Exception:
        return None

def computeGlobalResiduals(pars):
    # Check all 14 parameter bounds
    for val, (low, high) in zip(pars, param_bounds):
        if not (low <= val <= high): return (1e6,)

    total_fitness = 0.0
    
    # Extract shared parameters
    shared_kinetics = pars[0:7]
    G = pars[7]
    I_leak = pars[8]

    # Evaluate fitness for all 5 traces simultaneously using the shared parameters
    for i, col in enumerate(columns):
        # Irradiance is the only parameter changing per trace, starting at index 9
        irr = pars[9 + i]
        
        # Combine parameters to match what simulate_model expects
        current_pars = list(shared_kinetics) + [irr, G, I_leak]
        
        I_model = simulate_model(current_pars, t_fit)
        if I_model is None: return (1e6,)

        residual = (I_model - experimental_fit[col]) / normalisation[col]
        fitness = np.mean(weights_fit * residual**2)
        total_fitness += fitness

    # Return the average error across all traces as a single objective score
    return (total_fitness / len(columns),)

# ============================================================
# 3. DEAP GA SETUP
# ============================================================
if not hasattr(creator, "FitnessMin"):
    creator.create("FitnessMin", base.Fitness, weights=(-1.0,))
if not hasattr(creator, "Individual"):
    creator.create("Individual", list, fitness=creator.FitnessMin)

toolbox = base.Toolbox()
for i, (low, high) in enumerate(param_bounds):
    toolbox.register(f"attr_{i}", random.uniform, low, high)

# Dynamically gather all attributes for initialization
attrs = [getattr(toolbox, f"attr_{i}") for i in range(len(param_bounds))]
toolbox.register("individual", tools.initCycle, creator.Individual, attrs, n=1)
toolbox.register("population", tools.initRepeat, list, toolbox.individual)

def checkBounds(min_vals, max_vals):
    def decorator(func):
        def wrapper(*args, **kargs):
            offspring = func(*args, **kargs)
            for child in offspring:
                for i in range(len(child)):
                    child[i] = np.clip(child[i], min_vals[i], max_vals[i])
            return offspring
        return wrapper
    return decorator

mins, maxs = zip(*param_bounds)
toolbox.register("mate", tools.cxBlend, alpha=0.5)
toolbox.register("mutate", tools.mutGaussian, mu=0, sigma=0.15, indpb=0.2)
toolbox.decorate("mate", checkBounds(mins, maxs))
toolbox.decorate("mutate", checkBounds(mins, maxs))
toolbox.register("select", tools.selTournament, tournsize=3)

# ============================================================
# 4. OPTIMISATION LOOP (GA + LOCAL POLISH)
# ============================================================
def optimize_global(seed):
    print(f"\n{'='*50}\nOPTIMISING 9 MAIN PARAMETERS GLOBALLY\n{'='*50}")
    random.seed(seed)
    np.random.seed(seed)
    
    toolbox.register("evaluate", computeGlobalResiduals)
    pop = toolbox.population(n=150)
    hof = tools.HallOfFame(1)
    
    stats = tools.Statistics()
    stats.register("min_fit", lambda p: np.min([ind.fitness.values[0] for ind in p if ind.fitness.valid]))
    stats.register("best_ind", lambda p: list(tools.selBest(p, 1)[0])) 
    
    # Run GA 
    pop, logbook = algorithms.eaSimple(pop, toolbox, cxpb=0.85, mutpb=0.25, ngen=150, stats=stats, halloffame=hof, verbose=True)
    best_ga_params = hof[0]
    
    print(f"Refining Global Model with Local Optimizer (Nelder-Mead)...")
    
    res = minimize(
        lambda p: computeGlobalResiduals(p)[0], 
        best_ga_params, 
        bounds=param_bounds, 
        method='Nelder-Mead',
        options={'maxiter': 5000, 'maxfev': 10000, 'xatol': 1e-4, 'fatol': 1e-4}
    )
    
    best_final_params = res.x
    print(f"Final Polished Global Fitness: {res.fun:.8f}")
    
    # Extract evolutionary history
    history_records = []
    for record in logbook:
        row = {
            "Generation": record["gen"],
            "Fitness": record["min_fit"]
        }
        for name, val in zip(param_names, record["best_ind"]):
            row[name] = val
        history_records.append(row)
    
    return best_final_params, res.fun, history_records

# ============================================================
# 5. MAIN EXECUTION, DATA EXPORT & PLOTTING
# ============================================================
if __name__ == "__main__":
    pool = multiprocessing.Pool(processes=max(1, multiprocessing.cpu_count() - 1))
    toolbox.register("map", pool.map)

    # Trigger the Global Fit Strategy
    best_params, final_fitness, gen_history = optimize_global(RANDOM_SEED)

    pool.close()
    pool.join()

    # --------------------------------------------------------
    # Export Data to CSV
    # --------------------------------------------------------
    print(f"\n{'='*50}\nSAVING DATA TO CSV\n{'='*50}")
    
    df_timeseries = pd.DataFrame({"Time": t_data_full})
    summary_records = []
    
    # 1. Add all the GA generations to the summary records
    for row in gen_history:
        row["RMSE_pA_Overall"] = np.nan 
        row["R_squared_Overall"] = np.nan 
        summary_records.append(row)

    # Evaluate the Final Global Parameters across all 5 traces
    shared_kinetics = best_params[0:7]
    G = best_params[7]
    I_leak = best_params[8]
    
    all_residuals = []
    total_ss_res = 0.0
    total_ss_tot = 0.0
    
    for i, col in enumerate(columns):
        irr = best_params[9 + i]
        
        current_pars = list(shared_kinetics) + [irr, G, I_leak]
        
        I_model = simulate_model(current_pars, t_data_full)
        I_exp = Data[col].values
        residuals = I_exp - I_model 
        all_residuals.extend(residuals)
        
        ss_res = np.sum(residuals**2)
        ss_tot = np.sum((I_exp - np.mean(I_exp))**2)
        total_ss_res += ss_res
        total_ss_tot += ss_tot
        
        # Append to timeseries dataframe
        df_timeseries[f"Exp_{col}"] = I_exp
        df_timeseries[f"Fit_{col}"] = I_model
        df_timeseries[f"Res_{col}"] = residuals
        
    # Calculate global metrics
    overall_rmse = np.sqrt(np.mean(np.array(all_residuals)**2))
    overall_r2 = 1 - (total_ss_res / total_ss_tot) if total_ss_tot > 0 else np.nan
    
    # 2. Add the Final Polished row to the summary records
    final_row = {
        "Generation": "Final_Polished",
        "Fitness": final_fitness,
        "RMSE_pA_Overall": overall_rmse,
        "R_squared_Overall": overall_r2
    }
    for name, val in zip(param_names, best_params):
        final_row[name] = val
    summary_records.append(final_row)

    df_timeseries.to_csv("global_fitted_timeseries_data_14_parameters.csv", index=False)
    
    cols_order = ["Generation", "Fitness", "RMSE_pA_Overall", "R_squared_Overall"] + param_names
    pd.DataFrame(summary_records)[cols_order].to_csv("global_optimized_parameters_14_parameters.csv", index=False)
    
    print("Saved: global_fitted_timeseries_data_14_parameters.csv")
    print("Saved: global_optimized_parameters_14_parameters.csv")

    # --------------------------------------------------------
    # Final Plotting
    # --------------------------------------------------------
    plt.figure(figsize=(12, 8))
    colors = plt.cm.tab10(np.linspace(0, 1, len(columns)))
    
    for i, col in enumerate(columns):
        I_model = df_timeseries[f"Fit_{col}"]
        I_exp = df_timeseries[f"Exp_{col}"]
        
        plt.plot(t_data_full, I_exp, "o", color=colors[i], alpha=0.2, markersize=2, label=f"Exp {col}")
        plt.plot(t_data_full, I_model, "-", color=colors[i], linewidth=2.5, label=f"Fit {col}")
        
    plt.xlabel("Time (ms)", fontsize=14)
    plt.ylabel("Photocurrent (pA)", fontsize=14)
    
    plt.xlim(t_data_full.min() - 5, t_data_full.max() + 5)
    plt.legend(bbox_to_anchor=(1.02, 1), loc="upper left", borderaxespad=0.)
    plt.title("CHR2 Model Global Fit (9 Shared Params)", fontsize=16)
    plt.tight_layout()
    
    plt.savefig("global_voltage_trace_chr2_cell6_fitted_14_parameters.png", dpi=300, bbox_inches="tight")
    plt.show()
    
    print("\nFinished successfully.")
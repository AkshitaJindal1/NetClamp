# Genetic Algorithm for ChR2 Model Parameter Estimation

## Overview

This repository contains the Python-based Genetic Algorithm (GA) framework developed to estimate the biophysical parameters of a four-state Channelrhodopsin-2 (ChR2) model from experimentally recorded photocurrents.

The framework uses population-based optimisation to minimise the mismatch between model-predicted and measured photocurrent responses under different light intensities. By simultaneously fitting multiple photocurrent recordings, the algorithm estimates the kinetic and conductance parameters governing ChR2 activation, inactivation, and recovery.

This work forms part of **NetClamp**, a computational and experimental framework for developing model-based, all-optical closed-loop control of cellular electrical activity.

## ChR2 Model

The framework employs a four-state kinetic model of ChR2, comprising two closed states (C₁ and C₂) and two open states (O₁ and O₂). Light-dependent transitions between these states determine channel activation and photocurrent generation.

The model incorporates parameters describing channel transition rates, light sensitivity, maximal conductance, and baseline leak current.

## Genetic Algorithm

The genetic algorithm was implemented in Python 3.12 using the DEAP (Distributed Evolutionary Algorithms in Python) library (Fortin et al., 2012). DEAP provides a flexible framework for evolutionary optimisation, supporting population initialisation, fitness evaluation, selection, crossover, and mutation.

The framework integrates DEAP-based optimisation with numerical simulations of the four-state ChR2 model to estimate kinetic and conductance parameters by minimising the mismatch between experimentally recorded and model-predicted photocurrents.

The resulting parameter estimates provide a quantitative description of ChR2 photocurrent dynamics under different illumination conditions and support subsequent model-based light-control calculations.

# Dependencies

Python 3.12

DEAP – Genetic algorithm implementation

Additional numerical and data-processing libraries as required by the scripts

# Installation

pip install deap

# Reference

Fortin, F.-A., De Rainville, F.-M., Gardner, M.-A., Parizeau, M. and Gagné, C. (2012) ‘DEAP: Evolutionary algorithms made easy’, Journal of Machine Learning Research, 13(70), pp. 2171–2175. Available at: https://jmlr.org/papers/v13/fortin12a.html.

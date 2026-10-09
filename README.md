# NetClamp: Towards All-Optical Closed-Loop Control of Cellular Electrical Activity

## Overview

**NetClamp** is a experimental framework designed to enable all-optical, closed-loop control of cellular electrical activity through the integration of mathematical modelling, optical activity measurements, and spatially targeted optogenetic stimulation.

The framework aims to establish programmable interactions between electrically excitable cells by using model-derived light stimulation to reproduce defined electrophysiological responses. By combining optical recording and stimulation with computational models, NetClamp seeks to provide a scalable alternative to electrode-based approaches like dynamic clamp for investigating and controlling cellular electrical dynamics.

This repository contains the computational tools and analysis scripts developed as part of the NetClamp project, covering model parameterisation, optogenetic actuator characterisation, light waveform optimisation, and experimental validation.

## Key Components

The repository covers the following areas:

1. **Biophysical Modelling and Data Assimilation:** Estimation of ion-channel parameters from electrophysiological recordings using conductance-based models and optimisation techniques, including IPOPT.

2. **Optogenetic Actuator Modelling:** Mathematical characterisation of Channelrhodopsin-2 (ChR2) and enhanced Natronomonas pharaonis halorhodopsin (eNpHR3.0) to describe light-induced membrane currents.

3. **Parameter Inference:** Application of genetic algorithms and Bayesian inference to estimate model parameters from experimentally recorded photocurrents.

4. **Model-Based Optimal Control:** Computation of light stimulation waveforms designed to reproduce target excitatory and inhibitory postsynaptic current (EPSC/IPSC) profiles.

5. **Experimental Validation:** Analysis and validation of model-predicted photocurrent responses using real-time optical stimulation and patch-clamp electrophysiology.

6. **Optical Activity Detection and Spatial Stimulation:** Development and evaluation of optical recording, event detection, and patterned illumination methods to support future closed-loop integration.

## Research Objective

The long-term objective of NetClamp is to establish a fully optical feedback system capable of monitoring cellular electrical activity and delivering model-derived, spatially targeted stimulation in real time. Such a system could enable the investigation of programmable cellular interactions and the manipulation of electrical activity across networks of excitable cells.

## Current Development Status

The computational models, parameter estimation methods, optogenetic control strategies, and selected optical components have been developed and evaluated independently. Experimental validation has demonstrated the feasibility of model-based photocurrent control in GH3 cells. However, integration and validation of the complete all-optical closed-loop NetClamp system remain areas for future development.

## Applications

The framework is intended to support research in optogenetics, computational neuroscience, electrophysiology, model-based control, and the study of cellular and network electrical dynamics.

## Research Context

NetClamp was developed as part of doctoral research, funded by BBSRC South West Biosciences Doctoral Training Partnership (SWBio DTP), University of Exeter, UK, for investigating the experimental and computational foundations of model-based optical control of cellular electrical activity.

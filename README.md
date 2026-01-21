# R package `marp` Application (work in progress)

This repo contains scripts applying marp to real-world data, pipelines to transform said data and code for a simulation study to evaluate the methodology deployed by marp.

> ⚠️ This is an active research repository.  
> Code, interfaces, and outputs may change as methods and experiments evolve.

## File structure

- `demonstrations/`  
  Demonstrating marp on fuller2 (airplane glass failure), CHB-MIT Scalp EEG Database (seizure) data. 

- `pipelines/`  
  Pipelines used to transform data into a format for marp.

- `simulations/`  
  Code for a simulation study to evaluate the methods in marp. 

- `results/`  
  Output directory for figures, tables, and intermediate files (not tracked in Git). Currently empty.
  
## Data sources

Seizure data: https://physionet.org/content/chbmit/1.0.0/

Airplane glass data: https://datarepository.wolframcloud.com/resources/Sample-Data-Airplane-Glass/

## Usage

Run scripts from the project root, for example:


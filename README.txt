# main scripts to run simulations, in F://enhancer_dynamics/model_v2
# backed up to github.com/KSheu/ksheu.library2/

update_model.m
# updates model according to equations manually written in chromatin_model.xlsx

do_simulation.m
# runs ODE simulation based on data input signaling trajectories

do_simulation_simTF.m
# run simulation based on simulated input signaling trajectories. A scale factor is used to put the simulated input on the same scale as the microscopy data
# generates the input amplitude scan subpanel

do_simulation_simTF_totalactivity.
# generates the input total activity scan subpanel



# folder "paramscan" generates panels that scan the parameters of the model itself

do_simulation_simTF_paramscan.m
# scan of Hill coefficient and Kd


do_simulation_simTF_paramscan_scanfactor.m
# scans reverse and forward cooperativity factors

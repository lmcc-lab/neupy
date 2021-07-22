'''
Write your preferences for each module of neupy here.
'''

#Filtering nubase data
Experiment_only = False                     # boolean, When searching nubase, you can choose to include only experimentally verified data or all. Default is False
set_ground_state = True                     # boolean, Set all decays to ground state, ignoring possible higher energies after first decay if True.

#Decay Chain
max_decay_chain_depth = 15                  # int, The maximum decay chain depth. Default is 6
ignore_energetically_possible = False       # boolean, Ignore any theoretical branches in a decay chain. Default is False. If False, calculates what the branch ratio would be, unless there are multiple theoretical branches with no values.
min_branch_ratio = 0                        # int, The minimum branch ratio for a branch to be considered in a decay chain. Set -1 for all (including 0) branches. Defaults to 0

#Bateman
simple = False                              # boolean, Use only simple chains, that is, linear with beta decays only. Defaults to False

#Graph times and reactor events
graph_start = -3                            # float, Hrs or power to Hrs (10^(graph_start))
graph_end = 14                               # float, Hrs or power to Hrs (10^(graph_end))
resolution = 1000                           # int
logspace = True                             # boolean, Sets the time sampling to log if True. Uses graph_start and graph_end hrs as exponents, 10^a to 10^b
offset_t = 0                                # float, Hrs. Sets the time offset for plotting.

#neutpy chain contributions
percent_of_final = 100                       # float, Decimal as percent - percent of full result
top_num = 2612                                # int, Use the top # of contributing chains - number of most contributing chains. 2612 total chains
pref_percent = True                         # Prefer percentage value as the chain contribution decider if True. Otherwise use top number of contributions.

if simple == True:
    simpleTitle = 'Simple'
else:
    simpleTitle = 'NotSimple'
##
import os
import matplotlib.pyplot as plt

from lib.plot import RiskMap, Plotter

###########################################
###########################################





if __name__ == "__main__":

    print(f"\nPlotting in {os.getcwd()}")


    # random_prob_0.250
    name_topology = "airports"   # america europe airports
    # evname = "EARL" # EARL MATTHEW KARL GONZALO mock2 ciaran


    Map = RiskMap(name_topology)
    # Plot = Plotter(name_topology)

    ##      1) Leaflet map powergrids
    # Map.plot_leaflet(evname)

    ##      1 bis) Leaflet map airports
    # Map.plot_leaflet_airports()


    # list_event = ['INGRID', 'IRENE', 'EARL', 'KATE', 'SANDY', 'NATE', 'ISAAC', 'PAULA', 'MATTHEW', 'JOAQUIN', 'BILL', 'KATIA', 'HERMINE', 'ALEX', 'TOMAS', 'CRISTOBAL', 'IDA', 'KARL', 'ARTHUR', 'GONZALO', 'BERTHA']
    # [Map.plot_leaflet(fpath_risk,name_topology,name) for name in list_event]


    ##      2) Parametric plot
    # Plot.plot_heatmap2d()


    ##      3) Plot Map of the network infrstructure
    Map.plot_network()

    ##      4) Plot total risk for every storm
    # Map.plot_total_risk()


    ## Show or save
    save = False
    if not save: #args.save:
        plt.show()
    else:
        for figname, fig in Map.figdict.items():
            print("Saving {}...".format(figname))
            fig.savefig(
                "Figures/{}.pdf".format(figname), format='pdf', bbox_inches='tight', 
                pad_inches=0.01, transparent=True
            )

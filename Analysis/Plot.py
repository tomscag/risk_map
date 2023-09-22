""" plot figures """
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.colors as mcolors
from haversine import haversine, Unit

import folium
import h3
import h3pandas

import os
import glob
import time
from selenium import webdriver
from selenium.webdriver.firefox.options import Options



# Set fonts 
sizetext = 20
plt.rcParams['xtick.major.pad'] = '10'
plt.rc('text', usetex=True)
# plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'], 'size': sizetext})

csfont = {'fontname':'Comic Sans MS'}
hfont = {'fontname':'Helvetica'}
plt.rcParams['text.usetex'] = True   # Render text with Latex




#############################################
#############################################
#############################################

class Plotter():
    def __init__(self) -> None:
        self.figdict = {}   # For saving



    # static methods for risk map

    @staticmethod
    def fragility_model_storm(dist,force,const=1.092e-3,dist0=10):
        """ 
        Return the probability of a physical damage
        at distance (dist) and with wind (force)
        """
        if dist > 1800:
            return 0
        else:
            return  const*((force/65)**(8.02))/(dist+dist0)**2

    @staticmethod
    def load_topology_geodata(name_topology):
        fpath = f"../Data/Processed/Topologies/{name_topology}/{name_topology}.nodelist"
        return pd.read_csv(fpath,delimiter=" ",index_col=0,names=["label","lng","lat","geoid"],usecols=[0,1,2,3])

    @staticmethod  
    def damage_if_nodes_fail(fpath_risk,type="oad"):

        """ 
            OUTPUT
                data: (dataframe)
                    Return the damage (Risk) when that node fails, computed as 1-LCC
                        node_label  lng lat Risk
        """

        if type.lower() == "motter":
            print("Computing intrinsic risk with Motter model")
            toler_index      = "1"   # "0.0"  Tolerance index

            data = pd.read_csv(fpath_risk, sep=" ", header=0, index_col=0)
            data = data[["lng","lat",toler_index]]
            data = data.rename(columns={toler_index:"Risk"})

        elif type.lower() == "oad":
            print(f"Loading OAD results in {fpath_risk}")
            data = Plotter.load_topology_geodata(name_topology)
            LCC = pd.read_csv(fpath_risk,delimiter="\t",index_col=0,names=["LCC"]).squeeze()
            data["Risk"] = 1 - LCC
            data.dropna(inplace=True)
        return data


    @staticmethod
    def load_data_storm(name_storm,name_topology):

        _dir  = f'../Data/Processed/Storms/{name_topology}/storm_'

        data_storm = pd.read_csv(
                    _dir+name_storm+'.csv', delimiter=' ',
                    names=["Latitude", "Longitude", "wmo_wind.x", "PDM"],
                    header=0,dtype=float)        

        return data_storm
    
    @staticmethod
    def load_path_topology(name_topology):
        '''
        '''
        if name_topology=="america":
            return  "../Data/Processed/Topologies/america/powergrid_north_america.el", \
                    "../Data/Processed/Topologies/america/america.nodelist"

        elif name_topology=="europe":
            return "../Data/Processed/Topologies/europe/powergrid_europe.el", \
                    "../Data/Processed/Topologies/europe/europe.nodelist"
        
        elif name_topology=="airports":
            return "../Data/Processed/Topologies/airports/airports_world.edgelist",\
                    None  # To do
        
        else:
            print("Topology not recognized \n EXIT")
            return
    

    @staticmethod
    def add_storm_to_map(map,evname,name_topology):

        weight = 5
        color  = "#FFFF00"

        # Load data of the storm
        data_storm     = Plotter.load_data_storm(evname,name_topology)
        pos_storm = [(lat,lon) for lon,lat in  zip(data_storm["Longitude"], data_storm["Latitude"])]
        folium.PolyLine(locations=pos_storm, color=color,weight=weight).add_to(map)
        for lat,lon,pdm in zip(data_storm["Latitude"],data_storm["Longitude"],data_storm["PDM"]):
            folium.vector_layers.CircleMarker(location=(lat,lon),radius=pdm/10,
                                            fill=True,fillcolor=color).add_to(map)

        return map

    @staticmethod
    def export_map(map,bounds,width=1500, height=844,filename="test",outdir="./",delay=2.0):
        """ Save map as png using Selenium """
        options = Options()
        options.add_argument('-headless')
        
        sw, ne = bounds

        map.fit_bounds([sw, ne])
        map.save(filename+'.html')

        
        urlfn='file://'+os.getcwd()+'/'+filename+'.html'


        browser = webdriver.Firefox(options=options)
        browser.minimize_window()
        browser.set_window_size(width, height)
        browser.get(urlfn)

        #Give the map tiles some time to load
        time.sleep(delay)
        browser.save_screenshot(outdir+filename+'.png')
        browser.quit()


    ##################################
    ### Parametric plot 

    @staticmethod
    def fname_var_R0(r0,name_topology,NUM_NODES,NUM_SAMPLES,MAX_TSTEP):

        FILEPATH    = f"/home/tomsc/Projects/RiskMap/Analysis/Output_OAD/Simulation/{name_topology}_nodes_{NUM_NODES}_samples_{NUM_SAMPLES}_maxtime_{MAX_TSTEP}/"
        return FILEPATH + f"{name_topology}_r0_{r0:.5f}.dat"
        # return FILEPATH + f"powergrid_r0_{r0:.5f}.dat"



    def plot_heatmap2d(self,name_topology="america",NUM_NODES=500,NUM_SAMPLES=10,MAX_TSTEP=50):
        
        fig, ax = plt.subplots(figsize=(3,3))
        size_ticksnumber = 18 #plt.tick_params(labelsize=size_ticksnumber)
        size_axeslabel = 22 #plt.xlabel("$m$",{'fontsize': size_axeslabel})
        size_legend = 22 #plt.legend(fancybox=True, shadow=True, ncol = 3, numpoints=1,loc = 'upper center', fontsize = size_legend)
        size_text = 20 #plt.text(-0.16,1.05, r'$(A)$', horizontalalignment='left', verticalalignment='center',transform=ax.transAxes, fontsize= size_text) 
        size_ticksnumber_inset = 20
        size_axeslabel_inset = 20



        FILEPATH    = f"/home/tomsc/Projects/RiskMap/Analysis/Output_OAD/Simulation/{name_topology}_nodes_{NUM_NODES}_samples_{NUM_SAMPLES}_maxtime_{MAX_TSTEP}/"
        FILELIST    = glob.glob(FILEPATH+"*.dat")


        r0_list = sorted([float(item.split(".dat")[0].split("_")[-1]) for item in FILELIST])

        with open(FILELIST[0],"r") as fp:
            data    = fp.readlines()[6:]
            r1_list = [float(item.split("\t")[0].split()[0]) for item in data] # r1_list

        S = []
        for r0 in r0_list:
            file = Plotter.fname_var_R0(r0,name_topology,NUM_NODES,NUM_SAMPLES,MAX_TSTEP)
            with open(file,"r") as fp:               
                data = fp.readlines()
                name_topology = data[0].split()[-1]
                dynp_pINI = float(data[4].split()[-1])
                data = data[6:]
                row = [float(item.split("\t")[1].split()[0]) for item in data] # solution
                # S.append(row[0:97])
                S.append(row)

        arr = np.transpose(np.array(S, dtype=float))


        plt.tick_params(labelsize=1*size_ticksnumber) #if written below cax, it doesnt work
        plt.plot([1.0/(1-dynp_pINI),0.0],[0.0,1/(dynp_pINI*(1-dynp_pINI))], color='#dd181f', linewidth=1.5)


        interpolation = "bilinear" # none bilinear bicubic hanning
        im = plt.imshow(arr, extent=[np.min(r0_list),np.max(r0_list),np.min(r1_list),np.max(r1_list)], 
                    origin='lower', cmap='viridis', alpha=0.9, aspect='auto', interpolation=interpolation)

        cbar = plt.colorbar(im, cax = fig.add_axes([0.95, 0.12, 0.03, 0.66]), shrink=0.99, pad = 0.07)
        # cbar.ax.set_ylabel(r'$S$', rotation=0, fontsize = 25, labelpad=15)
        cbar.set_ticks([0,0.2,0.4,0.6,0.8])
        cbar.ax.set_title(r'$S$', {'fontsize': size_axeslabel})
        # cbar.ax.tick_params(labelsize=0.8*size_axeslabel)
        
        ax.set_xlabel(r"$\mathcal{R}_0$", {'fontsize': size_axeslabel})
        ax.set_xlabel(r"$\mathcal{R}_0$", **hfont)
        ax.set_ylabel(r"$\mathcal{R}_1$", {'fontsize': size_axeslabel})

        # ax.set_xlim((0,max(r0_list)))
        # ax.set_ylim((0,max(r1_list)))
        
        # xticks = [1, 2.5, 4]
        # plt.xticks([0.5,1.0,1.5,2.0])
        
        plt.tick_params(labelsize=size_ticksnumber)

        
        # Save
        self.figdict[f'heat_map'] = fig 

        return

    @staticmethod
    def compute_risk_from_event(evname,name_topology,fpath_risk,type):
        '''
            Compute the cumulative risk for each storm's snapshot 
            using the selected fragility model
            Returns the risk aggregated by geographical region
        '''
        risk_intrinsic_nodes = Plotter.damage_if_nodes_fail(fpath_risk,type)

        # Load storm's data
        data_storm     = Plotter.load_data_storm(evname,name_topology)

        risk_intrinsic_nodes["Prob"] = 0
        for snapshot in data_storm.iterrows():
            latlon_storm = (snapshot[1]['Latitude'], snapshot[1]['Longitude'])
            # pdm_storm = snapshot[1]['PDM'] # Potential Damage Multiplier
            wind_strength = snapshot[1]['wmo_wind.x'] # Wind strength
            distances  = [haversine(latlon_storm,(Lat,Lon), unit=Unit.KILOMETERS) 
                        for Lat,Lon in zip(risk_intrinsic_nodes.lat,risk_intrinsic_nodes.lng)]
            risk_intrinsic_nodes["Prob"] += [Plotter.fragility_model_storm(dist,wind_strength) for dist in distances]


        
        # Compute risk as damage times the cumulative probability
        scalefact = 1e0
        risk_intrinsic_nodes["RiskTot"] = [A*B*scalefact for A,B in zip(risk_intrinsic_nodes["Prob"],risk_intrinsic_nodes["Risk"])]

        norm_factor = 0.0006864062438303901  # To normalize risk to the strongest event
        risk_aggregated = risk_intrinsic_nodes.groupby("geoid")["RiskTot"].sum()
        risk_aggregated.drop("Other",inplace=True)
        risk_aggregated = risk_aggregated/norm_factor

        return risk_aggregated        


    #######################################
    ## Risk map plot (leaflet)

    def plot_leaflet(self,fpath_risk,name_topology,evname= "EARL",type="oad"):

        fig, axes = plt.subplots(figsize=(3,3))
        width, height  = (2000,1500)   # Pixels
        bounds         = ([16.4, -95.00],[55.5, -87.00])   # (south_west, north_east)     

        crs =  "EPSG3857"   # coordinate reference systems   EPSG3857 (default)  EPSG4326
        # Create map
        map = folium.Map(tiles="cartodbpositron",location=[39.50, -98.35], 
                         zoom_start=4, zoom_control=False, crs=crs)

        state_geo = f"../Data/Processed/Topologies/{name_topology}/{name_topology}-counties.geojson"

        risk = Plotter.compute_risk_from_event(evname,name_topology,fpath_risk,type)

        folium.GeoJson(
            state_geo,
            style_function = lambda feature: {
                'fillColor': Plotter.my_color_function(feature,risk),
                'color': '#080808',       #border color for the color fills
                'weight': .75,            #how thick the border has to be
                'opacity':1,
                'fillOpacity':.85,
                # 'dashArray': '5, 3'  #dashed lines length,space between them
            }
        ).add_to(map)

        map = Plotter.add_storm_to_map(map,evname,name_topology)
        
        map.save(evname+'.html')




    def my_color_function(feature,risk):
        """ Maps low values to green and hugh values to red."""
        bounds = np.linspace(0,1,6)
        bounds = np.logspace(-6,0,6)
        try:
            value = risk[feature['properties']['GEOID']]   # 
            if  value >bounds[0] and value <= bounds[1]:
                return '#2b83ba'   # blue
            elif value >= bounds[1] and value <= bounds[2]:
                return '#abd9e9'   # light blue
            elif value >= bounds[2] and value <= bounds[3]:
                return '#ffffbf'   # yellow
            elif value >= bounds[3] and value <= bounds[4]:
                return '#fdae61'   # orange
            elif value > bounds[4]:
                return '#d7191c'   # red
            else:
                return '#ffffff' # blue 
        except:
            return '#ffffff' #'#ffffff'   #808080
        


    ##########################################    
    ### Plot network infrastructure as a map

    def plot_network(self,name_topology):
        "Plot the networks"
        # import networkx as nx

        path_edgelist, path_nodelist = Plotter.load_path_topology(name_topology)
        # G = nx.read_edgelist(path_edgelist)
        # G = G.subgraph([str(item) for item in df_PowerGrid.index])
        data_nodes = pd.read_csv(path_nodelist,sep=" ",index_col=0, usecols=[0,1,2],names=["label","lon","lat"])
        data_edges = pd.read_csv(path_edgelist,sep=" ",names=["node1","node2"])
        fig, ax = plt.subplots(figsize=(6, 6))

        crs =  "EPSG3857"   # coordinate reference systems   EPSG3857 (default)  EPSG4326
        # Create map
        map = folium.Map(tiles="cartodbpositron",location=[39.50, -98.35], 
                         zoom_start=4.5, zoom_control=False, crs=crs)


        # Add edges to map
        data_edges.apply( lambda edge:  folium.PolyLine( (
                (data_nodes.loc[edge.node1].lat,data_nodes.loc[edge.node1].lon),
                (data_nodes.loc[edge.node2].lat,data_nodes.loc[edge.node2].lon),
                ),
                color = "#636363"
                ).add_to(map)    
                ,axis=1)

        # Add nodes to map
        data_nodes.apply(lambda point: folium.CircleMarker(location=[point.lat, point.lon],
                        radius=1,color="#252525", opacity=0.75,
                        weight=5
                        ).add_to(map),axis=1)

        if name_topology == "america":
            list_event  = ['INGRID', 'IRENE', 'EARL', 'KATE', 'SANDY', 'NATE', 'ISAAC', 'PAULA', 'MATTHEW', 'JOAQUIN', 'BILL', 'KATIA', 'HERMINE', 'ALEX', 'TOMAS', 'CRISTOBAL', 'IDA', 'KARL', 'ARTHUR', 'GONZALO', 'BERTHA']
            for event in list_event:
                data_storm = Plotter.load_data_storm(event,name_topology)
                data_storm.apply(lambda point: folium.CircleMarker(location=[point.Latitude, point.Longitude],
                    radius=0.05*point["wmo_wind.x"],color="#a50f15", opacity=0.75,
                    weight=5
                    ).add_to(map),axis=1)
                folium.PolyLine(tuple((a,b) for a,b in zip(data_storm.Latitude, data_storm.Longitude)) ,
                                color = "#de2d26",
                                ).add_to(map)



        map.save(name_topology+'.html')



        # plt.xlim([-135,-55])
        # plt.ylim([10,75])
        # if savefig:
        #     plt.savefig(f"./src/power_grid.pdf",format="pdf", bbox_inches="tight")
        
############################################
############################################
############################################




if __name__ == "__main__":

    print(f"\nPlotting in {os.getcwd()}")

    Pjotr = Plotter()
    
    name_topology = "europe"
    evname = "mock1" # EARL MATTHEW KARL GONZALO mock2
    r0 = 8
    r1 = 0.3
    fpath_risk = f"./Output_OAD/{name_topology}_r0_{r0}_r1_{r1}_samples_10_maxtime_2000.dat"

    ## Leaflet map
    # Pjotr.plot_leaflet(fpath_risk,name_topology,evname)


    ## Riskmap plot
    # Pjotr.plot_us_riskmap(fpath_risk,name_topology,evname)

    # [Pjotr.plot_us_riskmap(fpath_risk,name_topology,name) for name in ["EARL","ARTHUR","IRENE","ISAAC"]]


    # Parametric plot
    # Pjotr.plot_heatmap2d(name_topology="europe",NUM_NODES=1467,NUM_SAMPLES=75,MAX_TSTEP=2000)


    # Plot Map of the network infrstructure
    Pjotr.plot_network(name_topology)

    ## Show or save
    save = True
    if not save: #args.save:
        plt.show()
    else:
        for figname, fig in Pjotr.figdict.items():
            print("Saving {}...".format(figname))
            fig.savefig(
                "Figures/{}.pdf".format(figname), format='pdf', bbox_inches='tight', 
                pad_inches=0.01, transparent=True
            )


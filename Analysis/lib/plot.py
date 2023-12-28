""" plot figures """
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.colors as mcolors
from haversine import haversine, Unit
from lib.misc import (load_data_stressor,
                      load_topology_parameters,
                      load_topology_geodata,
                      fragility_model_earthquake,
                      fragility_model_storm)
from lib.input import Inputs

import folium
# import h3
# import h3pandas

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

class RiskMap():
    def __init__(self,name_topology) -> None:
        self.figdict = {}   # For saving
        self.name_topology = name_topology


        Inputso = Inputs()  
        inputs_dct= getattr(Inputso, self.name_topology)
        self.num_nodes = inputs_dct["num_nodes"]
        self.num_samples = inputs_dct["num_samples"]
        self.max_tstep = inputs_dct["max_tstep"]
        self.path_edgelist = inputs_dct["path_edgelist"]
        self.path_nodelist = inputs_dct["path_nodelist"]
        self.path_risk = inputs_dct["path_risk"]
        self.evname  = inputs_dct["event"]
        self.r0  = inputs_dct["r0"]
        self.r1  = inputs_dct["r1"]

        # Figure parameters (riskmap)
        self.width  = 2500 
        self.height = 1800
        self.delay  = 3.0


    
    def damage_if_nodes_fail(self,fpath_risk,type="oad"):

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
            # print(f"Loading OAD results in {fpath_risk}")
            data = load_topology_geodata(self.name_topology)
            LCC = pd.read_csv(fpath_risk,delimiter="\t",index_col=0,names=["LCC"]).squeeze()
            data["Risk"] = 1 - LCC
            data.dropna(inplace=True)
        return data


    def compute_fraction_at_risk(self,risk):
        '''
            Return the fraction of the nodes inside region
            with risk at least very low and high
                [TO DO: Refactor this function]
        '''
        data = load_topology_geodata(self.name_topology)
        C = data.groupby("geoid")["lat"].count()
        D = risk.index.to_series().apply( lambda item: C[str(item)])
        df=pd.concat([risk,D],axis=1).dropna()
        bounds = np.logspace(-6,0,6)
        df_filt1 = df[df> bounds[0]].dropna()
        df_filt2 = df[df> bounds[3]].dropna()
        f1 = df_filt1['geoid'].sum()/len(data)
        f2 = df_filt2['geoid'].sum()/len(data)
        print(f"Fraction of nodes with risk at least very low: {f1:.3f}\n")
        print(f"Fraction of nodes with risk at least high: {f2:.3f}\n")

        return f1,f2




    
    def add_storm_to_map(self,map):

        weight = 5
        color  = "#a50f15"

        # Load data of the storm
        data_storm     = load_data_stressor(self.evname,self.name_topology)
        pos_storm = [(lat,lon) for lon,lat in  zip(data_storm["Longitude"], data_storm["Latitude"])]

        # Add trajectory
        folium.PolyLine(locations=pos_storm, color=color,weight=weight).add_to(map)

        for lat,lon,pdm in zip(data_storm["Latitude"],data_storm["Longitude"],data_storm["PDM"]):
            folium.vector_layers.CircleMarker(location=(lat,lon),radius=pdm/10,
                                            fill=True,fillcolor=color,color=color).add_to(map)

        return map


    def export_map(self,map,filename="./test"):
        """ Save map as png using Selenium """
        options = Options()
        options.add_argument('-headless')
        

        ## Save html map
        map.save(filename+'.html')

        
        urlfn='file://'+os.getcwd()+'/'+filename+'.html'


        browser = webdriver.Firefox(options=options)
        browser.minimize_window()
        browser.set_window_size(self.width, self.height)
        browser.get(urlfn)
        time.sleep(self.delay)
        png = browser.get_screenshot_as_png()
        #Give the map tiles some time to load
        
        # browser.save_screenshot(filename+'.png')
        browser.quit()

        # Crop image 
        from PIL import Image
        from io import BytesIO
        im = Image.open(BytesIO(png)) # uses PIL library to open image in memory
        width, height = im.size

        x     = 0.2
        left  = x*width
        upper = x*height
        right = (1-x)*width
        lower = (1-x)*height
        # box = (250, 250, 750, 750) # left, upper, right, and lower
        box   = (left,upper,right,lower)
        im = im.crop(box) # defines crop points
        im.save(f'{filename}.png')


    ##################################
    ### Parametric plot 

    





    
    def compute_risk_from_event(self,evname,type):
        '''
            Compute the cumulative risk for each storm's snapshot 
            using the selected fragility model
            Returns the risk aggregated by geographical region
        '''
        fpath_risk = self.path_risk
        risk_intrinsic_nodes = RiskMap.damage_if_nodes_fail(self,fpath_risk,type)

        # Load storm's data
        data_storm     = load_data_stressor(evname,self.name_topology)

        risk_intrinsic_nodes["Prob"] = 0
        if self.name_topology == "america" or self.name_topology == "europe":
            for snapshot in data_storm.iterrows():
                latlon_storm = (snapshot[1]['Latitude'], snapshot[1]['Longitude'])
                # pdm_storm = snapshot[1]['PDM'] # Potential Damage Multiplier
                wind_strength = snapshot[1]['wmo_wind.x'] # Wind strength
                distances  = [haversine(latlon_storm,(Lat,Lon), unit=Unit.KILOMETERS) 
                            for Lat,Lon in zip(risk_intrinsic_nodes.lat,risk_intrinsic_nodes.lng)]
                risk_intrinsic_nodes["Prob"] += [fragility_model_storm(dist,wind_strength) for dist in distances]
        elif self.name_topology == "airports":
            for snapshot in data_storm.iterrows():
                latlon_storm = (snapshot[1]['Latitude'], snapshot[1]['Longitude'])
                magnitude = snapshot[1]['magnitude'] # Magnitude quake
                distances  = [haversine(latlon_storm,(Lat,Lon), unit=Unit.KILOMETERS) 
                            for Lat,Lon in zip(risk_intrinsic_nodes.lat,risk_intrinsic_nodes.lng)]
                risk_intrinsic_nodes["Prob"] += [fragility_model_earthquake(dist,magnitude) for dist in distances]



        
        # Compute risk as damage times the cumulative probability
        scalefact = 1e0
        risk_intrinsic_nodes["RiskTot"] = [A*B*scalefact for A,B in zip(risk_intrinsic_nodes["Prob"],risk_intrinsic_nodes["Risk"])]

        if self.name_topology == "america" or self.name_topology == "europe":
            norm_factor = 0.0006864062438303901  # To normalize risk to the strongest event
        elif self.name_topology == "airports":
            norm_factor = 0.00001
        risk_aggregated = risk_intrinsic_nodes.groupby("geoid")["RiskTot"].sum()
        risk_aggregated.drop("Other",inplace=True,errors='ignore')
        risk_aggregated = risk_aggregated/norm_factor

        return risk_aggregated        


    #######################################
    ## Risk map plot (powergrids)

    def plot_leaflet(self,evname= "EARL",type="oad"):
        self.evname = evname
        print(f"Plotting {self.evname}...")

        _,_, fpath_risk = load_topology_parameters(self.name_topology)

        fig, axes = plt.subplots(figsize=(3,3))
        width, height  = (2000,1500)   # Pixels
        bounds         = ([16.4, -95.00],[55.5, -87.00])   # (south_west, north_east)     

        tiles = "https://{s}.basemaps.cartocdn.com/light_nolabels/{z}/{x}/{y}.png"
        attr  = "<a href=http://www.openstreetmap.org/copyright>OpenStreetMap</a>"

        crs =  "EPSG3857"   # coordinate reference systems   EPSG3857 (default)  EPSG4326
        # Create map
        map = folium.Map(tiles=tiles,attr=attr,location=[35.50, -95.35], # [39.50, -98.35]
                         zoom_start=4.5, zoom_control=False, crs=crs)

        state_geo = f"../Data/Processed/Topologies/{self.name_topology}/{self.name_topology}-counties.geojson"

        risk = RiskMap.compute_risk_from_event(self,self.evname,type="oad")

        # Compute fraction ad risk
        RiskMap.compute_fraction_at_risk(self,risk)

        folium.GeoJson(
            state_geo,
            style_function = lambda feature: {
                'fillColor': RiskMap.my_color_function(feature,risk),
                'color': '#080808',       #border color for the color fills
                'weight': .75,            #how thick the border has to be
                'opacity':1,
                'fillOpacity':.85,
                # 'dashArray': '5, 3'  #dashed lines length,space between them
            }
        ).add_to(map)

        map = RiskMap.add_storm_to_map(self,map)

        # Add name to the map
        from folium.features import DivIcon
        folium.map.Marker(
            [22.00, -125.00],
            icon=DivIcon(
                icon_size=(1000,200),
                icon_anchor=(0,0),
                html=f'<div style="font-size:70pt">{self.evname}</div>',
                )
            ).add_to(map)

        RiskMap.export_map(self,map,filename=f"./Figures/self.{evname}")

        # Close map
        plt.close()


    #######################################
    ## Risk map plot (airports)

    def plot_leaflet_airports(self,type="oad"):
        self.evname = "earthquakes"
        self.name_topology = "airports" # Overwrite in case of mistakes
        
        print(f"Plotting airports map...")
        _,_, fpath_risk = load_topology_parameters(self.name_topology)


        tiles = "https://{s}.basemaps.cartocdn.com/light_nolabels/{z}/{x}/{y}.png"
        attr  = "<a href=http://www.openstreetmap.org/copyright>OpenStreetMap</a>"

        crs =  "EPSG3857"   # coordinate reference systems   EPSG3857 (default)  EPSG4326
        # Create map
        map = folium.Map(tiles=tiles,attr=attr,location=[35.50, -95.35], # [39.50, -98.35]
                         zoom_start=4.5, zoom_control=False, crs=crs)

        state_geo = f"../Data/Processed/Topologies/{self.name_topology}/{self.name_topology}-world-grid.json"

        risk = RiskMap.compute_risk_from_event(self,self.evname,type)

        # Compute fraction ad risk
        # RiskMap.compute_fraction_at_risk(self,risk)


        # Draw airports on the map
        path_edgelist, path_nodelist,_ = load_topology_parameters(self.name_topology)
        data_nodes = pd.read_csv(path_nodelist,sep=" ",index_col=0, usecols=[0,1,2],names=["label","lon","lat"])
        data_edges = pd.read_csv(path_edgelist,sep=" ",names=["node1","node2"])
        data_nodes.apply(lambda point: folium.CircleMarker(location=[point.lat, point.lon],
                        radius=0.70,color="#252525", opacity=0.75,
                        weight=5
                        ).add_to(map),axis=1)
        

        # Assign risk (color) to each cell
        folium.GeoJson(
            state_geo,
            style_function = lambda feature: {
                'fillColor': RiskMap.my_color_function(feature,risk),
                'color': '#080808',       # border color 
                'weight': .75,            # Width border
                'opacity':0.2,            # Border opacity
                'fillOpacity':.35,        # Fill opacity
            }
        ).add_to(map)


        RiskMap.export_map(self,map,filename=f"./Figures/{self.evname}")

        # Close map
        plt.close()





    def my_color_function(feature,risk):
        """ Maps low values to green and hugh values to red."""
        bounds = np.linspace(0,1,6)
        bounds = np.logspace(-6,0,6)
        try:
            value = risk[feature['properties']['GEOID']]   # GEOID
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
                return '#ffffff' # white 
        except:
            return '#ffffff' #'#ffffff'   #808080
        

    ##########################################
    ### Plot total risk

    def plot_total_risk(self):

        fpath_risk = self.path_risk
        plt.rc('font', size=30)  # Set font size for scientific notation
        events = ['INGRID', 'IRENE', 'EARL', 'KATE', 'SANDY', 'NATE', 'ISAAC', 'PAULA', 'MATTHEW', 'JOAQUIN', 'BILL', 'KATIA', 'HERMINE', 'ALEX', 'TOMAS', 'CRISTOBAL', 'IDA', 'KARL', 'ARTHUR', 'GONZALO', 'BERTHA']
        total_risk = pd.Series(index=events, dtype = float)
        Na = 16167 # America power grid nodes
        Ne = 13844 # Europe power grid nodes
        fig, ax = plt.subplots(figsize=(10,12))
        for evname in events:
            risk = RiskMap.compute_risk_from_event(self,evname,type="oad")
            total_risk.loc[evname] = risk.sum()/Na

        # Risk of the European mock events
        total_risk["MEDICANE 1"] = 3.06/Ne
        total_risk["MEDICANE 2"] = 6.28/Ne
        
        total_risk = total_risk.sort_values(ascending=True)
        total_risk = total_risk[-12:]

        plt.barh(total_risk.index,total_risk,facecolor = "#7570b3",alpha=0.75,edgecolor="#303030")
        plt.xticks(rotation = 0,fontsize=35)
        plt.yticks(fontsize=20)
        
        ax.set_xlabel(r"$R$", {'fontsize': 52})
        ax.set_ylabel(r"Hurricanes", {'fontsize': 40})
        ax.ticklabel_format(axis='x',style='scientific',scilimits=[-3,6])


        # Save
        self.figdict[f'total_risk'] = fig 


    ##########################################    
    ### Plot network infrastructure as a map

    def plot_network(self):
        "Plot the networks"
        tiles = "https://{s}.basemaps.cartocdn.com/light_nolabels/{z}/{x}/{y}.png"
        attr  = "<a href=http://www.openstreetmap.org/copyright>OpenStreetMap</a>"

        path_edgelist, path_nodelist = load_topology_parameters(self.name_topology)
        # G = nx.read_edgelist(path_edgelist)
        # G = G.subgraph([str(item) for item in df_PowerGrid.index])
        data_nodes = pd.read_csv(path_nodelist,sep=" ",index_col=0, usecols=[0,1,2],names=["label","lon","lat"])
        data_edges = pd.read_csv(path_edgelist,sep=" ",names=["node1","node2"])
        fig, ax = plt.subplots(figsize=(6, 6))

        crs =  "EPSG3857"   # coordinate reference systems   EPSG3857 (default)  EPSG4326
        # Create map
        map = folium.Map(tiles=tiles, attr = attr,location=[39.50, -98.35], 
                         zoom_start=4.5, zoom_control=False, crs=crs)


        # Add edges to map
        data_edges.apply( lambda edge:  folium.PolyLine( (
                (data_nodes.loc[edge.node1].lat,data_nodes.loc[edge.node1].lon),
                (data_nodes.loc[edge.node2].lat,data_nodes.loc[edge.node2].lon),
                ),
                weight=3,  # default 3, use 0.75 for airports
                color = "#636363"
                ).add_to(map)    
                ,axis=1)

        # Add nodes to map
        data_nodes.apply(lambda point: folium.CircleMarker(location=[point.lat, point.lon],
                        radius=1,color="#252525", opacity=0.75,
                        weight=5
                        ).add_to(map),axis=1)

        if self.name_topology == "america":
            list_event  = ['INGRID', 'IRENE', 'EARL', 'KATE', 'SANDY', 'NATE', 'ISAAC', 'PAULA', 'MATTHEW', 'JOAQUIN', 'BILL', 'KATIA', 'HERMINE', 'ALEX', 'TOMAS', 'CRISTOBAL', 'IDA', 'KARL', 'ARTHUR', 'GONZALO', 'BERTHA']
            for event in list_event:
                data_storm = RiskMap.load_data_stressor(event,self.name_topology)
                data_storm.apply(lambda point: folium.CircleMarker(location=[point.Latitude, point.Longitude],
                    radius=0.05*point["wmo_wind.x"],color="#a50f15", opacity=0.75,
                    weight=5
                    ).add_to(map),axis=1)
                folium.PolyLine(tuple((a,b) for a,b in zip(data_storm.Latitude, data_storm.Longitude)) ,
                                color = "#de2d26",
                                ).add_to(map)



        map.save(self.name_topology+'.html')



        # plt.xlim([-135,-55])
        # plt.ylim([10,75])
        # if savefig:
        #     plt.savefig(f"./src/power_grid.pdf",format="pdf", bbox_inches="tight")
        
############################################
############################################
############################################



class Plotter():
    def __init__(self,name_topology) -> None:
        self.figdict = {}   # For saving
        self.name_topology = name_topology


        Inputso = Inputs()  
        inputs_dct= getattr(Inputso, self.name_topology)
        self.num_nodes = inputs_dct["num_nodes"]
        self.num_samples = inputs_dct["num_samples"]
        self.max_tstep = inputs_dct["max_tstep"]
        self.path_edgelist = inputs_dct["path_edgelist"]
        self.path_nodelist = inputs_dct["path_nodelist"]
        self.path_risk = inputs_dct["path_risk"]
        self.evname  = inputs_dct["event"]
        self.r0  = inputs_dct["r0"]
        self.r1  = inputs_dct["r1"]

        # Figure parameters (riskmap)




    def fname_var_R0(self,r0):

        FILEPATH    = f"/home/tomsc/gdrive/Projects/RiskMap/Analysis/Output_OAD/Simulation/{self.name_topology}_nodes_{self.num_nodes}_samples_{self.num_samples}_maxtime_{self.max_tstep}/"
        return FILEPATH + f"{self.name_topology}_r0_{r0:.5f}.dat"
        # return FILEPATH + f"powergrid_r0_{r0:.5f}.dat"



    def plot_heatmap2d(self):

        fig, ax = plt.subplots(figsize=(3,3))
        size_ticksnumber = 18 #plt.tick_params(labelsize=size_ticksnumber)
        size_axeslabel = 22 #plt.xlabel("$m$",{'fontsize': size_axeslabel})
        size_legend = 22 #plt.legend(fancybox=True, shadow=True, ncol = 3, numpoints=1,loc = 'upper center', fontsize = size_legend)
        size_text = 20 #plt.text(-0.16,1.05, r'$(A)$', horizontalalignment='left', verticalalignment='center',transform=ax.transAxes, fontsize= size_text) 
        size_ticksnumber_inset = 20
        size_axeslabel_inset = 20



        FILEPATH    = f"/home/tomsc/gdrive/Projects/RiskMap/Analysis/Output_OAD/Simulation/{self.name_topology}_nodes_{self.num_nodes}_samples_{self.num_samples}_maxtime_{self.max_tstep}/"
        FILELIST    = glob.glob(FILEPATH+"*.dat")


        r0_list = sorted([float(item.split(".dat")[0].split("_")[-1]) for item in FILELIST])

        with open(FILELIST[0],"r") as fp:
            data    = fp.readlines()[6:]
            r1_list = [float(item.split("\t")[0].split()[0]) for item in data] # r1_list

        S = []
        for r0 in r0_list:
            file = RiskMap.fname_var_R0(self,r0)
            with open(file,"r") as fp:               
                data = fp.readlines()
                self.name_topology = data[0].split()[-1]
                dynp_pINI = float(data[4].split()[-1])
                data = data[6:]
                row = [float(item.split("\t")[1].split()[0]) for item in data] # solution
                # S.append(row[0:97])
                S.append(row)

        arr = np.transpose(np.array(S, dtype=float))


        plt.tick_params(labelsize=1*size_ticksnumber) #if written below cax, it doesnt work
        plt.plot([1.0/(1-dynp_pINI),0.0],[0.0,1/(dynp_pINI*(1-dynp_pINI))], color='#dd181f', linewidth=1.5)


        interpolation = "none" # none bilinear bicubic hanning
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



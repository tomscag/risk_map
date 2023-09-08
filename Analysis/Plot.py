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

    @staticmethod
    def fragility_model_storm(dist,force,const=1.092e-3,dist0=10):
        """ 
        Return the probability of a physical damage
        at distance (dist) and with wind (force)
        """
        return  const*((force/65)**(8.02))/(dist+dist0)**2

    @staticmethod
    def load_topology_geodata(name_topology):
        fpath = f"../Data/Processed/Topologies/{name_topology}/{name_topology}.nodelist"
        return pd.read_csv(fpath,delimiter=" ",index_col=0,names=["label","lng","lat"])

    @staticmethod  
    def load_risk_intrinsic_nodes(fpath_risk,type="oad"):

        """ 
            OUTPUT
                data: (dataframe)
                    Load dataframe in the format
                        node_label  lng lat
        """

        if type.lower() == "motter":
            print("Computing intrinsic risk with Motter model")
            toler_index      = "1"   # "0.0"  Tolerance index

            risk_intrinsic_nodes = pd.read_csv(fpath_risk, sep=" ", header=0, index_col=0)
            risk_intrinsic_nodes = risk_intrinsic_nodes[["lng","lat",toler_index]]
            risk_intrinsic_nodes = risk_intrinsic_nodes.rename(columns={toler_index:"Risk"})

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
    

    # static methods for risk map

    @staticmethod
    def aggregate_polygons(data_nodes,resol):
        """ Aggregate nodes into H3 polygons 
                data_nodes must have lng and lat columns
        """
        df_h3      = data_nodes.h3.geo_to_h3(resol)    # Assign points to hexagons

        # Aggregate data inside hexagons by summing (the risk is cumulative)
        data_nodes = df_h3.drop(columns=['lng', 'lat']).dropna(axis=0).groupby('h3_0'+str(resol)).sum() 
        gdf_h3     = data_nodes.h3.h3_to_geo_boundary()   # Extract vertices of each hexagon (geometry)

        df       = df_h3.drop(columns=['lng', 'lat']).dropna(axis=0).groupby('h3_0'+str(resol)).mean()
        return gdf_h3


    @staticmethod
    def create_map(gdf_h3,normalize=False):

        """ Create folium map """

        
        map = folium.Map(location=[39.50, -98.35], zoom_start=4)

        
        if normalize:
            dz = gdf_h3["Risk"]  
            dz = dz.apply(lambda x: x/max(dz))
        else:
            dz = gdf_h3["Risk"]
        fill_opacity = 0.75
        edges_color  = "#555555"
        png_enabled  = True
        cmap = mcolors.LinearSegmentedColormap.from_list('mycmap', ['#444444', '#FF0000']) # FF0000 = Red
        
        colors = cmap(dz)
        hex_colors = [mcolors.to_hex(color) for color in colors]
        hexagons = gdf_h3.index

        count = 0
        # Draw hexagons in the map
        for hexagon in gdf_h3.iloc:
            vertices = list(h3.h3_to_geo_boundary(hexagon.name))
            vertices.append(vertices[0]) # Add the first vertex to the end to close the polygon
            folium.vector_layers.Polygon(locations=vertices, color=edges_color, fill=True, 
                                        fill_color=hex_colors[count], fill_opacity=fill_opacity,
                                        png_enabled=png_enabled).add_to(map)
            count = count + 1
        
        return map

    @staticmethod
    def add_storm_to_map(map,data_storm):

        weight = 5
        color  = "#FFFF00"

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




    ########################################
    ## Riskmap plot 

    def plot_us_riskmap(self,fpath_risk,name_topology,evname= "EARL",type="oad"):

        fig, axes = plt.subplots(figsize=(3,3))
        risk_intrinsic_nodes = Plotter.load_risk_intrinsic_nodes(fpath_risk,type)
        
        # Figure parameters
        eventpath      = "./prob_failure_storms/"
        width, height  = (2000,1500)   # Pixels
        bounds         = ([16.4, -95.00],[55.5, -87.00])   # (south_west, north_east)        
        resol          = 3
        scalefact      = 6*10e4  # Motter:7.5  OAD:5

        # Load storm's data
        data_storm     = Plotter.load_data_storm(evname,name_topology)

        # Select the snapshot with strongest winds 
        latlng_storm_max = (data_storm.iloc[data_storm["wmo_wind.x"].idxmax()].Latitude, 
                            data_storm.iloc[data_storm["wmo_wind.x"].idxmax()].Longitude)
        force_storm_max  = data_storm.iloc[ data_storm["wmo_wind.x"].idxmax()]["wmo_wind.x"]
        distances      = [haversine(latlng_storm_max,(Lat,Lon), unit=Unit.KILOMETERS) 
                         for Lat,Lon in zip(risk_intrinsic_nodes.lat,risk_intrinsic_nodes.lng)]
        
        # Compute probability of failure of every node
        prob_failure_nodes = [Plotter.fragility_model_storm(dist,force_storm_max) for dist in distances]

        # Compute risk = risk_intrinsic_nodes * prob_failure_nodes
        risk_event_nodes = pd.Series([A*B*scalefact for A,B in  zip(risk_intrinsic_nodes.Risk, prob_failure_nodes)],name="Risk", index=risk_intrinsic_nodes.index)
        risk_event_nodes = pd.concat([risk_intrinsic_nodes[["lng","lat"]],risk_event_nodes],axis=1,join="outer")
        
        gdf_h3     = Plotter.aggregate_polygons(risk_event_nodes,resol)


        map = Plotter.create_map(gdf_h3,False)
        map = Plotter.add_storm_to_map(map,data_storm)
        Plotter.export_map(map,bounds, width, height,evname,"./",delay=2.0)

        # Save
        self.figdict[f'risk_map_{evname}'] = fig 


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
        
        # plt.text(0.025,0.06, r'{\fontfamily{ptm}\selectfont Non-percolating}', color='white', horizontalalignment='left',verticalalignment='center',fontsize=size_text,transform=ax.transAxes)
        # plt.text(0.95,0.9, r'{\fontfamily{ptm}\selectfont Percolating}', color='black', horizontalalignment='right',verticalalignment='center',fontsize=size_text,transform=ax.transAxes)
        plt.tick_params(labelsize=size_ticksnumber)
        
        # plt.text(-0.25, 1, r"$(a)$", horizontalalignment='left',verticalalignment='center',transform=ax.transAxes, size=0.9*size_axeslabel)
        # plt.title(f"{name_topology} nodes: {NUM_NODES}")
        
        # Save
        self.figdict[f'heat_map'] = fig 

        return





############################################
############################################
############################################




if __name__ == "__main__":

    print(f"\nPlotting in {os.getcwd()}")

    Pjotr = Plotter()
    
    ## Riskmap plot
    name_topology = "europe"
    evname = "mock2" # "EARL"
    r0 = 6
    r1 = 0.3
    fpath_risk = f"./Output_OAD/{name_topology}_r0_{r0}_r1_{r1}_samples_10_maxtime_2000.dat"
    Pjotr.plot_us_riskmap(fpath_risk,name_topology,evname)

    # [Pjotr.plot_us_riskmap(fpath_risk,name_topology,name) for name in ["EARL","ARTHUR","IRENE","ISAAC"]]


    # Parametric plot
    # Pjotr.plot_heatmap2d(name_topology="europe",NUM_NODES=1467,NUM_SAMPLES=75,MAX_TSTEP=2000)



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


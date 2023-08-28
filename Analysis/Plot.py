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
import time
from selenium import webdriver
from selenium.webdriver.firefox.options import Options

#############################################
#############################################
#############################################

class Plotter():
    def __init__(self) -> None:
        pass

    @staticmethod
    def fragility_model_storm(dist,force,const=1.092e-3,dist0=10):
        """ 
        Return the probability of a physical damage
        at distance (dist) and with wind (force)
        """
        return  const*((force/65)**(8.02))/(dist+dist0)**2

    

    @staticmethod  
    def load_risk_intrinsic_nodes(type="oad"):

        """ 
            This is the risk associated to each node
            obtained with the Motter model 
        """

        if type.lower() == "motter":
            print("Computing intrinsic risk with Motter model")
            toler_index      = "1"   # "0.0"  Tolerance index
            fname_risk     = "./Analysis/Output_Motter/Intrinsic_Risk_LngLat.dat" 

            risk_intrinsic_nodes = pd.read_csv(fname_risk, sep=" ", header=0, index_col=0)
            risk_intrinsic_nodes = risk_intrinsic_nodes[["lng","lat",toler_index]]
            risk_intrinsic_nodes = risk_intrinsic_nodes.rename(columns={toler_index:"Risk"})

        elif type.lower() == "oad":

            respath = "./Analysis/Output_OAD/lb_10_mu_1_al_0.3"

            print(f"Loading OAD results in {respath}")
            fname_risk     = "./Analysis/Output_Motter/Intrinsic_Risk_LngLat.dat" 
            risk_intrinsic_nodes = pd.read_csv(fname_risk, sep=" ", header=0, index_col=0)
            risk_intrinsic_nodes = risk_intrinsic_nodes[["lng","lat"]]

            
            scaling = 0.927136   # Scaling factor to get the relative size of LCC
            import glob
            dct = {}
            nodeflist = glob.glob(respath+"/node_*.dat")
            for nodefile in nodeflist:
                node = int(nodefile.split("_")[-2])            # Format node_1234_gillespie.dat
                if node in risk_intrinsic_nodes.index:   
                    with open(nodefile, 'r') as file:
                        item = file.readline()
                        lcc_gillespie = sum([float(f) for f in item.split()])/len(item.split())
                        dct[node] = lcc_gillespie
            dct = {key: value/scaling for key,value in dct.items()}
            lcc_gillespie  = pd.Series(dct)
            risk_intrinsic_nodes["Risk"] = 1 - lcc_gillespie
            risk_intrinsic_nodes = risk_intrinsic_nodes.dropna()


        return risk_intrinsic_nodes



    @staticmethod
    def load_data_storm(name_storm):

        _dir  = './Data/Processed/Storms/storm_'

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

    def plot_us_riskmap(self,evname= "EARL"):

        risk_intrinsic_nodes = Plotter.load_risk_intrinsic_nodes()
        
        # Figure parameters
        eventpath      = "./Analysis/prob_failure_storms/"
        width, height  = (2000,1500)   # Pixels
        bounds         = ([16.4, -95.00],[55.5, -87.00])   # (south_west, north_east)        
        resol          = 3
        scalefact      = 6*10e4  # Motter:7.5  OAD:5

        # Load storm's data
        data_storm     = Plotter.load_data_storm(evname)

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





############################################
############################################
############################################




if __name__ == "__main__":

    print(f"\nPlotting in {os.getcwd()}")

    Pjotr = Plotter()
    
    # Riskmap plot
    Pjotr.plot_us_riskmap("EARL")

    # [Pjotr.plot_us_riskmap(name) for name in ["EARL","ARTHUR","IRENE","ISAAC"]]

#pip install geopandas shapely pandas numpy geopy requests
#Imports
from cmu_graphics import *
import geopandas as gpd
from shapely.geometry import MultiPolygon
from shapely.ops import transform
from shapely.geometry import Point
from shapely.geometry import Polygon
import pandas as pd
import json
from geopy.distance import geodesic
import random
import requests
import numpy as np

###########################################################################################
#CODE FROM https://medium.com/@osah.dilshan/how-to-extract-gps-data-from-images-using-python-9c09254bc80e
import os
from PIL import Image
from PIL.ExifTags import TAGS

image_path = "/Users/hans/Downloads/CMUImages/IMG_1030.JPG"

def get_coordinate(image_path):

    # Open the image file
    image = Image.open(image_path)

    exif = {}
    if image._getexif() is not None:
        for tag, value in image._getexif().items():
            if tag in TAGS:
                exif[TAGS[tag]] = value

    if "GPSInfo" in exif:
        gps_info = exif["GPSInfo"]

        def convert_to_degrees(value):
            """
            Helper function to convert the GPS coordinates stored in the EXIF to degrees in float format.

            Args:
                value (tuple): The GPS coordinate as a tuple (degrees, minutes, seconds)

            Returns:
                float: The coordinate in degrees
            """
            d = float(value[0])
            m = float(value[1])
            s = float(value[2])
            return d + (m / 60.0) + (s / 3600.0)

        lat = convert_to_degrees(gps_info[2])
        lon = convert_to_degrees(gps_info[4])
        lat_ref = gps_info[1]
        lon_ref = gps_info[3]

        if lat_ref != "N":
            lat = -lat
        if lon_ref != "E":
            lon = -lon

        geo_coordinate = (lon, lat)

        return geo_coordinate

#################################################
# Json File Reading
with open('.src/geojson-maps.json', 'r') as f:
    json_data = json.load(f)

data = []
for feature in json_data["features"]:
    properties = feature["properties"]

    if properties.get("adm0_a3") == "SDS":  # I HATE SOUTH SUDAN
        properties["adm0_a3"] = "SSD"

    if "adm0_a3" in properties and "name" in properties and "pop_est" in properties and "subregion" in properties and "income_grp" in properties:
        data.append({
            "adm0_a3": properties["adm0_a3"],
            "name": properties["name"],
            "pop_est": properties["pop_est"],
            "subregion": properties["subregion"],
            "income_grp": properties['income_grp']
        })

df = pd.DataFrame(data)
exceptions = ["GRL", "ISL"]

filtered_df = df[(df["pop_est"] >= 2000000) | (df["adm0_a3"].isin(exceptions))]

# GeoPandas file reading
geojson_path = '.src/datahub.geojson'
world_data = gpd.read_file(geojson_path)
filtered_world_data = world_data[world_data['ISO_A3'].isin(filtered_df['adm0_a3'])]
filtered_world_data = pd.merge(
    filtered_world_data,
    filtered_df,
    left_on='ISO_A3',
    right_on='adm0_a3',
    how='left'
)

###########################################################################################
# Start of ChatGPT supported segment ######################################################
def geo_to_screen(lon, lat, width, height):
    R = 6378137  
    lat_rad = np.radians(lat)
    lon_rad = np.radians(lon)
    mercator_x = R * lon_rad
    mercator_y = R * np.log(np.tan(np.pi / 4 + lat_rad / 2))

    screen_x = int((mercator_x + 20037508.34) * (width / 40075016.68))
    screen_y = int((20037508.34 - mercator_y) * (height / 40075016.68))
    return screen_x, screen_y


def screen_to_geo(screen_x, screen_y, width, height):

    R = 6378137  

    mercator_x = screen_x * (40075016.68 / width) - 20037508.34
    mercator_y = 20037508.34 - screen_y * (40075016.68 / height)

    lat_rad = np.arctan(np.sinh(mercator_y / R))
    lat = np.degrees(lat_rad)

    lon = np.degrees(mercator_x / R)

    lat = max(min(lat, 90), -90)

    return lon, lat


# End of ChatGPT supported segment ########################################################
###########################################################################################

#Generate Shapes and Helper Functions
country_shapes = {}
for _, row in filtered_world_data.iterrows():
    country_code = row['adm0_a3']
    
    geom = row['geometry']
    if isinstance(geom, MultiPolygon):
        polygons = list(geom.geoms)
    else:
        polygons = [geom]

    screen_polygons = []
    for polygon in polygons:
        simplified_polygon = polygon.simplify(0.3) ###IMPORTANT FOR FASTER LOADING SPEED: simplify the higher float the more simple
        screen_coords = [geo_to_screen(x, y, 1200, 800) for x, y in simplified_polygon.exterior.coords] #width & height of app
        screen_polygons.append(screen_coords)
        country_shapes[country_code] = screen_polygons



def get_country_neighbors():
    response = requests.get("https://restcountries.com/v3.1/all")
    response.raise_for_status()
    
    countries_data = response.json()
    
    country_neighbors = {}
    
    for country in countries_data:
                
        # Each country's code
        country_code = country.get("cca3")
            
            # Neighboring countries' codes
        neighbors = list(country.get("borders", []))  # Convert to list instead of set
        
        
        if country_code in country_shapes:
            country_neighbors[country_code] = neighbors

    return country_neighbors

# Call the function
country_neighbors = get_country_neighbors()




def getPopulation(country_code):
    matching_country = filtered_world_data.loc[filtered_world_data['adm0_a3'] == country_code]
    if not matching_country.empty:
        return int(matching_country['pop_est'].values[0])
    else:
        return None

def getCountryBox(name): #Iterates through countries' polygon coordinates, finds lowest and highest x and y respectively

    name_polygons = country_shapes[name]
    leftTop = [float('inf'), float('inf')]
    rightBot = [0, 0]

    left=leftTop[0]
    top=leftTop[1]
    right=rightBot[0]
    bot=rightBot[1]

    for polygon in name_polygons:
        for x, y in polygon:
            if x < leftTop[0]:
                leftTop[0] = x
            if y < leftTop[1]:
                leftTop[1] = y
            if x > rightBot[0]:
                rightBot[0] = x
            if y > rightBot[1]:
                rightBot[1] = y

    return (leftTop,rightBot)


def getSubregionBoxDict():
    subregion_boxes_dict={}
    subregions_list = filtered_world_data['subregion'].unique()

    for subregion in subregions_list:
        leftTop = [float('inf'), float('inf')]
        rightBot = [0, 0]

        for _, row in filtered_world_data[filtered_world_data['subregion'] == subregion].iterrows():
            country_code = row['adm0_a3']
            if country_code in country_shapes.keys():

                name_polygons = country_shapes[country_code]

                for polygon in name_polygons:
                    for x, y in polygon:

                        if x < leftTop[0]:
                            leftTop[0] = x
                        if y < leftTop[1]:
                            leftTop[1] = y

                        if x > rightBot[0]:
                            rightBot[0] = x
                        if y > rightBot[1]:
                            rightBot[1] = y

        subregion_boxes_dict[subregion]=(leftTop,rightBot)
    return subregion_boxes_dict

subregion_boxes_dict=getSubregionBoxDict()


def getSubregionCountries():
    sub_countries_dict = {}
    for subregion in subregion_boxes_dict:
        subregion_countries_dict = {}
        countries_in_subregion = filtered_world_data[filtered_world_data['subregion'] == f'{subregion}']['adm0_a3'].tolist()
        sub_countries_dict.update({str(subregion):countries_in_subregion})


    return sub_countries_dict
sub_countries_dict=getSubregionCountries()


def drawCountry(name):
    name_polygons = country_shapes[name]
    for polygon in name_polygons:
        L=[]
        for x, y in polygon:
            L += [x] + [y]

    drawPolygon(*L,fill='lightGray', border='Black', borderWidth=1,
             opacity=100, rotateAngle=0, dashes=False, visible=True)

def find_nearest_country(mouseX, mouseY, country_shapes, app):
    click_point = Point(mouseX, mouseY)

    nearest_country = None
    min_distance = 0

    if not app.countriesIn:
        return None

    for country_code in app.countriesIn:
        if country_code in country_shapes:
            for screen_polygon in country_shapes[country_code]:
                shapely_polygon = Polygon(screen_polygon)

                if shapely_polygon.contains(click_point):
                    return country_code

                distance = shapely_polygon.distance(click_point)
                if distance < min_distance:
                    min_distance = distance
                    nearest_country = country_code

    return nearest_country

def pointInLand(mouseX, mouseY, country_shapes, app):
    click_point = Point(mouseX, mouseY)

    nearest_country = None
    min_distance = 0

    if not app.countriesIn:
        return None

    for country_code in app.countriesIn:
        if country_code in country_shapes:
            for screen_polygon in country_shapes[country_code]:
                shapely_polygon = Polygon(screen_polygon)

                if shapely_polygon.contains(click_point):
                    return True
                
    return False

def calculate_distance(lon1, lat1, lon2, lat2):
    point1 = (lat1, lon1)
    point2 = (lat2, lon2)
    return geodesic(point1, point2).km

###########################################################################################
#Helper Functions for MVC
def withinSubregion(app,mouseX,mouseY):
    app.subregionsIn=[]
    for subregion, (leftTop, rightBot) in subregion_boxes_dict.items():
            leftX, leftY = leftTop
            rightX, rightY = rightBot

            if leftX <= mouseX <= rightX and leftY <= mouseY <= rightY:
                app.subregionsIn.append(subregion)

    if len(app.subregionsIn) == 0:
        app.subregionsIn = []

def withinCountryinSub(app, mouseX, mouseY):
    app.countriesIn = []

    if len(app.subregionsIn) == 0:
        return app.countriesIn

    SubDict = getSubregionCountries()
    for subregion in app.subregionsIn:

        for country_code in SubDict[subregion]:
            if country_code in country_shapes.keys():

                leftTop, rightBot = getCountryBox(country_code)

                leftX, leftY = leftTop
                rightX, rightY = rightBot

                if leftX <= mouseX <= rightX and leftY <= mouseY <= rightY:
                    app.countriesIn.append(country_code)


def getIncomeGroup(country_code):
    matching_country = filtered_world_data.loc[filtered_world_data['adm0_a3'] == country_code]
    if not matching_country.empty:
        return int(matching_country['income_grp'].values[0][0])
    else:
        return None
    

def calculate_infectability(income_grp):

    infectability_map = {1: 0.01, 2: 0.02, 3: 0.03, 4: 0.04, 5: 0.05}
    return infectability_map.get(income_grp, 0.05)



# INFECTION FUNCTIONS
country_code_to_name = {}
country_code_to_name[None]=None
for _, row in filtered_world_data.iterrows():
    country_code = row['adm0_a3']
    country_name = row['name']

    country_code_to_name[country_code]=country_name

country_data = {}
for country_code in country_shapes.keys():
    income_grp=getIncomeGroup(country_code)
    country_data[country_code] = {
        "status": "healthy",
        "infected_population": 0,
        "total_population": getPopulation(country_code),
        "income_grp": income_grp,
        "infectability": calculate_infectability(income_grp),
        "name": country_code_to_name[country_code]
    }





def updateInfection(app):
    new_infections = []

    for country_code, data in country_data.items():
        if data['status']=="fully_infected" or data["status"] == "infected":
            # If the country is infected, advance infection
            if data["status"] == "infected":
                spread_in_country(app, country_code)

            if data['status']=="fully_infected" or data["status"] == "infected":
            # Attempt to infect neighboring countries
                for neighbor in get_neighbors(country_code):
                    if neighbor in country_shapes:
                        if country_data[neighbor]["status"] == "healthy":
                            if attempt_infection(neighbor):
                                new_infections.append(neighbor)

    # Update the new infections
    for country_code in new_infections:
        country_data[country_code]["status"] = "infected"
        country_data[country_code]["infected_population"] = 100  # Starting infection number


def spread_in_country(app, country_code):
    
    spreadData = country_data[country_code]

    growth_rate = spreadData["income_grp"] 

    if spreadData["infected_population"] < spreadData["total_population"]:
        spreadData["infected_population"] += int(spreadData["infected_population"] * growth_rate)

        # Transition to fully infected if >95% infected
        
        if spreadData["infected_population"] >= spreadData["total_population"] * 0.95:
            spreadData["status"] = "fully_infected"

def attempt_infection(neighbor):

    infection_chance = 0.2  # 20% chance for infection spread to neighbors
    return random.random() < infection_chance

def get_neighbors(country_code):
    return country_neighbors.get(country_code, [])

###########################################################################################
#Actual App Functions for MVC

def onAppStart(app):
    app.width=1200
    app.height=800
    app.UIy=550
    app.nearest_country=None
    app.population=None
    app.subregionsIn=[]
    app.countriesIn = []
    # app.coordinates = get_coordinate(image_path)
    app.coordinates = (132.31951944444444, 34.29629722222222)
    app.click_coordinates = None
    app.geo_click_coordinates= None
    app.distance = None
    app.target_dot = geo_to_screen(app.coordinates[0], app.coordinates[1], app.width, app.height)
    app.timer_delay = 100
    app.stepsPerSecond = 4

    
    


def drawGuessDot(app):
    if app.click_coordinates:
        clickX, clickY = app.click_coordinates
        drawCircle(clickX, clickY, 6, fill='red', border='black')
        
        targetX,targetY=app.target_dot
        drawCircle(targetX, targetY, 6, fill='blue', border='black')

        dx = targetX - clickX
        dy = targetY - clickY

        if abs(dx) > app.width / 2:

            if clickX > targetX:

                slope = dy / dx if dx != 0 else float('inf')

                x1_end = app.width
                x2_start = 0
                
                y1_end = clickY - slope * (x1_end - clickX)
                y2_start = targetY + slope * (targetX - x2_start)

                y2_start=(y1_end+y2_start)/2
                y1_end=(y1_end+y2_start)/2

                drawLine(clickX, clickY, x1_end, y1_end, fill='black')
                drawLine(x2_start, y2_start, targetX, targetY, fill='black')
            else:

                slope = dy / dx if dx != 0 else float('inf')

                x1_end = 0
                x2_start = app.width
                
                y1_end = clickY - slope * (x1_end - clickX)
                y2_start = targetY + slope * (targetX - x2_start)
                
                y2_start=(y1_end+y2_start)/2
                y1_end=(y1_end+y2_start)/2

                drawLine(clickX, clickY, x1_end, y1_end, fill='black')
                drawLine(x2_start, y2_start, targetX, targetY, fill='black')
        else:
            # Draw a direct line if no wrapping is needed
            drawLine(clickX, clickY, targetX, targetY, fill='black')



def drawCountries(app):
    for country_code in list(country_shapes.keys()):
        spreadData = country_data[country_code]

        name_polygons = country_shapes[country_code]
        for polygon in name_polygons:
            L = []
            for x, y in polygon:
                L += [x] + [y]
            if country_code == app.nearest_country:
                color = 'gray'
            if spreadData["status"] == "infected":
                color = "orange"
            elif spreadData["status"] == "fully_infected":
                color = "darkRed"
            elif spreadData["status"] == "eradicated":
                color = "green"
            else:
                color = 'lightGray'
            drawPolygon(*L, fill=color, border='Black', borderWidth=1,
                        opacity=100, rotateAngle=0, dashes=False, visible=True)
            

def onStep(app):
    updateInfection(app)


def redrawAll(app):
    drawCountries(app)
    drawGuessDot(app)

    drawRect(0,app.UIy,app.width,app.height-app.UIy,fill='linen')
    drawLabel(f"Country: {country_code_to_name[app.nearest_country]}",800,600,size=25)
    drawLabel(f"Population: {app.population}",800,625,size=25)
    # drawLabel(f"Subregion(s): {app.subregionsIn}",800,650,size=25)
    drawLabel(f"In Countries: {app.countriesIn}",800,675,size=25)
    drawLabel(f"Coordinates: {app.geo_click_coordinates}",300,600,size=25)
    drawLabel(f"Actual Coordinates: {app.coordinates}",300,625,size=25)
    drawLabel(f"distance: {app.distance}km",300,650,size=25)

def start_infection(country_code):
    if country_code in country_data:
        # Check if the country is currently healthy before starting an infection
        if country_data[country_code]["status"] == "healthy":
            country_data[country_code]["status"] = "infected"
            country_data[country_code]["infected_population"] = 10  # Initial infected population


def onMousePress(app, mouseX, mouseY):
    if mouseY<app.UIy:
        withinSubregion(app,mouseX,mouseY)
        withinCountryinSub(app,mouseX,mouseY)
        app.nearest_country = find_nearest_country(mouseX, mouseY, country_shapes, app)
        app.population = getPopulation(app.nearest_country)


        app.click_coordinates = (mouseX, mouseY)
        click_lon, click_lat = screen_to_geo(mouseX, mouseY, app.width, app.height)
        app.distance = int(calculate_distance(app.coordinates[0], app.coordinates[1], click_lon, click_lat))
        app.geo_click_coordinates = (int(click_lon), int(click_lat))
        if app.nearest_country:
            start_infection(app.nearest_country)



app.setMaxShapeCount(4000)

runApp()



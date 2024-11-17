#pip install geopandas shapely pandas numpy geopy requests

#Imports
from cmu_graphics import *
import geopandas as gpd
from shapely.geometry import MultiPolygon
from shapely.geometry import Point
from shapely.geometry import Polygon
import pandas as pd
import json
import requests
import random
import pyproj
###########################################################################################
# Json File Reading
with open('.src/geojson-maps.json', 'r') as f:
    json_data = json.load(f)

data = []
for feature in json_data["features"]:
    properties = feature["properties"]

    if properties.get("adm0_a3") == "SDS":  # I HATE SOUTH SUDAN
        properties["adm0_a3"] = "SSD"

    if "adm0_a3" in properties and "name" in properties and "pop_est" in properties and "subregion" in properties:
        data.append({
            "adm0_a3": properties["adm0_a3"],
            "name": properties["name"],
            "pop_est": properties["pop_est"],
            "subregion": properties["subregion"]
        })

df = pd.DataFrame(data)
exceptions = ["GRL", "ISL"]

filtered_df = df[(df["pop_est"] >= 2000000) & (df["adm0_a3"] != "SLE") & (df["adm0_a3"] != "BDI") & (df["adm0_a3"] != "LSO")]
filtered_df = filtered_df[filtered_df["subregion"].str.contains("Africa")]


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
# Start of ChatGPT generated/supported segment ######################################################
def geo_to_screen(lon, lat, width, height):
    # Africa bounding box (approximately)
    min_lon, max_lon = -20, 55    # Longitude range of Africa
    min_lat, max_lat = -35, 37    # Latitude range of Africa

    # Convert lat/lon to screen space using Mercator projection
    project = pyproj.Transformer.from_proj(
        proj_from=pyproj.Proj("epsg:4326"),  # WGS84 coordinate system
        proj_to=pyproj.Proj("epsg:3857")     # Web Mercator projection
    ).transform

    # Project the coordinates
    min_x, min_y = project(min_lat, min_lon)
    max_x, max_y = project(max_lat, max_lon)
    mercator_x, mercator_y = project(lat, lon)

    # Normalize based on Africa's coordinates and available screen space
    scale_x = width / (max_x - min_x)
    scale_y = height / (max_y - min_y)

    # Preserve aspect ratio to avoid distortion
    scale = min(scale_x, scale_y)

    screen_x = int((mercator_x - min_x) * scale)
    screen_y = int((max_y - mercator_y) * scale)
    return screen_x, screen_y

def screen_to_geo(screen_x, screen_y, width, height):
    # Convert screen space back to lat/lon using Mercator projection
    inverse_project = pyproj.Transformer.from_proj(
        proj_from=pyproj.Proj("epsg:3857"),  # Web Mercator projection
        proj_to=pyproj.Proj("epsg:4326")     # WGS84 coordinate system
    ).transform

    mercator_x = screen_x * (40075016.68 / width) - 20037508.34
    mercator_y = 20037508.34 - screen_y * (40075016.68 / height)
    lon, lat = inverse_project(mercator_y, mercator_x)
    return lon, lat


# End of ChatGPT generated/supported segment ########################################################
###########################################################################################

#Generate Shapes and Helper Functions
country_shapes = {}
for _, row in filtered_world_data.iterrows():
    country_name = row['name']
    geom = row['geometry']
    if isinstance(geom, MultiPolygon):
        polygons = list(geom.geoms)
    else:
        polygons = [geom]

    screen_polygons = []
    for polygon in polygons:
        simplified_polygon = polygon.simplify(0.1) ###IMPORTANT FOR FASTER LOADING SPEED: simplify the higher float the more simple
        screen_coords = [geo_to_screen(x, y, 1200, 550) for x, y in simplified_polygon.exterior.coords]        
        screen_polygons.append(screen_coords)
        country_shapes[country_name] = screen_polygons

def getPopulation(country_name):
    matching_country = filtered_world_data.loc[filtered_world_data['name'] == country_name]
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
            country_name = row['name']
            if country_name in country_shapes.keys():

                name_polygons = country_shapes[country_name]

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
        countries_in_subregion = filtered_world_data[filtered_world_data['subregion'] == f'{subregion}']['name'].tolist()
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

def find_nearest_country(mouse_x, mouse_y, country_shapes, app):
    click_point = Point(mouse_x, mouse_y)

    nearest_country = None
    min_distance = 0

    if not app.countriesIn:
        return None

    for country_name in app.countriesIn:
        if country_name in country_shapes:
            for screen_polygon in country_shapes[country_name]:
                shapely_polygon = Polygon(screen_polygon)

                if shapely_polygon.contains(click_point):
                    return country_name

                distance = shapely_polygon.distance(click_point)
                if distance < min_distance:
                    min_distance = distance
                    nearest_country = country_name

    return nearest_country

country_codes = set(filtered_world_data['adm0_a3'])
def get_country_neighbors():
    response = requests.get("https://restcountries.com/v3.1/all")
    response.raise_for_status()
    
    countries_data = response.json()
    
    country_neighbors = {}
    
    for country in countries_data:
                
        country_code = country.get("cca3")
            
        neighbors = list(country.get("borders", []))  # Convert to list instead of set
        
        
        if country_code in country_codes:
            country_neighbors[country_code] = neighbors

    return country_neighbors



country_neighbors = get_country_neighbors()
country_neighbors['MOZ']=country_neighbors.get('MOZ')+['MDG']
country_neighbors['MDG']=country_neighbors.get('MDG')+['MOZ']

print(country_neighbors.keys())

country_code_to_name = {}
country_code_to_name[None]=None

country_name_to_code = {}
country_name_to_code[None]=None

for _, row in filtered_world_data.iterrows():
    country_code = row['adm0_a3']
    country_name = row['name']

    country_code_to_name[country_code]=country_name
    country_name_to_code[country_name]=country_code


territories={1: ['EGY', 'LBY', 'TUN', 'DZA', 'MAR', 'SDN'],  # North Africa
    2: ['CIV', 'BEN', 'TGO', 'GHA', 'SEN', 'NGA', 'GMB', 'MLI', 'BFA', 'NER', 'GNB', 'GHA', 'MRT'],  # West Africa
    3: ['CAM', 'GAB', 'CAF', 'COD','GAB', 'COG','CMR'],  # Central Africa
    4: ['ETH', 'KEN', 'TZA', 'SOM', 'UGA', 'RWA', 'ERI', 'SSD'],  # East Africa
    5: ['ZAF', 'BWA', 'NAM', 'ZWE', 'AGO'],  # Southern Africa
    6: ['MDG'],  # Island Countries
    7: ['MOZ', 'LBR', 'ZMB','MWI'],  # Central and Southern Regions
             }

region_colors = {
        1: "lightBlue",  # North Africa
        2: "lightGreen",  # West Africa
        3: "lightYellow",  # Central Africa
        4: "lightCoral",   # East Africa
        5: "lightPink",    # Southern Africa
        6: "black",  # Horn of Africa
        7: "red",  # Great Lakes Region
    }

def get_neighbors(country_code):
    return country_neighbors.get(country_code, [])

def get_center(name_polygons):
    x_coords = []
    y_coords = []
    
    for polygon in name_polygons:
        for x, y in polygon:
            x_coords.append(x)
            y_coords.append(y)
    
    mean_x = sum(x_coords) / len(x_coords) if x_coords else None
    mean_y = sum(y_coords) / len(y_coords) if y_coords else None
    
    return (mean_x, mean_y)


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

        for country_name in SubDict[subregion]:
            if country_name in country_shapes.keys():

                leftTop, rightBot = getCountryBox(country_name)

                leftX, leftY = leftTop
                rightX, rightY = rightBot

                if leftX <= mouseX <= rightX and leftY <= mouseY <= rightY:
                    app.countriesIn.append(country_name)


def get_random_half_countries(country_shapes):
    country_list = list(country_shapes.keys())
    random.shuffle(country_list) 
    half_count = len(country_list) // 2  
    return country_list[:half_count]


class Player:

    def __init__(self,startingCountries,color):
        self.active=False
        self.color=color
        self.owned=startingCountries
        self.phases=['Reinforcement','Attack','Fortification']
        self.phaseIndex=0
        self.gamePhase=[]
    
    def drawArmies(self):
        for country in self.owned:
            pass
    
    def fortify(self):
        pass

    
        

class Game:
    def __init__(self):
        pass

    def start(self):
        starting1 = set(get_random_half_countries(country_shapes))
        starting2 = set(country_shapes) - starting1
        player1 = Player(starting1, 'green')
        player2 = Player(starting2, 'blue')
        self.players = [player1, player2]
    

###########################################################################################
#Actual App Functions for MVC

def onAppStart(app):
    app.width=1200
    app.height=800
    app.UIy=550
    app.nearest_country='Congo'
    app.population=None
    app.subregionsIn=[]
    app.countriesIn = []
    app.neighbors= []
    app.tView=False
    starting1=set(get_random_half_countries(country_shapes))
    starting2=set(country_shapes.keys())-starting1
    player1=Player(starting1,'green')
    player2=Player(starting2,'blue')

    app.players=[player1,player2]
    app.activePlayer=app.players[0]



def drawCountries(app):
    circleCenter=[]
    for country_name in list(country_shapes.keys()):
        name_polygons = country_shapes[country_name]
        

        for polygon in name_polygons:
            L=[]
            for x, y in polygon:
                L += [x] + [y]


            if app.tView:

                region = None
                for region_id, countries in territories.items():
                    if country_name_to_code[country_name] in countries:
                        region = region_id
                        break
                
                # If a region is found, color the country accordingly
                if region is not None:
                    color = region_colors.get(region, "lightGray")
                else:
                    color = "lightGray"  # Default color if no region is found
                        


            else:       
                if country_name==app.nearest_country:
                    color='dimGray'
                elif country_name_to_code[country_name] in app.neighbors:
                    color='red'
                else:
                    color='lightGray'
            
            drawPolygon(*L,fill=color, border='Black', borderWidth=1,
                opacity=100, rotateAngle=0, dashes=False, visible=True)
        
        x,y=get_center(name_polygons)
        circleCenter.append((x,y))
    
    if not app.tView:
        for x,y in circleCenter:
            print(country_name in app.players[0].owned)
            if country_name in app.players[0].owned:
                drawCircle(x,y,10,fill="green")
                drawLabel("1",x,y,size=18,bold=True)
            else:
                drawCircle(x,y,10,fill="red")
                drawLabel("1",x,y,size=18,bold=True)

def onKeyPress(app,key):
    if key=='t':
        app.tView=not app.tView

def redrawAll(app):
    drawCountries(app)


    drawRect(0,app.UIy,app.width,app.height-app.UIy,fill='linen')
    drawLabel(f"Country: {country_name_to_code[app.nearest_country]}",650,600,size=25)
    drawLabel(f"Population: {app.population}",650,625,size=25)
    drawLabel(f"Neighbor(s): {app.neighbors}",650,650,size=25)
    drawLabel(f"In Countries: {app.countriesIn}",650,675,size=25)

def onMouseMove(app, mouseX, mouseY):
    if mouseY<app.UIy:
        withinSubregion(app,mouseX,mouseY)
        withinCountryinSub(app,mouseX,mouseY)
        app.nearest_country = find_nearest_country(mouseX, mouseY, country_shapes, app)
        app.population = getPopulation(app.nearest_country)
        app.neighbors=get_neighbors(country_name_to_code[app.nearest_country])

app.setMaxShapeCount(4000)

runApp()



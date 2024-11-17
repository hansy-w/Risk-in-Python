#pip install geopandas shapely pandas numpy geopy requests

#Imports
from cmu_graphics import *
import geopandas as gpd
from shapely.geometry import MultiPolygon, GeometryCollection
from shapely.ops import transform
from shapely.geometry import Point
from shapely.geometry import Polygon
import pandas as pd
import json
import numpy as np
###########################################################################################
# Json File Reading
with open('./src/geojson-maps.json', 'r') as f:
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

filtered_df = df[(df["pop_est"] >= 2000000) | (df["adm0_a3"].isin(exceptions))]

# GeoPandas file reading
geojson_path = './src/datahub.geojson'
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


# End of ChatGPT generated/supported segment ##############################################
###########################################################################################

territories = {
    1: ['United States of America'],
    2: ['Mexico', 'Guatemala', 'Honduras', 'El Salvador', 'Nicaragua', 'Costa Rica', 'Panama','Jamaica', 'Haiti', 'Dominican Rep.', 'Cuba', 'Puerto Rico'],
    3: ['Colombia', 'Venezuela', 'Ecuador'],
    4: ['Brazil'],
    5: ['Argentina', 'Chile', 'Uruguay', 'Paraguay'],
    6: ['Peru', 'Bolivia'],
    7: ['United Kingdom', 'Ireland', 'Iceland'],
    8: ['France', 'Belgium', 'Netherlands', 'Switzerland'],
    9: ['Germany', 'Austria', 'Czechia'],
    10: ['Spain', 'Portugal'],
    11: ['Italy', 'Greece', 'Albania', 'Croatia', 'Bosnia and Herz.'],
    12: ['Norway', 'Sweden', 'Denmark'],
    13: ['Finland', 'Lithuania', 'Latvia', 'Estonia'],
    14: ['Poland', 'Hungary', 'Slovakia', 'Slovenia'],
    15: ['Romania', 'Bulgaria', 'Serbia', 'North Macedonia'],
    16: ['Ukraine', 'Moldova', 'Belarus'],
    17: ['Russia'],
    18: ['Kazakhstan', 'Uzbekistan', 'Turkmenistan', 'Kyrgyzstan', 'Tajikistan'],
    19: ['Azerbaijan', 'Armenia', 'Georgia','Turkey', 'Cyprus'],
    20: ['Iran', 'Iraq', 'Syria', 'Lebanon', 'Jordan', 'Israel'],
    21: ['Saudi Arabia', 'Yemen', 'Oman', 'United Arab Emirates', 'Qatar', 'Kuwait'],
    22: ['Egypt', 'Libya', 'Sudan', 'S. Sudan'],
    23: ['Morocco', 'Algeria', 'Tunisia', 'Mauritania'],
    24: ['Nigeria', 'Niger', 'Chad'],
    25: ['Cameroon', 'Central African Rep.', 'Congo', 'Dem. Rep. Congo', 'Gabon'],
    26: ['Ethiopia', 'Somalia', 'Eritrea'],
    27: ['Kenya', 'Uganda', 'Tanzania', 'Rwanda', 'Burundi'],
    28: ['Angola', 'Zambia', 'Malawi'],
    29: ['Mozambique', 'Zimbabwe', 'Botswana'],
    30: ['South Africa', 'Namibia', 'Lesotho', 'Eswatini'],
    31: ['India', 'Sri Lanka','Bangladesh', 'Nepal', 'Bhutan'],
    32: ['Pakistan', 'Afghanistan'],
    33: ['China', 'Taiwan'],
    34: ['Japan', 'South Korea', 'North Korea'],
    35: ['Mongolia'],
    36: ['Vietnam', 'Laos', 'Cambodia'],
    37: ['Thailand', 'Myanmar'],
    38: ['Philippines', 'Malaysia', 'Singapore','Indonesia', 'Papua New Guinea'],
    39: ['Australia', 'New Zealand'],
    40: ['Mali', 'Burkina Faso', "Côte d'Ivoire", 'Ghana', 'Benin', 'Togo'],
    41: ['Senegal', 'Gambia', 'Guinea', 'Sierra Leone', 'Liberia'],
    42: ['Madagascar'],
    43: ['Canada'],
    44:['Greenland']
}

def map_countries_to_territories(territories):
    country_to_territory = {}
    for territory_num, countries in territories.items():
        for country in countries:
            country_to_territory[country] = territory_num

    return country_to_territory
country_to_territory = map_countries_to_territories(territories)

region_colors = {
    "Europe": "blue",
    "North America": "green",
    "Central America": "yellow",
    "South America": "red",
    "Africa": "orange",
    "Oceania": "purple",
    "Asia": "gray"
}

# Map territories to regions
territory_regions = {
    "Europe": [7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19],
    "North America": [1, 43, 44],
    "Central America": [2],
    "South America": [3, 4, 5, 6],
    "Africa": [22, 23, 24, 25, 26, 27, 28, 29, 30,40,41],
    "Oceania": [38, 39, 42],
    "Asia":[20, 21, 31, 32, 33, 34, 35, 36, 37]
}

colors_to_rgb={
    "blue":(65, 105, 225),
    "green":(34, 139, 34),
    "yellow":(255, 228, 181),
    "red":(178, 34, 34),
    "orange":(255, 99, 71),
    "purple":(138, 43, 226),
    "gray":(47, 79, 79),
    "lightGray":(211, 211, 211)
}

# Helper function to adjust brightness for gradients
def adjust_brightness(color, factor):
    r, g, b = colors_to_rgb[color]
    r = min(255, int(r * factor))
    g = min(255, int(g * factor))
    b = min(255, int(b * factor))
    return rgb(r,g,b)

#Generate Shapes and Helper Functions
country_shapes = {}
for _, row in filtered_world_data.iterrows():
    country_name = row['name']
    geom = row['geometry']
    
    if isinstance(geom, MultiPolygon):
        polygons = list(geom.geoms)
    elif isinstance(geom, Polygon):
        polygons = [geom]
    elif isinstance(geom, GeometryCollection):
        # If it's a GeometryCollection, collect all polygons within it
        polygons = [g for g in geom.geoms if isinstance(g, (Polygon, MultiPolygon))]
    else:
        continue  # Skip if not a Polygon or MultiPolygon

    screen_polygons = []
    for polygon in polygons:
        simplified_polygon = polygon.simplify(0.3)  # Simplify for faster loading
        screen_coords = [geo_to_screen(x, y, 1200, 800) for x, y in simplified_polygon.exterior.coords]
        screen_polygons.append(screen_coords)
        
    country_shapes[country_name] = screen_polygons

def getPopulation(country_name):
    matching_country = filtered_df.loc[filtered_df['name'] == country_name]
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
    subregions_list = filtered_df['subregion'].unique()

    for subregion in subregions_list:
        leftTop = [float('inf'), float('inf')]
        rightBot = [0, 0]

        for _, row in filtered_df[filtered_df['subregion'] == subregion].iterrows():
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
        countries_in_subregion = filtered_df[filtered_df['subregion'] == f'{subregion}']['name'].tolist()
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



###########################################################################################
#Actual App Functions for MVC

def onAppStart(app):
    app.width = 1200
    app.height = 800
    app.UIy = 550
    app.nearest_country = None
    app.population = None
    app.subregionsIn = []
    app.countriesIn = []
    app.selected_territory = None




def drawCountries(app):
    for territory_id, country_list in territories.items():
        # get region color
        region = None
        for reg_name, reg_territories in territory_regions.items():
            if territory_id in reg_territories:
                region = reg_name
                break

        base_color = region_colors.get(region, "lightGray")
        
        # adjust brightness 
        gradient_factor_step = 0.8 / len(country_list)
        for i, country_name in enumerate(country_list):
            if country_name not in country_shapes:
                continue
            
            # Calculate  color 
            gradient_factor = 1.2 - (i * gradient_factor_step)
            country_color = adjust_brightness(base_color, gradient_factor)

            for polygon in country_shapes[country_name]:
                L = []
                for x, y in polygon:
                    L += [x, y]
                drawPolygon(*L, fill=country_color, border="black", borderWidth=1,
                            opacity=100, rotateAngle=0, dashes=False, visible=True)



def redrawAll(app):
    drawCountries(app)


    drawRect(0,app.UIy,app.width,app.height-app.UIy,fill='linen')
    drawLabel(f"Country: {app.nearest_country}",650,600,size=25)
    drawLabel(f"Population: {app.population}",650,625,size=25)
    # drawLabel(f"Subregion(s): {app.subregionsIn}",650,650,size=25)
    # drawLabel(f"In Countries: {app.countriesIn}",650,675,size=25)

def onMouseMove(app, mouseX, mouseY):
    if mouseY < app.UIy:
        withinSubregion(app, mouseX, mouseY)
        withinCountryinSub(app, mouseX, mouseY)
        app.nearest_country = find_nearest_country(mouseX, mouseY, country_shapes, app)

        app.population = getPopulation(app.nearest_country)

        if app.nearest_country in country_to_territory:
            app.selected_territory = country_to_territory[app.nearest_country]  
            app.countriesIn = territories[app.selected_territory]  
        else:
            app.selected_territory = None  # Deselect if no country is clicked

211, 211, 211
app.setMaxShapeCount(4000)

runApp()



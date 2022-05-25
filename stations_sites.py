import geopandas as gpd
import pandas as pd
import json
import numpy as np
from typing import List, Tuple
import folium
# from screened_sites import low_coverage_sites


def add_geometry(threshold:float=0) -> None:
    """
    Geocode json list of State's properties and output the .geojson file. 
    Use a threshold in [0,1] to be applied on geocoding score
    The State's properties json file is open data from data.gouv.fr
    """
    coordinates_dict = {}
    with open("inventaire-immobilier-de-letat.json") as f:
        buildings = json.load(f)

    buildings_df = pd.DataFrame([building["fields"] for building in buildings])
    print(f"Number of buildings {len(buildings_df)}")

    for dep in np.concatenate(
        (
            np.arange(1, 20),
            np.arange(21, 96),
            [97, 99],
            np.arange(971, 979),
            [986, 987, 988, 989],
        )
    ):
        coordinates_dict[dep].rename(
            columns={"latitude": "lat", "longitude": "lon"}, inplace=True
        )
    coordinates_df = pd.concat(coordinates_dict.values(), axis=0)
    print("Coordinates")
    geocoded_df = pd.merge(
        buildings_df, coordinates_df, left_index=True, right_index=True
    )
    geocoded_df = geocoded_df[geocoded_df["result_score"] >= threshold]
    geocoded_df.drop(columns=["result_score"]).to_csv("geocoded.csv")

    geobuildings = gpd.GeoDataFrame(
        geocoded_df, geometry=gpd.points_from_xy(geocoded_df.lon, geocoded_df.lat)
    )

    print("GEOBUILDINGS")
    print(geobuildings.head())
    geobuildings.to_file("inventaire-immobilier-de-letat_geocoded.geojson")
    return


def load_stations_coordinates() -> List[List[float]]:
    """
    Stations geojson is open data from : https://www.infoclimat.fr/opendata/stations_xhr.php?format=geojson

    Returns:
        List[List[float]]: _description_
    """
    stations = gpd.read_file("stations.geojson")
    return [[p.y, p.x] for p in stations.geometry]


def load_candidate_site_coordinates() -> List[List[float]]:
    """
    The files here is basically a geocoding of the dictionary of low_coverage_sites that we have previously screened
    """
    coords = []
    for _,p in gpd.read_file("ALREADY_SAVED_INTERESTING_SITES.geojson").iterrows():
        if not pd.isna(p.lat) and not pd.isna(p.lat):
            coords.append([p.lat, p.lon])
    return coords


def load_candidate_addresses() -> List[str]:
    """
    The files here is basically a geocoding of the dictionary of low_coverage_sites that we have previously screened
    Aligned with ``load_candidate_site_coordinates``

    Returns:
        Lits[str]: list of addresses text fields combining city names, street names and building function
    """
    addresses = []
    for _,p in gpd.read_file("ALREADY_SAVED_INTERESTING_SITES.geojson").iterrows():
        content = ""
        if not pd.isna(p.lat) and not pd.isna(p.lat):
            if p.ville:
                content += p.ville + "_"
            if p.adresse:
                content += p.adresse + "_"
            if p.fonction:
                content += p.fonction
            addresses.append(content)
    return addresses


def save_map() -> None:
    """ 
    Save a .html where ticks with stations are put on open street view
    """
    site_coords = load_stations_coordinates()
    cand_coords = load_candidate_site_coordinates()
    cand_towns = load_candidate_addresses()
    m = folium.Map(location=[48.86, 2.33])
    tooltip = ""

    for coord in site_coords:
        folium.Marker(
            coord,
            popup="",
            icon=folium.Icon(color="red"),
            tooltip=tooltip,
        ).add_to(m)
    for k, coord in enumerate(cand_coords):
        folium.Marker(
            coord,
            popup=f"<i>{cand_towns[k]}</i>",
            icon=folium.Icon(color="orange"),
            tooltip=tooltip,
        ).add_to(m)
    m.save("index.html")

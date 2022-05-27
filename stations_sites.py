import geopandas as gpd
import pandas as pd
import json
import numpy as np
from typing import List, Tuple
import folium
from screened_sites import low_coverage_sites

DEPARTEMENTS = np.concatenate(
    (
        np.arange(1, 20),
        np.arange(21, 96),
        [97, 99],
        np.arange(971, 979),
        [986, 987, 988, 989],
        ["2A", "2B"],
    )
)

def add_geometry(threshold: float = 0) -> None:
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

    for dep in DEPARTEMENTS:
        # the departement buildings must have been previously geocoded with API ADRESSE : https://adresse.data.gouv.fr/api-doc/adresse
        coordinates_dict[dep] = pd.read_csv(f"GEOCODED_{dep}.csv", index_col=0)[
            ["latitude", "longitude", "result_score"]
        ]
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


def load_stations_coordinates(departement: str = "ALL") -> List[List[float]]:
    """
    Stations geojson is open data from : https://www.infoclimat.fr/opendata/stations_xhr.php?format=geojson

    Returns:
        List[List[float]]: list of [latitude, longitude] of the stations
    """
    stations = gpd.read_file("stations.geojson")
    if departement in [str(k) for k in range(1, 10)]:
        departement = "0" + departement
    if departement.upper() != "ALL":
        stations = stations[stations.departement == str(departement).upper()]
    return [[p.y, p.x] for p in stations.geometry]


def load_candidate_site_coordinates(departement: str = "ALL") -> List[List[float]]:
    """
    The files here is basically a geocoding of the dictionary of low_coverage_sites that we have previously screened
    """
    coords = []
    sites = gpd.read_file("ALREADY_SAVED_INTERESTING_SITES.geojson")
  #  if departement in [str(k) for k in range(1, 10)]:
  #      departement = "0" + departement
    if departement.upper() != "ALL":
        sites = sites[sites.dept == str(departement).upper()]
    for _, p in sites.iterrows():
        if not pd.isna(p.lat) and not pd.isna(p.lat):
            coords.append([p.lat, p.lon])
    return coords


def load_candidate_addresses(departement: str = "ALL") -> List[str]:
    """
    The files here is basically a geocoding of the dictionary of low_coverage_sites that we have previously screened
    Aligned with ``load_candidate_site_coordinates``

    Returns:
        Lits[str]: list of addresses text fields combining city names, street names and building function
    """
    addresses = []
    sites = gpd.read_file("ALREADY_SAVED_INTERESTING_SITES.geojson")
   # if departement in [str(k) for k in range(1, 10)]:
   #     departement = "0" + departement
    if departement.upper() != "ALL":
        sites = sites[sites.dept == str(departement).upper()]
    for _, p in sites.iterrows():
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


def geocode_interesting_sites() -> None:
    """
    Geocode the list of interesting sites that have "manually" screened
    At this point there is no regularization of city names, so this is really a basic sandbox version
    """
    properties = gpd.read_file("inventaire-immobilier-de-letat_geocoded.geojson")

    interesting_sites = np.concatenate(list(low_coverage_sites.values()))
    properties["major-interest"] = properties["ville"].map(
        lambda ville: ville in interesting_sites
    )
    major_interest = properties[properties["major-interest"]]
    major_interest.to_file("ALREADY_SAVED_INTERESTING_SITES.geojson")


def save_map(departement: str = "ALL") -> None:
    """
    Save a .html where ticks with stations are put on open street view

    Args:
        departement (str): number of the departement, as a string. Examples: '06' (Nice), '2B' (Haute-Corse), '75' (Paris), or "ALL" for whole France
    """
    site_coords = load_stations_coordinates(departement)
    if not site_coords:
        return
    cand_coords = load_candidate_site_coordinates(departement)
    cand_towns = load_candidate_addresses(departement)
    if departement.upper() == "ALL":
        m = folium.Map(location=[48.86, 2.33])
    else:
        m = folium.Map(location=site_coords[0])
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
    m.save(f"index_{departement}.html")

if __name__ == "__main__":
    save_map("ALL")
  #  for dep in DEPARTEMENTS:
  #      save_map(str(dep))
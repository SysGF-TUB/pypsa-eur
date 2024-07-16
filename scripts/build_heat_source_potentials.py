import geopandas as gpd
from typing import List
from _helpers import set_scenario_config

"""
Maps heat source potential data clustered onshore regions. Heat source potentials are based on Manz et al. 2024: "Spatial analysis of renewable and excess heat potentials for climate-neutral district heating in Europe".
The `ClusteredRegionData` class encapsulates the process of mapping energy data to geographical regions represented by clustered regions and converting this data into a consistent energy unit.

Inputs:
- `data/clustered_regions.geojson`: GeoJSON file containing the geometries of clustered regions.
- `data/energy_data.geojson`: GeoJSON file containing energy data to be mapped to the clustered regions.

Outputs:
- `results/mapped_energy_data.geojson`: GeoJSON file containing the energy data mapped to the clustered regions, with energy data converted to the default unit if necessary.

Relevant Settings
-----------------

.. code:: yaml
    sector:
        district_heating:
            heat_sources:
                industrial_excess_heat: 
                    source_data:
                    data_column:
                    temperature_celsius: 
                    source_type: 
                    unit: 

References
----------
- `Manz et al. 2024: "Spatial analysis of renewable and excess heat potentials for climate-neutral district heating in Europe" <https://www.sciencedirect.com/science/article/pii/S0960148124001769>`.

"""


class ClusteredRegionData:
    """
    Manages and processes geographical data related to clustered regions, focusing on mapping and converting energy data associated with these regions.

    Parameters
    ----------
    clustered_regions : gpd.GeoDataFrame
        A GeoDataFrame containing the geometries of clustered regions.
    data : gpd.GeoDataFrame
        A GeoDataFrame containing the data to be mapped to the clustered regions. This data is automatically converted to the coordinate reference system (CRS) of the clustered_regions upon initialization.
    data_column : str
        The name of the column in 'data' GeoDataFrame that contains the energy data to be processed.
    energy_unit : str, optional
        The unit of the energy data provided in 'data'. Defaults to "MWh".
    default_energy_unit : str, optional
        The default unit of energy to which the data will be converted if necessary. Defaults to "MWh".
    mapped : bool, optional
        A flag indicating whether the data has been mapped to the clustered regions. Defaults to False.
    converted : bool, optional
        A flag indicating whether the energy data has been converted to the default energy unit. Defaults to False.

    Attributes
    ----------
    clustered_regions : gpd.GeoDataFrame
        Stores the geometries of clustered regions.
    data : gpd.GeoDataFrame
        Stores the data to be mapped to the clustered regions after conversion to the CRS of clustered_regions.
    data_column : str
        Indicates the column in 'data' containing the energy data.
    energy_unit : str
        Indicates the unit of the energy data in 'data'.
    default_energy_unit : str
        Indicates the default unit of energy for conversion.
    _mapped : bool
        Internal flag to indicate if data has been mapped.
    converted : bool
        Indicates if the energy data has been converted to the default unit.

    Methods
    -------
    data_in_regions
        Returns a GeoDataFrame containing the data mapped to the clustered regions. If the data has not been mapped yet, it performs the mapping. If the energy data has not been converted to the default energy unit, it performs the conversion.

    _map_to_clustered_regions
        (Private) Maps the energy data to the clustered regions based on their geometries. Called internally by the data_in_regions property.

    _convert_energy
        (Private) Converts the energy data to the default energy unit. Intended to be implemented to perform unit conversion and is called internally by the data_in_regions property.

    Usage
    -----
    Encapsulates the process of mapping energy data to geographical regions represented by clustered regions and converting this data into a consistent energy unit.
    """
    def __init__(
        self,
        clustered_regions: gpd.GeoDataFrame,
        data: gpd.GeoDataFrame,
        data_column: str,
        energy_unit: str = "MWh",
        default_energy_unit: str = "MWh",
        mapped: bool = False,
        converted: bool = False,
    ) -> None:
        """
        Initializes the ClusteredRegionData object with the necessary data and configurations.

        Parameters
        ----------
        clustered_regions : gpd.GeoDataFrame
            A GeoDataFrame containing the geometries of clustered regions.
        data : gpd.GeoDataFrame
            A GeoDataFrame containing the data to be mapped to the clustered regions. This data is automatically converted to the coordinate reference system (CRS) of the clustered_regions upon initialization.
        data_column : str
            The name of the column in 'data' GeoDataFrame that contains the energy data to be processed.
        energy_unit : str, optional
            The unit of the energy data provided in 'data'. Defaults to "MWh".
        default_energy_unit : str, optional
            The default unit of energy to which the data will be converted if necessary. Defaults to "MWh".
        mapped : bool, optional
            A flag indicating whether the data has been mapped to the clustered regions. Defaults to False.
        converted : bool, optional
            A flag indicating whether the energy data has been converted to the default energy unit. Defaults to False.

        Attributes
        ----------
        clustered_regions : gpd.GeoDataFrame
            Stores the geometries of clustered regions.
        data : gpd.GeoDataFrame
            Stores the data to be mapped to the clustered regions after conversion to the CRS of clustered_regions.
        data_column : str
            Indicates the column in 'data' containing the energy data.
        energy_unit : str
            Indicates the unit of the energy data in 'data'.
        default_energy_unit : str
            Indicates the default unit of energy for conversion.
        _mapped : bool
            Internal flag to indicate if data has been mapped.
        converted : bool
            Indicates if the energy data has been converted to the default unit.

        Returns
        -------
        None
        """

        self.clustered_regions = clustered_regions
        self.data = data.to_crs(clustered_regions.crs)
        self.data_column = data_column

        self._mapped = mapped
        if mapped:
            self._data_in_regions = data

        self.energy_unit = energy_unit
        self.default_energy_unit = default_energy_unit
        self.converted = False

    @property
    def data_in_regions(self) -> gpd.GeoDataFrame:
            """
        Retrieves or calculates the GeoDataFrame containing the data mapped to the clustered regions.

        The method checks if the data has already been mapped to the clustered regions. If not, it performs the mapping by calling the `_map_to_clustered_regions` method. It also checks if the energy data has been converted to the default energy unit. If not, it performs the conversion by calling the `_convert_energy` method.

        Returns
        -------
        gpd.GeoDataFrame
            A GeoDataFrame containing the data mapped to the clustered regions, with energy data converted to the default unit if necessary.
        """
        if self._mapped:
            pass
        else:
            self._data_in_regions = self._map_to_clustered_regions()
            self._mapped = True
        if not self.converted:
            self.converted = True
            self._data_in_regions = self._convert_energy()

        return self._data_in_regions

    def _map_to_clustered_regions(self):
        """
        Maps the data to the clustered regions based on their geographical locations.

        This method performs a spatial join between the `data` GeoDataFrame and the `clustered_regions` GeoDataFrame. The join is done based on the 'within' predicate, meaning that data points are mapped to a clustered region if they are geographically within that region.

        After mapping, it aggregates the data by the 'name' field of the clustered regions, summing up the values of the `data_column` for each region. This aggregated data is then set as a new column in a copy of the `clustered_regions` GeoDataFrame.

        Returns
        -------
        gpd.GeoDataFrame
            A GeoDataFrame containing the original geometries of the clustered regions, with an additional column for the aggregated data from the `data_column`.
        """
        data_in_regions = gpd.sjoin(
            self.data, self.clustered_regions, how="right", predicate="within"
        )

        # Initialize an empty list to store the merged GeoDataFrames
        ret_val = self.clustered_regions.copy()

        # Aggregate the data
        ret_val[self.data_column] = (
            data_in_regions.groupby("name")[self.data_column]
            .sum()
            .reset_index(drop=True)
        )

        return ret_val

    def combine(self, others: dict, own_key: str) -> gpd.GeoDataFrame:
        """
        Combines the heat potentials from multiple ClusteredRegionData instances into a single GeoDataFrame.

        This method takes the current instance's data mapped to clustered regions and combines it with the data from other ClusteredRegionData instances provided in the `others` dictionary. Each instance's data is identified by a unique key in the resulting GeoDataFrame.

        Parameters
        ----------
        others : dict
            A dictionary where keys are strings representing the column names for each ClusteredRegionData instance's data in the combined GeoDataFrame, and values are the ClusteredRegionData instances to be combined.
        own_key : str
            The column name for the current instance's data in the combined GeoDataFrame.

        Returns
        -------
        gpd.GeoDataFrame
            A GeoDataFrame containing the combined data from the current instance and the instances provided in `others`. Each instance's data is in its own column, identified by the corresponding key.

        Raises
        ------
        ValueError
            If any value in `others` is not an instance of ClusteredRegionData.
            If the geometries of the clustered regions in any of the ClusteredRegionData instances do not match.
        """

        ret_val = self.data_in_regions.copy()
        ret_val.rename(columns={self.data_column: own_key}, inplace=True)

        for key, val in others.items():
            ...
            if not isinstance(val, ClusteredRegionData):
                raise ValueError(
                    "All elements in the list must be of type ClusteredRegionData"
                )
            if any(
                val.clustered_regions.geometry.to_numpy()
                != self.clustered_regions.geometry.to_numpy()
            ):
                raise ValueError(
                    "All elements in the list must have the same clustered_regions"
                )

            ret_val[key] = val.data_in_regions[val.data_column]

        return ret_val

    def _convert_energy(
        self,
        conversion: dict = {"Wh": 1, "kWh": 1e3, "MWh": 1e6, "GWh": 1e9, "TWh": 1e12},
    ) -> gpd.GeoDataFrame:
        """
        Convert the energy unit of the data to the desired energy unit.

        Parameters:
            conversion (dict): A dictionary that maps energy units to conversion factors.
                The default conversion dictionary is {"Wh": 1, "kWh": 1e3, "MWh": 1e6, "GWh": 1e9, "TWh": 1e12}.

        Returns:
            gpd.GeoDataFrame: The converted data in a GeoDataFrame.

        Raises:
            ValueError: If the energy unit is not found in the conversion dictionary.
        """
        try:
            self._data_in_regions[self.data_column] = (
                self._data_in_regions[self.data_column]
                * conversion[self.energy_unit]
                / conversion[self.default_energy_unit]
            )
            return self._data_in_regions
        except KeyError:
            raise ValueError("Energy unit not found in the conversion dictionary")


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_heat_sources",
            simpl="",
            clusters=48,
        )

    clustered_regions = gpd.read_file(snakemake.input.regions_onshore)

    set_scenario_config(snakemake)

    clustered_heat_sources_dict = {
        heat_source: ClusteredRegionData(
            clustered_regions=clustered_regions,
            data=gpd.read_file(
                snakemake.params.heat_sources[heat_source]["source_data"]
            ),
            data_column=snakemake.params.heat_sources[heat_source]["data_column"],
            energy_unit=snakemake.params.heat_sources[heat_source]["unit"],
        )
        for heat_source in snakemake.params.heat_sources
    }

    heat_source_potentials = list(clustered_heat_sources_dict.values())[0].combine(
        {
            key: val
            for key, val in clustered_heat_sources_dict.items()
            if key != list(clustered_heat_sources_dict.keys())[0]
        },
        own_key=list(clustered_heat_sources_dict.keys())[0],
    )

    heat_source_potentials.index = heat_source_potentials.name
    heat_source_potentials = heat_source_potentials.drop(columns="name")
    heat_source_potentials.to_csv(snakemake.output["heat_source_potentials"])

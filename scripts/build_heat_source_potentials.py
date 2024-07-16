import geopandas as gpd
from typing import List
from _helpers import set_scenario_config


class ClusteredRegionData:

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
        This function maps the heat potentials to the Clustered regions


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
        This function combines the heat potentials in the regions

        """

        ret_val = self.data_in_regions.copy()
        ret_val.rename(columns={self.data_column: own_key}, inplace=True)

        for key, val in others.items():
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
        This function converts the energy unit of the data to the desired energy unit

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

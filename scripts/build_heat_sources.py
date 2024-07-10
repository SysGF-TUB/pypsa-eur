import geopandas as gpd
from typing import List
from _helpers import set_scenario_config


class ClusteredRegionData:

    def __init__(
        self,
        clustered_regions: gpd.GeoDataFrame,
        data: gpd.GeoDataFrame,
        data_columns: List[str],
        mapped: bool = False,
    ) -> None:
        self.clustered_regions = clustered_regions
        self.data = data.to_crs(clustered_regions.crs)
        self.data_columns = data_columns

        self._mapped = mapped
        if mapped:
            self._data_in_regions = data

    @property
    def data_in_regions(self) -> gpd.GeoDataFrame:
        if self._mapped:
            return self._data_in_regions
        else:
            self._data_in_regions = self._map_to_clustered_regions()
            self._mapped = True
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

        for column in self.data_columns:
            # Aggregate the data
            ret_val[column] = (
                data_in_regions.groupby("name")[column].sum().reset_index(drop=True)
            )

        return ret_val
    
    def combine(self, others:dict, own_key: str):
        """
        This function combines the heat potentials in the regions

        """
        
        ret_val = self.data_in_regions.copy()
        ret_val.rename(columns={col: f"{own_key}_{col}" for col in self.data_columns}, inplace=True)
        data_columns = [f"{own_key}_{col}" for col in self.data_columns]

        for key, val in others.items():
            if not isinstance(val, ClusteredRegionData):
                raise ValueError("All elements in the list must be of type ClusteredRegionData")
            if any(val.clustered_regions.geometry.to_numpy() != self.clustered_regions.geometry.to_numpy()):
                raise ValueError("All elements in the list must have the same clustered_regions")

            for col in val.data_columns:
                ret_val[f"{key}_{col}"] = val.data_in_regions[col]
                data_columns.append(f"{key}_{col}")

        return ClusteredRegionData(clustered_regions=self.clustered_regions, data=ret_val, data_columns=data_columns, mapped=True)

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
            data=gpd.read_file(snakemake.params.heat_sources[heat_source]["source_data"]),
            data_columns=snakemake.params.heat_sources[heat_source]["data_columns"],
        )
        for heat_source in snakemake.params.heat_sources
    }

    clustered_heat_sources_combined = list(clustered_heat_sources_dict.values())[0].combine(
        {
            key: val for key, val in clustered_heat_sources_dict.items() if key != list(clustered_heat_sources_dict.keys())[0]
        },
        own_key=list(clustered_heat_sources_dict.keys())[0],
    )
    clustered_heat_sources_combined.data_in_regions.to_csv(snakemake.output["heat_sources"])
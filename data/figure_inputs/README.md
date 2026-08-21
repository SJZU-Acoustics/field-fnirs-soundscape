# Figure inputs for Figure 1

Two inputs used only by `code/fig1_design.R`. Neither is a measurement of this
study, so they ship with the code rather than in the Mendeley deposit.

## `p27_basemap_osm.csv` — OpenStreetMap base geometry (2,956 rows)

Map geometry around the site transect (water, park, buildings, roads, footways)
as WGS-84 polygon/path vertex sequences, one row per vertex with `layer`,
`feature_id` and `seq` columns. Frozen from OpenStreetMap by the working
pipeline's `build_site_basemap.py` (Overpass API query, August 2026).

© OpenStreetMap contributors, available under the **Open Database License
(ODbL)** — see <https://www.openstreetmap.org/copyright>. Redistribution of
this CSV as a Produced Work is permitted with attribution; if you adapt it and
publish the result, ODbL share-alike applies to the adapted database.

## `p27_montage_layout.csv` — digitised fNIRS montage (30 rows)

Source/detector positions (head-surface coordinates), channel assignments,
regions and the retained/excluded signal-quality status of the 10 channels.
Digitised from the data collector's montage record figure by the working
pipeline's `digitise_fnirs_montage.py`; authors' own work.

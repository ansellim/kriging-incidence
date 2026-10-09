# Data and provenance

All inputs required by the offline analysis are included. Coordinates in the site table identify study locations, not individual case addresses.

| File | Source and role |
|---|---|
| `raw/assignment2.csv` | Thirty-site table transcribed from SPH5105 Individual Assessment 2 and checked against the assessment PDF. The table has 56 cases in total. No disease name, denominator, observation period, or catchment definition was supplied. |
| `cache/singapore_land.geojson` | [URA Master Plan 2019 Subzone Boundary (No Sea)](https://data.gov.sg/datasets/d_8594ae9ff96d0c708bc2af633048edfb/view). Administrative land proxy and planning-area labels; not a surveyed current coastline. |
| `cache/community_clubs.geojson` | [People's Association Community Clubs (2024)](https://data.gov.sg/datasets/d_f706de1427279e61fe41e89e24d440fa/view). Original cached government layer. |
| `cache/community_clubs.csv` | Derived club names, addresses, and WGS84 coordinates from the same PA layer, used by the analysis. |
| `cache/mrt_exits.geojson` | [LTA MRT Station Exits](https://data.gov.sg/datasets/d_b39d3a0871985372d7e1637193335da5/view). Cached catalog version updated September 16, 2026; the catalog describes data from August 2025. |
| `cache/facility_verification.csv` | Dated checks against the official onePA community-club directory, with URLs, address agreement, accessibility notes, and limitations. Checks performed October 7 and 8, 2026. |
| `cache/source_catalog.json` | Provider identifiers, vintage notes, feature counts, and cache date. |

The government layers were accessed on October 7, 2026 and are available under the [Singapore Open Data Licence](https://data.gov.sg/open-data-licence). This repository acknowledges URA, the People's Association, and LTA as the sources and has no official status or agency endorsement. The assessment table's provenance is distinct from the government layers; no independent reuse license is asserted for it.

## Site table fields

| Column | Meaning |
|---|---|
| `id` | Site identifier, 1 to 30. |
| `latitude`, `longitude` | WGS84 coordinates in decimal degrees (EPSG:4326). |
| `x`, `y` | SVY21 / Singapore TM coordinates in meters (EPSG:3414). |
| `green_space` | Green-space percentage; the buffer size and cover definition are unspecified. |
| `cases` | Nonnegative integer case count, 0 to 6. |

The analysis uses projected `x` and `y` for distance calculations and checks them against the supplied latitude and longitude. It dissolves the URA polygons for clipping, filters community clubs, removes historical or unrelated venues, and treats outdated construction flags according to the dated verification record. These choices are documented in the report and code. Cached inputs are preserved so later government updates do not silently change the results.

# Kristineberg (Gullmar fjord) sea-surface temperature

Hourly sea-surface (and air) temperature from the **Kristineberg Marine
Research Station** automatic weather station, on the Gullmar fjord, Swedish
west coast (Skagerrak) — the field-temperature record for the brown shrimp
(*Crangon crangon*) panel of manuscript Figure 5.

## Files

- `kristineberg_sea_temp_hourly.csv.gz` — the analysis-ready file. Hourly means
  of the continuous sea-temperature sensor, 2008–2024.
  - `datetime_utc` — hour (UTC).
  - `year` — calendar year.
  - `sea_temp_c` — sea-surface temperature (°C), hourly mean of the 5-minute
    `SeaTemp` sensor.
  - `air_temp_c` — air temperature (°C), hourly mean.
- The raw 5-minute, all-variable archive (`../../../data-raw/kristineberg/`
  `kristineberg_field_temp_5min_all.csv.gz`) and the script that builds both
  files (`data-raw/build_kristineberg_temps.R`) live outside the package.

## Source

University of Gothenburg / Sven Lovén Centre, Kristineberg:
<https://www.weather.mi.gu.se/kristineberg/en/data.shtml>. Data were fetched via
the station's `result.php` archive form (5-minute CSV, all 12 selectable
fields), stacked, cleaned, and averaged to hourly means. Processing follows
P. Pottier's earlier column mapping (`no data` → `NA`; drop lab-sensor spikes
with `surface_water_temp < 100`).

## Data-use caveat

The station's data page states these measurements **lack formal quality
assurance and may contain errors**, and asks users to **contact the station's
technical administrator before any use in a publication**. This must be done
before Figure 5 is finalised for submission.

## Use in Figure 5 (important)

The Gullmar fjord surface only just reaches this shrimp population's thermal
limits: across all 17 years (2008–2024) the maximum hourly sea-surface
temperature in this file is **25.4 °C**, and the warmest summer by peak hourly
temperature (**2021**) reaches **24.8 °C**, while the corrected shrimp fit gives
**T_crit ≈ 24.9 °C** and **CTmax₁ₕ ≈ 31.6 °C**. A literal trace therefore
accumulates little or no heat injury. The manuscript heat-injury figure instead
uses the warmest summer on record by peak temperature (**2021**) with an
explicit **+3 °C warming projection** so the trace clearly crosses T_crit —
labelled as a projection in the figure caption, not as observed lethal exposure.

(These T_crit / CTmax values are the corrected fit after the 2026-06 shrimp
data-processing fix: `Mortality_after_trial` is a proportion, so survivor counts
are `round((1 - mortality) * N)`. The earlier figures used 26.1 °C / T_crit
27.7 °C, which were artefacts of the mis-processed mortality column and of an
older temperature extract; do not reuse those numbers.)

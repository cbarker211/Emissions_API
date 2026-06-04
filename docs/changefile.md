For a full list of changes, please visit https://github.com/cbarker211/Emissions_API

# v2.6.0

This updates adds launch emissions in atmospheric layers, fixes bugs, and streamlines the data.

## Launches
- Launch emissions now output in atmospheric layers (0-15, 15-50, 50-80, >80).
- Fixed issue with Falcon Heavy emissions, which were previously included for boosters only and always included regardless of whether they landed or not.
- Added booster propellant mass for Atlas D.
- Added Chang Zheng 2C SM upper stage.

## Re-entries
- Removed re-entries where the object exploded in orbit, and re-entries with re-entry destinations of "lost" and "Z".
- Removed some debris objects, most of which had zero mass anyway.
- Three new re-entries added in 2021 and 2025.
- Classified more objects as fairings where name includes "PLF", slightly decreasing emissions.
- Sorted out stage classification for some upper stages which were incorrectly being assigned 0 mass.
- Added Chang Zheng 2C SM upper stage.

## QOL additions with no impact on emission totals
- Minor QOL fixes to clean up the scripts, including removing discosweb related functions and removing unused packages.
- Combined date and time fields into single datetime objects.
- Raul data now loaded as a pandas dataframe.
- Updated re-entries with a location of '-1'. This was used a placeholder to denote named regions where I needed to assign the geolocation.
- OV1-7P re-classified as payload.

# v2.5.0

This major update includes an extension of re-entry emissions back to 1957. 

## Re-entries

- Only checking missing stages with an apogee of >100 km.
- Skipping objects from explosions.
- Not cross-checking against DISCOSweb or Aerospace Corp.
- Not including objects from failed launches pre-1957.
- Not checking for missing fairings.
- Not updating timings for re-entries with missing timings.
- Minor changes to re-entry mass for 2020-2024 from updates to rocket stage masses.

## Refactoring

- Now loading all information for launches and rockets using GCAT, and merging into single files across the space age.
- This has updated the timings, locations, and launch site names, but this won't affect the emissions. The coordinate changes are only small and will be in the same grid box.
- The names and variants are slightly different in GCAT, so the new names needed to be added to launch event altitudes and update_rocket_launch_data.py script.
- GCAT also uses slightly different designations for failed launches, which are now accounted for.
- Removed lat lon from re-entry outputs, as these aren't even used on the website.

## Bug fixes to launch emissions

Changing the names to the GCAT names uncovered missing launch event altitudes and incorrect propellant data which was fixed. Overall there were changes in emissions for:
- Changes pre-2020:
    - Pegasus XL/HAPS
    - Minotaur I
    - Proton-M/Briz-M, Protom-M/DM-2
    - PSLV CA
    - H-IIA 202 1 variant
    - Minotaur IV variants
    - Soyuz ST-A/B Fregat
    - Epsilon variants.
    - Atlas V XXX 'C' variants
    - Vega 1 variant
    - H-IIA 204 2 variant
    - Delta 4H A variant
    - GSLV Mk III a variant
    - Electron
    - GSLV Mk II b/c variant
    - LVM3
    - Ariane 5ECA+ variant
    - Antares 230+ variant
    - CZ 2C and variants 
    - CZ 3B/C and variants 
    - Minotaur I/IV and associated variants
    - Epsilon.
- Changes post-2020: 
    - CZ 2C and variants 
    - CZ 3B/C and variants 
    - Minotaur-I
    - CZ 7A (above only)
    - CZ 8
    - Shuang Quxian 1 3 variant
    - 2023-201 Soyuz-2-1B has changed to Soyuz-2-1A in GCAT
    - 2024-240 Chang Zheng 5/YZ-2 has changed to Chang Zheng 5B/YZ-2.
    - New Starship launches 2024-U01, 2024-U03, 2024-U04, 2024-U06
    - Kairos and Ariane 62 data is now from GCAT so the emissions are slightly different.
- GCAT only has Fregat, no Fregat-M or Fregat-MT (this only affects the third stage which is above the model).

# v2.0.0

This major update contains an extension of launch emissions for 1957-2024, as well as code refactoring and bug fixes.

## Emissions

### Changes affecting emission totals
- Reworked how failed launches are handled. Fixed a bug for rocket launches where stage 1 failed during ascent. These were previously zeroed **below** the cutoff altitude instead of **above**. This affects one launch in 2020 and two in 2021.
- Fixed a bug where SECO of Starship was assigned incorrectly, leading to an underestimate of emissions (affects one launch in 2023).
- Minor adjustments to syntax throughout. No changes to results.
- Fixed a bug where the first 9 levels of the fine_grid_mass array were wrongly set to zero during a failed launch on 12/09/2020, affecting launches in 2020 only from 12/09/20 onwards (49 launches). This meant that the launch emissions at the surface were being underestimated, and will particularly affect NOx and CO2 (largest emissions close to surface).
- Tweaks to Falcon launches to correctly assign the re-entry method for each individual stage (ocean platform, ground landing pad, expended). Previously all were assumed to be reusuable.
- Slight tweak to altitude of FEI for Pegasus.

### QOL additions with no impact on emission totals
- Huge improvement in speed (4x) due to refactoring of code.
    - Minimised np.sum calls.
    - Input files loaded as pandas dataframes.
    - Moved the calculation of the grid to outside the main day-month loop. This will need to be changed back if using meterological box heights.
    - Refactored the automatic assignment of launch event altitudes.
    - Refactored the calculations of fine_grid_mass_stage_{i} (fp changes to results).
- Removed hardcode of launch IDs for lat/lon assignment for Naro Space Center and China Sea Launch. Now automatic based on name, very small change to coordinates of one launch in 2023.
- Small adjustments to prepare the scripts for launches back to 1957.
- Added ignition altitudes for air-launched rockets back to 1957 (Pegasus variants).
- Edited output files to a fixed 6 decimal places.
- Added emission totals above the model.

## Launches
- 2021-076 and 2022-056 now marked as containing megaconsellation payloads due to the Chinese National Constellation.
- Changing launch site names so that launches from the same spaceport (eg Hainan Commercial Spaceport,Wenchang Satellite Launch Center) are set as the same, so we can group them on the tracker. Coordinates remain the same so no change to emissions.

## Rockets

### Changes affecting emission totals 
- Minor fix to Long March (CZ) 3B/YZ-1 upper stage masses.
- Fixed bug where GSLV Mk II and H-III 22 had the stages assigned incorrectly. 
- DISCOSweb has updated 2020-008,2020-015,2020-018,2020-020,2020-031,2020-068,2020-075 from Soyuz-2-1B Fregat-M to Soyuz-2-1B Fregat.
- Fixed bug where booster/stage fuel types were reversed for GSLV Mk II.
- Fixed incorrect fuel types for Long March (CZ) 2C in 2023/2024.
- Added missing stage mass for Long March (CZ) 5/YZ-2.
- Added missing fuel/stage mass for Long March (CZ) 7/YZ-1A.
- Atlas V N22 and Atlas V 042 (variant 1) corrected booster to use AJ-60A, not GEM63. 
- Removed 4th stage of Shavit, unsure whether this is included?

### QOL additions with no impact on emission totals
- Added masses and proxies for missing mass info in GCAT. 
    - Where launch mass is given but not dry mass, we assume a mass ratio of 10.
    - E.g. for a launch mass of 100, dry mass is 10.
- Rocket attributes now include a 5th stage, -1 stage and variant.
- Updated Shavit to include Shavit, Shavit 1 and Shavit 2.
- Zhuque-2 renamed to Zhuque-2E from 2024 onwards.
- Falcon 9 v1.2 renamed to Falcon 9 FT5 throughout.

# v1.5.0

This version contained an update to include 2023 and 2024 emissions.

- Various updates to speed up the scripts.
- Added BECO for Angara A5 Orion.
- Databases from Jonathan McDowell's GCAT are now loaded in automatically from the URL instead of from a downloaded static file.
- Added reentry black carbon, HCl, and Cl.
- Added new megaconstellations to list.
    - Lingxi (from Yinhe)
    - Protosat-1 (from E-Space)
    - Tranche (US SDA)
    - Amazon Kuiper
     - Chinese National Constellation (Guowang/Xingwang/Guangwang/Hulianwang)
    - Qianfan
- Fixed re-entry location of launches from 'Naro Space Center', 'China Sea Launch', 'Newquay, Spaceport Cornwall', and 'Mojave Air and Space Port'.

# v1.0.0 

This version is detailed in https://www.nature.com/articles/s41597-024-03910-z and contains emissions from 2020-2022. 
The original repository can be found here https://github.com/cbarker211/Satellite-Megaconstellation-Emission-Inventory-Development.
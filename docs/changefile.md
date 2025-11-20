For a full list of changes, please visit https://github.com/cbarker211/Emissions_API

# v2.0.0

This major update contains an extension of emissions for 1957-2024, as well as code refactoring and bug fixes.

## Emissions

### Bug fixes affecting emission totals
- Reworked how failed launches are handled. Fixed a bug for rocket launches where stage 1 failed during ascent. These were previously zeroed **below** the cutoff altitude instead of **above**. This affects one launch in 2020 and two in 2021.
- Fixed a bug where SECO of Starship was assigned incorrectly, leading to an underestimate of emissions (affects one launch in 2023).
- Minor adjustments to syntax throughout. No changes to results.
- Fixed a bug where the first 9 levels of the fine_grid_mass array were wrongly set to zero during a failed launch on 12/09/2020, affecting launches in 2020 only from 12/09/20 onwards (49 launches). This meant that the launch emissions at the surface were being underestimated, and will particularly affect NOx and CO2 (largest emissions close to surface).
- Tweaks to Falcon Heavy launches to improve handling of booster/stage1 after use. Fixes launches where all boosters and core are expended (one in 2023, one in 2024).
- Slight tweak to altitude of FEI for Pegasus.

### QOL fixes with no impact on emission totals
- Moved the calculation of the grid to outside the main day-month loop to increase speed. This will need to be changed if using meterological box heights.
- Refactored the automatic assignment of launch event altitudes to improve readability and prevent errors. 
- Completely refactored the calculations of fine_grid_mass_stage_{i}. Floating point changes to results.
- Removed hardcode of launch IDs for lat/lon assignment for Naro Space Center and China Sea Launch. Now automatic based on name, very small change to coordinates of one launch in 2023.
- Small adjustments to prepare the scripts for launches back to 1957.
- Added ignition altitudes for air-launched rockets back to 1957 (Pegasus variants and NOTS EV1).

## Launches
- 2021-076 and 2022-056 now marked as containing megaconsellation payloads due to the Chinese National Constellation.
- Changing launch site names so that launches from the same spaceport (eg Hainan Commercial Spaceport,Wenchang Satellite Launch Center) are set as the same, so we can group them on the tracker. Coordinates remain the same so no change to emissions.

## Rockets
- Added masses and proxies for missing mass info in GCAT.
- Rocket attributes now include a 5th stage and variant.
- Minor fix to Long March (CZ) 3B/YZ-1 upper stage masses.
- Updated Shavit to include Shavit, Shavit 1 and Shavit 2.
- Fixed bug where GSLV Mk II and H-III 22 had the stages assigned incorrectly. 
- DISCOSweb has updated 2020-008,2020-015,2020-018,2020-020,2020-031,2020-068,2020-075 from Soyuz-2-1B Fregat-M to Soyuz-2-1B Fregat.
- Zhuque-2 renamed to Zhuque-2E from 2024 onwards.
- Fixed bug where booster/stage fuel types were reversed for GSLV Mk II.
- Fixed incorrect fuel types for Long March (CZ) 2C in 2023/2024.
- Added missing stage mass for Long March (CZ) 5/YZ-2.
- Added missing fuel/stage mass for Long March (CZ) 7/YZ-1A.
- Atlas V N22 and Atlas V 042 (variant 1) corrected booster to use AJ-60A, not GEM63. 
- Removed 4th stage of Shavit, unsure whether this is included?

# v1.5.0

This version contained an update to include 2023 and 2024 emissions.

- Fixed inc
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
- Fixed re-entry location of launches from 'Naro Space Center', 'China Sea Launch', 'Newquay, Spaceport Cornwall', and 'Mojave Air and Space Port.

# v1.0.0 

This version is detailed in https://www.nature.com/articles/s41597-024-03910-z and contains emissions from 2020-2022. 
The original repository can be found here https://github.com/cbarker211/Satellite-Megaconstellation-Emission-Inventory-Development.
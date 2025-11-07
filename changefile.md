For a full list of changes, please visit https://github.com/cbarker211/Emissions_API

# v1.0.0

This update contains code refactoring and bug fixes.

- Reworked how failed launches are handled. Fixed a bug for rocket launches where stage 1 failed during ascent. These were previously zeroed **below** the cutoff altitude instead of **above**. This affects one launch in 2020 and two in 2021.
- Small adjustments to prepare the emissions calculation script for launches back to 1957. No changes to results.
- Moved the calculation of the grid to outside the main day-month loop to increase speed. This will need to be changed if using meterological box heights. No changes to results.
- Refactored the automatic assignment of launch event altitudes to improve readability and prevent errors. Fixed a bug where SECO of Starship was assigned incorrectly, leading to an underestimate of emissions (affects one launch in 2023).
- Completely refactored the calculations of fine_grid_mass_stage_{i}. Floating point changes to results.
- Minor adjustments to syntax throughout. No changes to results.
- Removed hardcode of launch IDs for lat/lon assignment for Naro Space Center and China Sea Launch. Now automatic based on name, very small change to coordinates of one launch in 2023.
- Fixed a bug where the first 9 levels of the fine_grid_mass array were wrongly set to zero during a failed launch on 12/09/2020, affecting all launches from 12/09/20 onwards (49 launches). This meant that the launch emissions were being underestimated, and will particularly affect NOx and CO2 (largest emissions close to surface).

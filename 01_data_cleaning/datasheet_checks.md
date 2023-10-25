# Things to check in data sheets

Lomatium data cleaning notes. Things in this list should be checked in the raw data.

## October 24 2023

### 2021 demo multi-umbel plants

Entered as having more than one umbel (usually two, often four) with one umbel diameter listed. Not listed in the multi-umbel dataset that year. Unclear if this was umbel number being mis-entered or if simply the max umbel was listed (according to the 2021 protocol).

After checking, go back to the `combine_all_data.R` script (where the demo and multi-umbel datasheets are reconciled).

##### tag 3135 (plot 12, 7G)
##### tag 3143 (plot 9, 14-0.95)
##### tag 3163 (plot 2, 11D)
##### tag 3352 (plot 2, 19H)
##### tag 3360 (plot 6, 12C)
##### tag 3543 (plot 1, 2F)
##### tag 3567 (plot 1, 3I)
##### tag 3584 (plot 2, 7F)
##### tag 3675 (plot 13, 13I)
##### tag 3736 (plot 15, 14B)
##### tag 3814 (plot 1, 0J)
##### tag 3868 (plot 15, 15A)
##### tag 3978 (plot 13, 2F)

Less important, but there are a very small number of cases of SOS plants in 2022 in demo that didn't get entered. Maybe check the datasheet for notes while you're there.

### 2020 demo plant umble numbers

In 2020 a bunch of plants were recorded as having >10 umbels (no diameter, no stalk height). Very clearly some kind of mistake. Perhaps this should be number of flowers and was a mistake made in the field? Maybe some data entry error?

Check the following on the 2020 data sheets:

##### tag 3749 (plot 2, 6B)
##### tag 3751 (plot 2, 8A)
##### tag 3786 (plot 4, 6H)
##### tag 3787 (plot 4, 3G)
##### tag 3788 (plot 5, 0F)
##### tag 3747 (plot 6, 18A)
##### tag 3784 (plot 6, 18F)
##### tag 3785 (plot 6, 15G)
##### tag 3790 (plot 6, 5I)
##### tag 3791 (plot 6, 9J)

When done go to the `combine_all_data.R` script (the part assessing plants `in.demo` but not `in.mumb`).


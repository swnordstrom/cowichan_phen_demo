# Data cleaning assumptions/decisions

## Demo cleaning

##### `umble.no` field greater than 1 but umble diameter listed

## Demo and seed comparison

##### Sampling in 2021

Were plants sampled for seed production in 2021 conditioned on not having umbel death/florivory? Looking at the percentages of zeros (by plot) across the years, it looks like there are more zeros in 2022-2023 than 2021. This could totally be natural but, if estimating the probability of umbel failure (in a zero-inflated poisson model) then we'd want to know this. 

##### Lumping plants based on notes

For plants where it goes missing (NP/NPNT/leaf count zero) and a different plant appears nearby, the year after the first absence, I'll collect them if there's a note saying to do so. E.g., plants 3456 and 7587 (Plot 6, 13A/B). 7587 is listed as at location A, but I bet it's on the border.

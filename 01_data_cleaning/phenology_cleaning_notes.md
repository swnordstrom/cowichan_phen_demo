# Phenology (resurvey) data cleaning/processing notes

### Goals

Objective of this process is to get an initiation of flowering date (thinking right now this variable will be called `y` in the IPM formulation) for use in models predicting seed set (and other transitions, likely).

The data for 2021-2023 has individual-level records, with each plant monitored approximately once-weekly (still assessing this). The cycle of flower and seed production is discretized into stages (unfortunately the stages are not consistent - a different set of stages was used in 2021 compared to 2022-2023). For each plant on each date (presumably each date with a non-zero number of features in at least one stage), the number of objects on the plant in each stage (e.g., bud) was recorded. 

Note that we are interested for now only in the initiation of flowering for each plant. We could possibly extend this to the initiation of flowering for each inflorescence (such that multiple inflorescences might produce multiple flowering days for each plant) but we'll save that for later. For now, if we only want the initiation of flowering for each plant; this is simply the earliest date where the plant has anything in the earliest stage we deem possible for flowers to shed and receive pollen.

### Stage descriptions

In 2021 the following stages were used:

* closed bud
* unfurling bud
* partially-open umbel
* open and flat umbel (umblets are not fully spread and make a flat surface)
* open (umblets are spread, possibly making a curved surface)
* open and losing petals
* seeding
* "dead" - brown, desicated, not producing seeds

Note that all of the stages except "dead" form a sequence, whereas "death" can follow any stage.

In 2022 the following stages were used:

* buds
* pods (this is distinct from pods, but I am not entirely sure how - it must be a visual distinction)
* flowers
* flowers losing petals
* seeding
* dead
* eaten

As above it's likely that "dead" or "eaten" can occur after any stage.

The stages used in 2023 are the same as in 2022, although frustratingly on the data sheet the "pods" and "flowers" stages are swapped in the order of columns. Here, "pods" appears between "flowers" and "flowers losing petals" - I'm going to assume this is a mistake.

##### Decisions to make

* Which stage in 2021 should we classify as "receiving pollen"? The 2021 metadata file says that plants at the "partial" stage can "probably be pollinated". If not this, then umbels in the "flat" stage surely are capable of shedding and receiving pollen

* Stages in 2022-2023: I'm going to assume that "flowering" is the first stage where mating can occur. This assumes that the "pods" stage occurs before flowering and that plants in that stage are not yet reproductive

### Processing steps

First, need to get unique identifiers for each plant. Some tags are used multiple times, including possibly in the same plot. As such we should include plot and coordinates in the ID as well. I went ahead and did this; in 2022 the coordinates were not recorded, but there were no tags used more than once in the same plot (as far as I can tell). I'll store this info in the `plaid` column.

Next, I'll need to discretize the sampling dates into periods. I did this by visual inspection. note, though, that this approach does lose us some data: there are two plots in 2022 and one plot in 2023 that appeared to be missed for one sampling period. Two other irregularities: one plot was surveyed twice within one survey period in 2022 and two plots in 2023 were surveyed at the very end of the study in an additional survey period, but no other plots were studied.



# Social status mediates harvest-induced selection on behavioral type and predictability in Coyotes (_Canis latrans_)
### Nick A. Gulotta¹, Joseph W. Hinton², Michael J. Chamberlain¹
### 
####  ¹ Warnell School of Forestry and Natural Resources, University of Georgia, Athens, Georgia, USA. 
####  ² Wolf Conservation Center, South Salem, New York, USA

## Abstract
Human harvest is a pervasive selective force, yet how social structure shapes harvest-induced selection on behavior remains poorly understood. Using coyotes (_Canis latrans_), we evaluated whether harvest-induced selection acts on both mean behavioral expression (behavioral type) and residual intra-individual variation (behavioral predictability), and whether selection differs between dominant resident and subordinate transient individuals occupying contrasting life-history states. Selection on behavioral type operated in opposite directions according to social status. Residents occurring closer to hunter-associated landcover experienced reduced survival, whereas transients occurring closer to the same landcover experienced greater survival. Selection on behavioral predictability also differed by social status. Greater unpredictability increased survival among transients, whereas residents experienced stabilizing selection favoring intermediate levels of predictability. Our findings demonstrate that social status fundamentally alters both the direction and form of harvest-induced selection, highlighting the importance of accounting for life-history state when evaluating the behavioral and evolutionary consequences of harvest.

## General Project File Structure

```
├── README.md                                                          <- The top-level README including general project descriptions
|
├── Behavioral data
│   ├──TriState_daily_Landcover_FINAL.csv                              <- Behavioral data ready for modeling
|
├── Survival data
│   ├── Low_Arctic_productivity_gyrf.csv                               <- Survival data ready for modeling
|
├── Rds files
│   ├── hardwood_riiv_FINAL.zip                                        <- Fitted model file for hardwood landcover
│   ├── shrub_riiv_FINAL.zip                                           <- Fitted model file for shrub landcover
|
├── R script files
│   ├── DHGLM_BehaviorSurvival_Analysis.R                              <- Fitted model file for hardwood landcove



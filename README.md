# hba-nl
Pipeline commissioning for LOFAR2 HBA-NL

This github collects code related to HBA-NL pipeline commissioning efforts for LOFAR2.0

## Flagging fractions
Look at the fractions of data flagged to assess the overall quality.
This currently is focused on the calibrator - is it of sufficient quality?
Can look at both global statistics for an observation, e.g., mean/median of flagged data fraction.
Can also set a threshold for flagged data fraction per station and then see if a sufficient number of stations are of sufficient quality.
### Useage
``` 
import metrics.flagging_fraction as ff
ff.examine_global_flagging_metrics()
```
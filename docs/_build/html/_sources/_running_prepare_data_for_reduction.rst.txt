Prepare Data for Reduction
**************************

If you did a good job preparing and doing the observation this should be an easy
step, either way, keep in mind the following steps.

- Remove all *focus* sequence.
- Remove all *target acquisition* or *test* frames.
- Using your observation's log remove all unwanted files.
- Make sure all data has the same gain (``GAIN``) and readout noise (``RDNOISE``)
- Make sure all data has the same Region Of Interest or ROI (``ROI``).

The pipeline does not modify the original files unless there are problems with
fits compliance, is never a bad idea to keep copies of your original data though.



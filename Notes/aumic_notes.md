# AU Mic Research Notes: Feb 2017--
---

## To Do
-   dig around in bad flare timebin to find the source of the bad visibilities-- plot baseline length vs. amp, antenna number vs amp

## Questions
-   if clean process if fully automated, should final image rms be the same over multiple cleans w same parameters?

<br>

---
---

<br>

## Image Centering: 2/26/17

While comparing images made with different date combinations (i.e. removing August date because of poor quality), I noticed that the disk was offset from the image center for certain combinations. We have decided that this is caused by the non-homogeneous pointing centers (due to proper motion) of the three datasets. When `tclean` is called on a collection of datasets with different pointing centers, the pointing center of the first of the datasets is chosen as the origin of the image, and all datasets are combined in the *uv*-domain, with their phase offsets preserved. The resulting sky-domain image is both offset from the image center and a false representation of the disk.

To fix this issue, I tried using `contcat` with its *dirtol* parameter set to a high value (2"). As long as the pointing/phase centers of the datasets do not differ from each other by more than *dirtol*, the datasets are combined as if they all share the pointing center of the first dataset and all is well. A quick test indicates that this method is succesful.

However, I ran into another problem while attempting the `concat` method. The `26mar2014_aumic_spw0.corrected_weights.ms` dataset is missing `table.f8_TSM1` (all other datasets have this table), and because of this `concat` fails when applied to this dataset. 




**Fixvis:** `If the phase center is changed, the corresponding modifications are applied to the visibility columns given by the parameter "datacolumn" which is by default set to "all" (DATA, CORRECTED, and MODEL).`

-   All dates have 257, 257 `CRPIX`
    -   24 June: 311.291115417 -31.3424694444
-   24 June fixvis: `phasecenter = J2000 20h45m09.8677s -031d20m32.89s`
    -   not many precision points?
    -   311.2911154166666 -31.342469444444443 in deg... essentially the same as the image center
-   Center pixel info for all dates:
    -   aumic_18aug_usermask_natural:
        -   311.2910612917 -31.34236676111
    -   aumic_26mar_usermask_natural
        -   311.2910179167 -31.34232222222
    -   aumic_24jun_usermask_natural
        -   311.291115417 -31.3424694444
    -   aumic_usermask_natural
        -   311.2910612917 -31.34236676111
    -   aumic_marjune_usermask_natural
        -   311.291115417 -31.3424694444
March + June (the culprit):
    -   311.291115417 -31.3424694444
    -   257 257 center pixel

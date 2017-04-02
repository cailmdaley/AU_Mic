# AU Mic Research Notes: Feb 2017--
---

## To Do
-   check pixel location for FITS files (convention)
    - check how rotate rotates pixel
-   do percentage noise check
-   make august clean
-   compare gaussian fit to actual slice of data
-   recheck beam FWHM calculation
-   remove 10 min from spw3 only; restore other windows

## Questions
- $\chi^2$ went up for June dates, I suspect because of different weighting parameters? values are:
    - 2.735890291283209,
    - 3.1727607499540902,
    - 2.8817286537130427,
    - 2.0514430147058822,
    is this a problem?
- OK if imsize*cell != image size in arcsec?

<br>

---
---

<br>

## 3/21/17: Pixel location:

- ctrpix remains the same if I make image 257 pixels

CRPIXn from FITS standard:
>The value field shall contain a floating point number, identifying the location of a reference point along axis n, in units of the axis index.  This value is based upon a counter that runs from 1 to NAXISn with an increment of 1 per pixel.  The reference point value need not be that for the center of a pixel nor lie within the actual data array.  Use comments to indicate the location of the index point relative to the pixel.

From STSCI:
>  When the data matrix represents a digital image, transformation between the data matrix and the physical picture requires knowledge of where in the pixel -- center or corner -- the data point is. Historically, astronomers have generally assumed that the index point in a FITS file represents the center of a pixel. This interpretation is endorsed by GC. It differs from the common practice in computer graphics of treating the center of a pixel as a half-integral point. GC note that the pixel in a FITS file is commonly regarded as a volume element in physical space, which might be viewed from different perspectives via transposition and rotation. Under such operations, only the center of an element remains invariant. Pending adoption of a standard convention by the astronomical community, FITS writers should use appropriate comments in the comment field of the card image or the COMMENT keyword to inform readers of the file which convention is being used. Once the community has accepted a convention, a single comment noting that the convention is being used will be sufficient.

---
## 3/21/17: Flare date and bad spws

Recently I realized that the time window we split out to fix the bad spw in the June date was exactly the time window of the flare. This makes me somewhat suspicious, and Meredith and I decided I should do some more digging, espcially considering all the work we put into making the flare data useable.

The `plotms` of amp vs. time for spw3 (the bad one)  and spw1 (well behaved) are roughly the same--both show a huge spike in the last (flare) time window. This leads me to believe that it's not the flare itself that's screwing up spw3; if this were the case, we should see the same thing for spw1.

-   Antenna 1 and 2 are almost constantly 'on' in last time window, as opposed to dashed in previous windows?
-   same for baseline, phase
-   weights get very low for in flare window for both spws
- everything I tried seems to match for both spws...

This is a little confusing, since spw1 had a pretty nice $\chi^2$; but we did remove that flare time window for *all* spws...

---
<br>

## 2/26/17---3/20/17: Image Centering

While comparing images made with different date combinations (i.e. removing August date because of poor quality), I noticed that the disk was offset from the image center for certain combinations. We have decided that this is caused by the non-homogeneous pointing centers (due to proper motion) of the three datasets. When `tclean` is called on a collection of datasets with different pointing centers, the pointing center of the first of the datasets is chosen as the origin of the image, and all datasets are combined in the *uv*-domain, with their phase offsets preserved. The resulting sky-domain image is both offset from the image center and a false representation of the disk.

To fix this issue, I tried using `contcat` with its *dirtol* parameter set to a high value (2"). As long as the pointing/phase centers of the datasets do not differ from each other by more than *dirtol*, the datasets are combined as if they all share the pointing center of the first dataset and all is well. A quick test indicates that this method is succesful.

However, I ran into another problem while attempting the `concat` method. The `26mar2014_aumic_spw0.corrected_weights.ms` dataset is missing `table.f8_TSM1` (all other datasets have this table), and because of this `concat` fails when applied to this dataset. However, recreating the `.ms` file from the corresponding `.uvf` file seems to have fixed this problem. I'm starting over with a `cleans` directory, and have added the suffix `_old` to the originals.

---
***My troubleshooting notes:***

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

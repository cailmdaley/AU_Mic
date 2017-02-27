# Things to do

- [X] image flare date with flare
- [$\checkmark$] check fixvis inputs, check if we were using them correctly
- [$\checkmark$] check headers of individual files compared to combined (coordinates)
- [X] put star at actual pointing center, NOT center of image
- [X] rewrite table? fix sig figs

## Image Centering: 2/26/17
**Fixvis:** `If the phase center is changed, the corresponding modifications are applied to the visibility columns given by the parameter "datacolumn" which is by default set to "all" (DATA, CORRECTED, and MODEL).`

- All dates have 257, 257 `CRPIX`
  - 24 June: 311.291115417 -31.3424694444
- 24 June fixvis: `phasecenter = J2000 20h45m09.8677s -031d20m32.89s`
    - not many precision points?
    - 311.2911154166666 -31.342469444444443 in deg... essentially the same as the image center
- Center pixel info for all dates:
  - aumic_18aug_usermask_natural:
    - 311.2910612917 -31.34236676111
  - aumic_26mar_usermask_natural
    - 311.2910179167 -31.34232222222
  - aumic_24jun_usermask_natural
    - 311.291115417 -31.3424694444
  - aumic_usermask_natural
    - 311.2910612917 -31.34236676111
  - aumic_marjune_usermask_natural
    - 311.291115417 -31.3424694444
  - Question: if clean process if fully automated, should final image rms be the same over multiple cleans w same parameters? 

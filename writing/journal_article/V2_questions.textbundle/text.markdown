# 
![](assets/1EF6948F-DF8C-4EC5-B45C-2D8FBFE860EA.png)


## General things

### title

The word "Hidden" doesn't seem like the best choice to convey that this mass is responsible for exciting the vertical structure. Not sure what to suggest-- perhaps something more clinical like   
  
"The Mass of Stirring Bodies in the AU Mic Debris Disk inferred from Resolved Vertical Structure”.

### phase noise?

JC: Given the resolved thickness of the disk is rather narrow, can you discount uncalibrated phase noise as a significant impact? I suspect the answer is “no” at Band 6, but it may be worth address in the paper since it is not part of the model fitting.  
  
gain calibrator reduces phase noise… still some remaining though

- acknowledge as possibility

- look at ALMA/CASA technical manuals for phase noise stuff

- ask John what we should do

### observing wavelength?

mean of four centroid frequencies is 222.1375 Ghz; this implies a wavelength of 1.34958 mm…  
  
(if 1.3mm, make sure to change aumic_imaged)  
   
be explicit about how it was calculated

- change aumic_imaged

- fix sentence in observations

### significance of detection

- 0.031 +/- 0.005 implies 6 sigma  
- comparison with skinny model implies 4 sigma  
  
DW: ”4-6\sigma"? either pick one, e.g.  6 (as discussed in S4.3) and add something like "for our fiducial model", or use both but also indicate that this significance is model dependent?  
  
fix in introduction and conclusion!

- say 4 sigma

### should i be more trigger-happy with subsections in the discussion?

- yes

## Abstract

### can’t rule out Neptune analog?

MW: I was concerned about the conclusion that you can rule out a gas giant or Neptune analogue, and when I got to 5.3 I think this is not right (more below).

- near outer edge of disk

## Introduction

### Boccaletti Features

- Relevance the paper as a whole?
  Kevin:  
  I think this level of detail on the scattered light features isn’t really necessary in the introduction. The scattered light features are interesting, but don’t have much to do with the general goal of this project of trying to measure the scale height, and hence the mass of perturbing objects. You can condense this down into one or two sentences, and then discuss the features in more detail later when you look for them in the sub-mm data.   
    
  - some of Eugene’s comments/additions nicely tie our results to the discussion on the progenitors of the dust  
  - basically, are the moving features relevant enough to present in the introduction (possibly implying that they have relevance to the conclusions of the paper), or should they be saved for the results section where we say we don’t see the moving features?

	- ok to leave in introduction

- is Eugene’s language too biased towards his own work when discussing Boccaletti features?
  Notably, the debris around AU Mic exhibits a so-far-unique time variability at scattered light wavelengths.  
  \cite{boccaletti15} and \cite{boccaletti18} identify several local intensity maxima offset from the disk midplane on the southeast side of the disk.   
  These features are moving quickly away from the star along the disk midplane at projected velocities that are not consistent with Keplerian rotation; in fact, the outermost features appear to be unbound from the star (see also \citet{sezestre17} who provide kinematic fits and invoke   
  a dust source of unspecified nature).   
  \citet{chiang&fung17} propose that these fast-moving features  
  are dust particles repelled by a time-variable stellar wind that triggers  
  dust avalanches when the wind blows strongest.   
  These avalanches are seeded by the debris remaining from the recent disruption of a $\sim$400-km sized progenitor.   
  In this paper we will provide an independent constraint on the presence of comparably sized planetesimals.  
   

	- all ok

### spectral type dependence of Thebault results?

“Thus longer-wavelength (λ ≥ 50 µm) observations are required to measure disk scale height as determined by dynamical stirring alone, since the grains dominating the emission at these wavelengths are large enough to be impervious to the eﬀects of radiation pressure”  
  
JC: This must depend on the spectral type of the host star. You might want to give a bit more background on how the limit varies with SpT.  
  
Thebault: As for the orbital evolution computation, the important parameter is the value of β(s) for estimating the radiation pressure force. We thus chose to scale all particle sizes by their β(s) ∝ 1/s value… We consider here s min and s max , such that β(s min ) = 0.4 (close to the blow-out cut-oﬀ size) and β(s max ) = 0.025. This size range covers however most of the crucial grain population that dominates the ﬂux in scattered light.  
  
Cail: I couldn’t find any dependence on stellar properties in the Thebault paper. I think this is because particle sizes are simply described by their beta values (see below). Where these beta values come from isn’t clear…   
  
I would propose to leave things as is.

## Observations

### Flare

Cail: The real question here is: do we want to give more attention to the flare in this paper, or save it for somewhere/someone else? There’s clearly much more analysis we could do w/r/t the flare, but it would definitely take a little time.  
  
I would propose plotting the flare and looking at the polarization vs time (both should be very quick).  If something interesting pops up in the polarization, it might be worth expanding the flare discussion (maybe giving it its own subsection?)

- plot flare light curve?

- ‘deblend’ with uniform weighting?
  “Faithfully disentangling the two components proved diﬃcult, both because stellar and disk emission are blurred together in the central beam and because the stellar ﬂux varied signiﬁcantly across the three nights of observation (Table 3)”  
    
  JC: Can you deblend them a bit by using uniform weighting?  
    
  Cail: I don’t really know what this means, but I suspect that it won’t be better than fitting a point source directly to the visibilities. I would leave things as-is.

- look at polarization
  DW:  
  Regarding the flare, did you look at the polarization? In the Prox Cen data, there was a clear signature of polarization in a difference in XX and YY vs. time. It's not worth getting distracted by the flare, and maybe this analysis doesn't belong at all in this paper, or maybe it could be in an appendix if you have done the fits for XX and YY already.

### Miriad

SA: You claim to use MIRIAD (p4) in the calibration, but it doesn't seem like you do until you measure some emission.  My hope is it was only used for that!  Even so, why not just use the CASA viewer to do the same thing?  
  
Cail: Sean is right, I don’t think we ended up using MIRIAD at all in the data reduction. However, I would propose not re-measuirng the disk flux in CASA—I know you love cgcurs!

### self-calibration?

SA: The calibration discussion makes no mention of self-calibration.  Was that tried?  If not, why not?  If so, and it didn't help, it would be good to say so and report the interval(s) attempted.  [With a peak SNR ~ 20, I would think you could get a boost from self-cal.]  
  
Cail: I got nothing here. I don’t really know what self-calibration is, let alone why we didn’t do it. I guess I’ll ignore this comment unless you suggest otherwise..

## Results

### AU Mic images

DW: Figure 1, consider to separate the panels in order to make it easier for readers to lift the left panel only from the paper, to show in talks  
  
Cail: This comment made me realize that the tapered clean image has little-to-no function. Maybe we should take it out? I’ll only add some spacing between the panes, unless you think I should take out the tapered clean altogether.

### Boccaletti Plots

- Did Eugene misunderstand?
  From top to bottom, the three panes show i) the disk spine surface brightness with errors given by the rms noise, ii) the disk spine deviation from the midplane, and iii) the disk FWHM as a function of projected separation from the star.   
  The Gaussian traced by the dotted line in the first pane shows the projected width in the radial direction of the naturally weighted synthesized beam of the combined dataset.  
  %EIC: Following sentence ("The fact that ...") is not obvious to me just by looking at the figure.  
  %       In fact, the beam size as indicated by the Gaussian drawn on the top panel looks like it has a FWHM  
  %       of ~5 AU > the beam-subtracted FWHM plotted in the third panel. (I was pleased to see this  
  %       acknowledged in the main text, and to read about the various reasons why this is not necessarily  
  %       cause for alarm.)  
  %       I think the last sentence just needs to be reworded --- to acknowledge that the vertical beam  
  %       is technically larger than what is plotted, but that this is OK for various reasons discussed in  
  %       the main text under section 4.3. I gave it a stab.  
  %  The fact that the image-domain vertical height of the disk is in excess of the beam contribution implies that our data spatially resolve the vertical structure of the disk.  
  Although the last panel indicates a vertical disk FWHM that is technically smaller than the vertical FWHM of the beam, it still represents a significant detection of vertical structure, for reasons explained in \S \ref{subsection: vertical analysis}.}   
    
  Cail: This one’s tricky. Right off the bat, it seems like Eugene took the size of the beam in the radial direction to be the size in the vertical direction. That being said, the radial & vertical FWHMs aren’t that different so I don’t think this disqualifies his comments. It seems that he’s concerned by the fact that the beam-subtracted disk FWHM is smaller than the beam FWHM itself (this does strike me as odd, but I trust the math).  
    
  I would propose explicitly stating that the beam radial FWHM plotted in the top pane is *not* the beam vertical FWHM. Also comparing the beam-subtracted FWHM to the data resolution as Eugene suggested, and explaining that before subtracting the beam FWHM, the Gaussian FWHM was 8 au (or whatever it was exactly). However, I’m not really sure how to handle the sentence Eugene commented out and the one he replaced it with. On the one hand, the beam-subtracted disk FWHM is smaller than the beam, which might be grounds for calling it ‘unresolved’? On the other hand, the apparent disk FWHM (~8 au) exceeds that of the beam (~ 5 au), indicating that the disk is resolved. But wouldn’t this be true when convolving the beam with *any* Gaussian vertical profile with non-zero standard deviation? 

- uncertainty / fit information
  MW: What I am missing from this is a feeling for the typical uncertainties in these profiles, as well as how good a fit a Gaussian profile is to these vertical profiles.  
    
  Cail: I used astropy’s gaussian fitter, which doesn’t report uncertainties, to make these plots. It probably wouldn’t take too long to find & implement a fitter that does, and it would be really nice to have uncertainties on all three panes (this would also probably give a better estimate of the uncertainty for the top panel than just using the rms noise). To Mark’s second comment about whether the vertical profile is actually Gaussian, I could add another inlay to the top pane showing a sample vertical slice (or maybe the mean profile from all the slices?).  
    
  In the interest of getting the second draft out quickly, I would propose only adding a note stating that the vertical profile is well described by a Gaussian. But let me know if you’d like me to do either of the things suggested above—I think they’re good ideas.

- rotate by MCMC PA?
  DW was wondering how 0 elevation/deviation from the midplane was defined in the second pane. As I recall, I just did this by eye in the past. This is a little anachronistic, but I would propose using the MCMC best-fit PA (and explaining that we did so).

### Disk symmetry and eccentricity

HS:  
Also ion page 7, I suggest adding one sentence between ‘millimeter wavelength.’ and ‘Additional’ that explains what the absence of any significant brightness difference between the two limbs implies - and maybe even cite Margaret's paper on this. - Margaret, do you want to send one or two sentence about this?  
  
MP:  
Hilke mentioned the lack of significant brightness differences between the ansae and implications for the disk eccentricity (beginning of section 3). You have some info on the disk eccentricity from that but it depends on the disk's argument of periastron (if the argument of periastron is in the sky plane, the eccentricity would be the ratio of the fluxes minus unity, or 0.03 +/- 0.06, so consistent with zero). I could make a plot of the implied eccentricity as a function of the argument of periastron if you want to include a discussion of that - it would be something like "there is a 68% (or 95%, or 99%) chance the eccentricity is less than <a certain value>".  
  
Cail: This all seems like a great idea, especially as it doesn’t require any work on our end! Unless I hear otherwise, I’ll email Margaret asking her to make this plot, and what paper I should cite on the matter.  


### Disk flux discepancy

JC:  
From Table 1, the star(without the flare)  seems to have a flux of 0.45 mJy.   
  
However, you report a total disk+star flux of 4. 97 mJy, and a disk-only flux of 4.8 mJy. Thus I would infer a stellar flux of ~ 0.17 mJy.   
  
Is the apparent discrepancy because the 3 sigma contour is not enclosing all of the flux?  
  
  
  
  
Cail: This one has me a little worried. The pre-flare flux of 0.45 mJy presented in the table is very much not consistent with the MCMC June flux of 0.220 +/- 0.020 mJy. My guess is that we overestimated the pre-flare stellar emission when fitting a point source to the visibilities.   
  
Possibly a better way to do this would be to get the total flux from the star+disk in each time bin, and then subtract the median MCMC disk flux of 4.80 +/- 0.17 mJy from that value. To get the disk+star flux, should I just use the shortest baseline in the data? Should I extrapolate to a baseline of zero instead? If so, how do I so so in a reproducible way?

### Gas

Cail: These are probably all questions for Zach but I’d rather not bother him now that he’s free from school..  If I don’t hear from you, I’ll write to him and ask.  

- velocity / channel range?
  Cail: looks to be about -10 to 10 km/s in V_LRSK from the plotted spectrum, but I couldn’t find an explicit statement in the main text..

- spatial area?
  Cail: the caption of Zach’s AU Mic spectrum figure says that it was made with an 8” x 8” box, but I couldn’t find anything with regards to the flux upper limit itself.

- cospatial with dust?
  “For a given excitation temperature, an upper limit on the total gas mass can be inferred from the upper limit on integrated ﬂux. We ﬁnd the upper limit on the total gas mass to range between xxx and xxx for excitation temperatures between 10 K and 250 K.”  
    
  Cail: For the calculation of the gas mass upper limit, Kevin wonders if the gas was assumed to be cospatial with dust. My guess is that no such assumption was made. I guess the range of excitation temperatures considered implies some range of distances from the star? Unless I’m missing something (quite possible) I would propose to leave things as is.

## Analysis

### model formalism

- scale height definition
  “The disk vertical structure is set by the aspect ratio h = H(r)/r. At a given radius r, the vertical density proﬁle is assumed to be Gaussian with a standard deviation equal to the scale height H(r).”  
    
  MW: Can you call it “a CONSTANT aspect ratio”? As is I thought that H(r) is an arbitrary function of r whereas in fact it is linear because h is constant.   
    
  SA: I was a bit confused about the formulation of the model, particularly with respect to the vertical distribution.  You define h = H(r) / r, but maybe it would pay to explicitly state that h is independent of r.  Presuming that's correct, why would h be independent of r?  Maybe there's a theoretical reason, but its worth describing that.  
    
  Cail: This seems to have confused a whole bunch of people. I think it would be much clearer to say something like:  
    
  “At a given radius r, the vertical density proﬁle is assumed to be Gaussian with a standard deviation equal to the scale height H(r). The scale height is given by H(r) = h*r, where the aspect ratio h is a constant.” I know you said this is backwards, but it’s so much less confusing.. Also, I have no idea how to theoretically justify our choice of a linear scale height.   
    
  Unless you suggest otherwise, I’ll introduce the scale height/aspect ratio as a proposed above, and just say that the assuming h to be a linear function of r is standard for debris disk models.

- Justification for Gaussian vertical structure
  DW:  
  A basic assumption that underlies the analysis is that a Gaussian profile in the vertical dimension describes the collection of dust orbits. This is not so obvious. Some justification from dynamics would be comforting. I know for our Kuiper Belt that the vertical distribution can be approximated very well by the sum of two Gaussians, see Brown 2001, AJ, 121, 2804.  
    
  SA:  
  Speaking of which, is there a theoretical reason that the vertical distribution of particles would have a Gaussian density profile?  If so, a description of that would be helpful.  
    
  Cail: This definitely seems like something worth justifying. I’ll read/cite the paper David suggested and see if I can find a citation for extrasolar debris disks from there, but if you can think of any citations it would be appreciated!

- definition of disk center
  SA:  
  How did you decide on the disk center?  Why not fit the center position too?  
    
  Cail: Unless you suggest otherwise: I’ll add a sentence explaining that the disk center is set by the star, and that it wasn’t left as a free parameter because we already fit for the stellar location.

### fitting / MCMC

- Ray tracing and optical depth?
  “Previous studies of the scale height of debris disks have demonstrated a degeneracy between vertical structure, radial structure, and viewing geometry (e.g., Milli et al. 2014). In light of this, we adopt a modeling approach that combines appropriate ray tracing methods with a Markov Chain Monte Carlo (MCMC) algorithm in order to investigate the degree to which these known degeneracies impact our ability to measure the vertical structure of the disk.”  
    
  MW: Also, a small point is that you imply in 4 that the degeneracy leads you to adopt appropriate ray tracing methods. But here we are looking at optically thin thermal emission for which sophisticated ray tracing methods are not needed.  
    
  Cail: This seems like a place where we could be a little more precise, but I won’t do anything unless you suggest something (I’m over my head)!

- how to explain prior bounds
  “We assume uniform priors for all parameters; the dust mass was sampled in logarithmic space, formally equivalent to assuming a log uniform prior. No prior bounds were placed on the logarithm of the dust mass. Priors placed on the stellar ﬂux and disk inner radius, width, and aspect ratio ensured that these parameters were always greater than zero.”  
    
  SA:  
  I hope you wouldn't get this, but I hear rumors of ApJ editors that are nitpicking on priors.  What you've done is fine (one could argue that choosing non-informative priors is not really Bayesian, but whatever), but there is a formally incorrect statement: you can't have a uniform prior with no bounds (its mathematically impossible), so just say you have some safely large bounds.  
    
  Cail:  
  This is tricky, because when I say ‘unbounded’ I mean that the any value between -inf and inf was allowed. Unless I hear otherwise, I’ll say something like ‘Bounds placed on <whatever parameter> were chosen to be large enough so as not to impose any prior assumptions” but these feels just a tiny bit disingenuous.. maybe we should be explicit and say that we left the bounds at +/- infinity, but that these should be thought of as ‘safely large bounds’ rather than unbounded.

- How to better discuss the inclination priors?
  %EIC: "Leaving the inclination unbounded does not affect the preferred inclination ..."  
  %      This sounds funny and confusing for a couple reasons.  
  %      -- Leaving the inclination unbounded would seem to be the best possible thing  
  %      to do because "unbounded" implies "unbiased"; we don't want to affect the outcome in some unfair way.  
  %      -- Technically the inclination is bounded --- it is confined to be between 0 and 90 deg.  
  %      I think this needs to be reworded.  
  Leaving the inclination unbounded does not affect the preferred inclination; the MCMC chain simply converges to two positions symmetric about \ang{90}.  
    
  SA:  
  On p11, I was confused by the statement about the inclination that "MCMC chain simply converges to two positions symmetric about 90 degrees."  What does that mean?  
    
  Cail: I wouldn’t even ask about this, but I’ve really been struggling to find a clear way of expressing this. I’d propose replacing the offending sentence with the following, but please let me know if you can think of any improvements!  
    
  “Because AU Mic is so close to edge-on, the preferred inclination falls very close to the 90º prior upper bound. To ensure that the proximity of the solution to the edge of parameter space does not affect the posterior distribution, we investigated the effect of allowing the inclination to vary above 90º. Doing so produced a inclination distribution with two symmetric modes on either side of 90º. When the resulting ’unbounded’ inclination distribution was reparameterized such that all values fell between 0º and 90º, the original ‘bounded’ inclination distribution was recovered, indicating that the placing a 90º prior upper bound has no effect on the inclination posterior.”

- chi^2 value too high?
  “This model formalism resulted in a best-ﬁt χ 2 value of 626163.598 (reduced χ 2 = 2.053)”  
    
  JC: This is pretty high value and formally the model is rejected and the uncertainties should be suspect. Presumably this is to the residuals around the star that you show later.  
    
  Cail: I don’t think there’s anything in particular we can do about this, right? I guess we could add a sentence explaining that the high chi^2 value is due to deviations of the surface density profile from a simple power law, but I won’t do so unless you say so.

- add definition of lnprob?
  DW: Table 3, if the Likelihood is going to reported here, then I think its definition should be in the text in S4.1.  
    
  Cail: I’m happy to add the definition to the text, but I wonder if we shouldn’t just remove the ‘ln Likelihood’ row from the table… The only purpose it really serves is to provide a point of comparison between the two model formalisms listed in the table, and there are of course better metrics.. I’ll take it out unless i hear otherwise. 

- corner plot labeling
  %EIC: "Leaving the inclination unbounded does not affect the preferred inclination ..."  
  %      This sounds funny and confusing for a couple reasons.  
  %      -- Leaving the inclination unbounded would seem to be the best possible thing  
  %      to do because "unbounded" implies "unbiased"; we don't want to affect the outcome in some unfair way.  
  %      -- Technically the inclination is bounded --- it is confined to be between 0 and 90 deg.  
  %      I think this needs to be reworded.  
  Leaving the inclination unbounded does not affect the preferred inclination; the MCMC chain simply converges to two positions symmetric about \ang{90}.

## Discussion

### compare flux density instead of dust mass?

DW:  
The "dust mass" here traces back to the flux density and assumed opacity, so maybe it is really the flux density measurements that should be compared?.  
  
Cail:  
My only would objection would be the impact of observing wavelength? I guess that we could scale all the fluxes to the same wavelength if we assume blackbody emission, but this seems like a mess. I would propose to leave things as is.

### radial structure

- log-log SB profile?
  DW:  
  Figure 6, since all the fits are power-laws, would it make more sense for this log-linear plot into a log-log plot?  
    
  Cail: I don’t really like this idea, especially as a logarithmic x-axis would smoosh together the separations we care most about (~40 au). I’d propose to do nothing.

- constrain ‘bump’ separation with SED?
  MW: One comment for this section is to pose the question of whether we can use the SED to distinguish between the bump at 10au projection being caused by emission at 10au or 40au? That is, would there be too much mid-IR emission from a ring at 10au at that level?  
    
  Cail: This is a cool idea, but seems beyond the scope/purpose of this paper. I’d propose to do nothing.

### vertical structure

- do these sentences sound out of place?
  “Because the measurements quoted above are determined from the observed vertical thickness of AU Mic’s disk, they can be aﬀected by the radial structure and viewing geometry of the disk as well as scattering eﬀects. As such, parametric modeling provides a more reliable way to assess AU Mic’s vertical structure. Krist et al. (2005) report a FWHM…”  
    
  Kevin:  
  The first two sentences seem a little out of place here. It sounds like you are setting up your parametric modeling, but then you go on to discuss other results from the literature. Are you referring to your own parametric modeling? Or are you referring to modeling done by Krist et al.?  
    
  Cail: I think these sentences sound fine and do a good job of distinguishing between apparent and modeled FWHMs. As such, I’d propose to do nothing. However, if Kevin was confused, maybe other people will too?

- compare to protoplanetary scale heights?
  SA: It might be worthwhile to make a few flippant comments comparing the mm scale height for AU Mic to primordial disks.  For example, HL Tau seems to have an even lower aspect ratio (h ~ 0.01; see Pinte et al. 2016).  There's a few edge-on cases with similar, ~few percent, aspect ratios.  I certainly recognize the physical differences in these scenarios, but its interesting to think about this in an evolutionary context.  
    
  Cail: This is a cool idea! I’ll start with the HL Tau reference, but can you think of any other citations/disks?

### uncertainty propagation with stellar mass?

Cail:  
  
DW wonders about propagating uncertainties on the stellar mass. This doesn’t seem very feasible, considering that the mass we use is an average of several values reported in the literature (see below). The only thing I can think of is estimating our own uncertainty, which seems bad.. I’d propose leaving things as is.  
  
From Schuppler, where I got 0.5 M_sun:  
“We assumed a stellar mass of 0.5 M_sun, motivated by the mean of the wide mass range given in the literature (0.3-0.6 M_sun , Plavchan et al. 2009; Houdebine & Doyle 1994).”  
  
this is hard because we’re cobbling together an estimate from several sources…

### 1-D version of previous plot?
![](assets/BB046249-B646-499A-BFD4-3F618C5638D6.png)

## Conclusions

### Plavchan’s Planet

- things aren’t so obvious to Eugene
  %EIC: begin  
  %EIC: maybe "obvious" to you, but not obvious to me. How is a ring at 10 AU "explained" by a planet?  
  %       I deleted the speculation, but left intact the reference to Peter Plavchan.  
  %If a secondary ring of dust does exist in the AU Mic debris disk, planets are of course an obvious explanation;   
  P.~Plavchan et al.~(in preparation) propose a Jovian-mass exoplanet candidate interior to \SI{1}{au}, but AU Mic's stellar activity make it difficult to confirm the radial velocity detection.  
  %Even if the candidate is confirmed, it is unlikely that a planet so close to the star would be responsible for a ring at \SI{10}{au}.  
  %EIC: Sorry, I disagree strongly with the last sentence. Perhaps you can explain to me why a confirmed  
  %       detection of a planet (particularly one interior to 1 AU) would explain AU Mic's fast-moving features?  
  %       I am happy to discuss why I do not think such an explanation is viable, if you would like (skype/phone/email).  
  %That being said, a confirmed detection could provide a promising origin for AU Mic's fast-moving features.  
  %EIC: end  
    
  Cail: This one made me chuckle. To Eugene’s first objection, I thought it was common knowledge that planets can sculpt rings through pressure ‘bumps,’ etc. I guess I should be a little more conservative and say that planets are often invoked to explain rings, rather than saying that a planet provides an obvious explanation. I should also probably include a citation.. maybe your debris disk review paper, or Mark Wyatt’s 1999 thesis for a more theoretical approach. Can you think of a more relevant paper?  
    
  To Eugene’s second point: In Plavchan’s email to me about the planet, he said: “Such a planet would be an excellent explanation for the moving structures observed further out in the disk with direct imaging. While that's a very compelling story to tell, I haven't convinced myself yet that the planet is genuine.” Going back to look at the Sezestre (2017), the best-fit planet separation is ~8 au, but separations <= 5 au are still acceptable. Thus it doesn’t seem so crazy to link Plavchan’s planet to the moving features. I suspect Eugene is disagreeing so strongly because he doesn’t believe the conclusions of Sezestre?  
    
  I’d propose rewriting the sentences Eugene doesn’t like in a more conservative manner (as outlined above), carefully outlining the logic and doing a better job of citing things. I’m not particularly interested in getting involved in an email discussion with Eugene, so ideally Eugene can just add his comments to the next draft.. but maybe I should email him about it directly, since he seems so opinionated? I won’t do this last part unless you think it’s necessary. 

- move to discussion of ring?
  Cail: Someone (I think Kevin?) suggested moving the discussion of Plavchan’s planet to, well, the discussion. Unless you suggest otherwise, this is what I’ll do!

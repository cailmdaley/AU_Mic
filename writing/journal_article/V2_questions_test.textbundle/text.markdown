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

“Thus longer-wavelength ($\lambda \geq 50$ µm) observations are required to measure disk scale height as determined by dynamical stirring alone, since the grains dominating the emission at these wavelengths are large enough to be impervious to the eﬀects of radiation pressure”  
  
**JC:** This must depend on the spectral type of the host star. You might want to give a bit more background on how the limit varies with SpT.  
  
**Cail:** I couldn’t find any dependence on stellar properties in the Thebault paper. I think this is because particle sizes are simply described by their beta values (see below). These beta values of course carry information about the star's radiation pressure (and probably spectral type), but where Thebault gets these values from isn’t clear… I would propose to leave things as is.

**Thebault:** "As for the orbital evolution computation, the important parameter is the value of $\beta$ (s) for estimating the radiation pressure force. We thus chose to scale all particle sizes by their $\beta(s) \propto$ 1/s value… We consider here s min and s max , such that $\beta$ (s min ) = 0.4 (close to the blow-out cut-oﬀ size) and $\beta$ (s max ) = 0.025. This size range covers however most of the crucial grain population that dominates the ﬂux in scattered light." 
  
 *Note [page 1]:* John is right and you are right.  Thebault gives values in terms of beta, and you give a value in terms of a wavelength, which is effectively in terms of a grain size — so you are indeed implicitly assuming a spectral type.  I think everyone will be happy if you just say how you settled on 50um (I mean, you must have assumed a stellar luminosity or something somewhere along the way, right?).  That would involve you mentioning the beta constraints from the Thebault paper that you quote here.  


## Observations

### Flare

**Cail:** The larger question here is: do we want to give more attention to the flare in this paper, or save it for somewhere/someone else? There’s clearly much more analysis we could do w/r/t the flare, but it would definitely take a little time. I would propose plotting the flare and looking at the polarization vs time (both should be very quick).  If something interesting pops up in the polarization, it might be worth expanding the flare discussion (maybe giving it its own subsection?)

 *Note [page 1]:* I don’t think it’s really worth plotting the flare (the table is more useful to anyone who wants to work on the flare data later, and a plot with no interpretation won’t add anything to the paper).  But sure, check the polarization vs. time if it’s quick — that would give a new dimension to the data that we don’t already include and which we already know is likely to be highly relevant.  

- plot flare light curve?    

- ‘deblend’ with uniform weighting?  
  “Faithfully disentangling the two components proved diﬃcult, both because stellar and disk emission are blurred together in the central beam and because the stellar ﬂux varied signiﬁcantly across the three nights of observation (Table 3)”  
    
  **JC:** Can you deblend them a bit by using uniform weighting?  
    
  **Cail:** I don’t really know what this means, but I suspect that it won’t be better than fitting a point source directly to the visibilities. I would leave things as-is.  
  
 *Note [page 2]:* I agree.  Uniform weighting won’t meaningfully change things.  The problem is that the disk is edge-on, so all the resolution in the world won’t allow you to disentangle the disk from the star in the central beam. You could probably make a slightly more clear statement to that effect in the relevant section of the paper.  

- look at polarization  
  **DW:** Regarding the flare, did you look at the polarization? In the Prox Cen data, there was a clear signature of polarization in a difference in XX and YY vs. time. It's not worth getting distracted by the flare, and maybe this analysis doesn't belong at all in this paper, or maybe it could be in an appendix if you have done the fits for XX and YY already.

### Miriad

**SA:** You claim to use MIRIAD (p4) in the calibration, but it doesn't seem like you do until you measure some emission.  My hope is it was only used for that!  Even so, why not just use the CASA viewer to do the same thing?  
  
**Cail:** Sean is right, I don’t think we ended up using MIRIAD at all in the data reduction. However, I would propose not re-measuirng the disk flux in CASA—I know you love cgcurs!
 *Note [page 2]:* I don’t love cgcurs, but I don’t actually know how to do the same thing in CASA off the top of my head, and of course it doesn’t matter.  So just make sure the paper is clear that you only use MIRIAD to measure the integrated flux.  

### self-calibration?

**SA:** The calibration discussion makes no mention of self-calibration.  Was that tried?  If not, why not?  If so, and it didn't help, it would be good to say so and report the interval(s) attempted.  [With a peak SNR ~ 20, I would think you could get a boost from self-cal.]  
  
**Cail:** I got nothing here. I don’t really know what self-calibration is, let alone why we didn’t do it. I guess I’ll ignore this comment unless you suggest otherwise..

 *Note [page 2]:* Oh, Sean is right that we should probably have tried self-calibration — it’s not something I’ve really acclimated to yet because usually debris disks are too faint for it to be useful (AU Mic is an exception). But, at this point I don’t think it’s worth going back and redoing everything unless the referee asks for it.  It could give us an SNR boost.  I would probably ignore it, but I’m not proud of ignoring it.  
(The basic idea behind self-calibration is that after doing the gain cal with a quasar, you do another set of gain cals that uses the signal from the source itself, assuming it’s bright enough to detect in a single loop — you do need a bright source to do it, but if you can do it, it’s obviously an improvement over using a quasar that was observed asynchronously with the source and which has its own peculiar variations due to its slightly different position in the sky.)

## Results

### AU Mic images

**DW:** Figure 1, consider to separate the panels in order to make it easier for readers to lift the left panel only from the paper, to show in talks  
  
**Cail:** This comment made me realize that the tapered clean image has little-to-no function. Maybe we should take it out? I’ll only add some spacing between the panes, unless you think I should take out the tapered clean altogether.
 *Note [page 2]:* I’d be on board with taking out the tapered clean image.  As you say, it has basically no function, and as David says, people really just want to lift the “best” image to show in talks so it’s worthwhile to make your money plot as simple as possible. 

### Boccaletti Plots

- Did Eugene misunderstand?  
  "From top to bottom, the three panes show i) the disk spine surface brightness with errors given by the rms noise, ii) the disk spine deviation from the midplane, and iii) the disk FWHM as a function of projected separation from the star.   
  The Gaussian traced by the dotted line in the first pane shows the projected width in the radial direction of the naturally weighted synthesized beam of the combined dataset.  
  %**EIC:** Following sentence ("The fact that ...") is not obvious to me just by looking at the figure.  
  %       In fact, the beam size as indicated by the Gaussian drawn on the top panel looks like it has a FWHM  
  %       of ~5 AU > the beam-subtracted FWHM plotted in the third panel. (I was pleased to see this  
  %       acknowledged in the main text, and to read about the various reasons why this is not necessarily  
  %       cause for alarm.)  
  %       I think the last sentence just needs to be reworded --- to acknowledge that the vertical beam  
  %       is technically larger than what is plotted, but that this is OK for various reasons discussed in  
  %       the main text under section 4.3. I gave it a stab.  
  %  The fact that the image-domain vertical height of the disk is in excess of the beam contribution implies that our data spatially resolve the vertical structure of the disk.  
  Although the last panel indicates a vertical disk FWHM that is technically smaller than the vertical FWHM of the beam, it still represents a significant detection of vertical structure, for reasons explained in \S \ref{subsection: vertical analysis}.}"
    
  **Cail:** This one’s tricky. Right off the bat, it seems like Eugene took the size of the beam in the radial direction to be the size in the vertical direction. That being said, the radial & vertical FWHMs aren’t that different so I don’t think this disqualifies his comments. It seems that he’s concerned by the fact that the beam-subtracted disk FWHM is smaller than the beam FWHM itself (this does strike me as odd, but I trust the math).  
    
  I would propose explicitly stating that the beam radial FWHM plotted in the top pane is *not* the beam vertical FWHM. Also comparing the beam-subtracted FWHM to the data resolution as Eugene suggested, and explaining that before subtracting the beam FWHM, the Gaussian FWHM was 8 au (or whatever it was exactly). However, I’m not really sure how to handle the sentence Eugene commented out and the one he replaced it with. On the one hand, the beam-subtracted disk FWHM is smaller than the beam, which might be grounds for calling it ‘unresolved’? On the other hand, the apparent disk FWHM (~8 au) exceeds that of the beam (~ 5 au), indicating that the disk is resolved. But wouldn’t this be true when convolving the beam with *any* Gaussian vertical profile with non-zero standard deviation? 
 *Note [page 3]:* Yes, Eugene misunderstood, and your sentence is correct (his is not).  I think the best way to clarify would be to add to the previous sentence: “iii) the disk FWHM as a function of projected separation from the star[, after subtraction of the beam FWHM in the vertical direction.]”  And yes, make it explicit that the beam radial FWHM differs from the vertical FWHM.  

- uncertainty / fit information  
  **MW:** What I am missing from this is a feeling for the typical uncertainties in these profiles, as well as how good a fit a Gaussian profile is to these vertical profiles.  
    
  **Cail:** I used astropy’s gaussian fitter, which doesn’t report uncertainties, to make these plots. It probably wouldn’t take too long to find & implement a fitter that does, and it would be really nice to have uncertainties on all three panes (this would also probably give a better estimate of the uncertainty for the top panel than just using the rms noise). To Mark’s second comment about whether the vertical profile is actually Gaussian, I could add another inlay to the top pane showing a sample vertical slice (or maybe the mean profile from all the slices?).  
    
  In the interest of getting the second draft out quickly, I would propose only adding a note stating that the vertical profile is well described by a Gaussian. But let me know if you’d like me to do either of the things suggested above—I think they’re good ideas.
 *Note [page 4]:* I think your idea about adding another inlay to the top pane showing a sample vertical profile is a nice one to give a sense of how representative the Gaussian profile is.  The mean profile from all slices is also a nice idea, but runs into the problem that the slices aren’t independent because of the beam.  Maybe a median instead?  I’d have to think about this a bit more, but I think a median is probably OK.   Either way, keep this fairly low priority — adding a little explanation is also fine.  

- rotate by MCMC PA?  
  **Cail:** DW was wondering how 0 elevation/deviation from the midplane was defined in the second pane. As I recall, I just did this by eye in the past. This is a little anachronistic, but I would propose using the MCMC best-fit PA (and explaining that we did so).
 *Note [page 4]:* Sure, that sounds like a good plan.

### Disk symmetry and eccentricity

**HS:** Also ion page 7, I suggest adding one sentence between ‘millimeter wavelength.’ and ‘Additional’ that explains what the absence of any significant brightness difference between the two limbs implies - and maybe even cite Margaret's paper on this. - Margaret, do you want to send one or two sentence about this?  
  
**MP:** Hilke mentioned the lack of significant brightness differences between the ansae and implications for the disk eccentricity (beginning of section 3). You have some info on the disk eccentricity from that but it depends on the disk's argument of periastron (if the argument of periastron is in the sky plane, the eccentricity would be the ratio of the fluxes minus unity, or 0.03 +/- 0.06, so consistent with zero). I could make a plot of the implied eccentricity as a function of the argument of periastron if you want to include a discussion of that - it would be something like "there is a 68% (or 95%, or 99%) chance the eccentricity is less than <a certain value>".  
  
**Cail:** This all seems like a great idea, especially as it doesn’t require any work on our end! Unless I hear otherwise, I’ll email Margaret asking her to make this plot.
 *Note [page 5]:* Go for it!  Sounds fun.  Also a bit of a headache because of the edge-on geometry, but if Margaret is up for tackling it, great!


### Disk flux discepancy

**JC:** From Table 1, the star(without the flare)  seems to have a flux of 0.45 mJy. However, you report a total disk+star flux of 4. 97 mJy, and a disk-only flux of 4.8 mJy. Thus I would infer a stellar flux of ~ 0.17 mJy. Is the apparent discrepancy because the 3 sigma contour is not enclosing all of the flux?  
  
**Cail:** This one has me a little worried. The pre-flare flux of 0.45 mJy presented in the table is very much not consistent with the MCMC June flux of 0.220 +/- 0.020 mJy. My guess is that we overestimated the pre-flare stellar emission when fitting a point source to the visibilities.   
  
Possibly a better way to do this would be to get the total flux from the star+disk in each time bin, and then subtract the median MCMC disk flux of 4.80 +/- 0.17 mJy from that value. To get the disk+star flux, should I just use the shortest baseline in the data? Should I extrapolate to a baseline of zero instead? If so, how do I so so in a reproducible way?
 *Note [page 5]:* Hm, yes, easy to overestimate the stellar flux, since it’s so hard to disentangle the disk and star without the MCMC analysis.  I generally like your idea of using the total flux and then subtracting the disk flux from the MCMC fit.  However, one worry is that John is right that the 3-sigma contour may not actually enclose all of the flux, in which case you’d be underestimating the stellar flux by doing what you suggest.  Another worry is that it’s going to be basically impossible to get the star+disk flux from a simple fit to the visibilities given the complex flux distribution of the source.  Neither method you suggest is really going to be reliable.  I wonder if you can just do an image-domain fit to the central point source for each time bin, and then subtract the disk flux in that central beam from the MCMC fit — do you think that might work?

### Gas

**Cail:** These are probably all questions for Zach but I’d rather not bother him now that he’s free from school.. But if you don't know off the top of your head, I’ll write to him and ask!

- velocity / channel range?  
  **Cail:** looks to be about -10 to 10 km/s in V_LRSK from the plotted spectrum, but I couldn’t find an explicit statement in the main text..
 *Note [page 5]:* I’m pretty sure that’s right. 

- spatial area?  
  **Cail:** the caption of Zach’s AU Mic spectrum figure says that it was made with an 8” x 8” box, but I couldn’t find anything with regards to the flux upper limit itself.
 *Note [page 5]:* Yes, the area within which the spectrum was integrated was 8x8”, so I think that’s the spatial area you mean (though I’m not sure exactly what you’re asking — I’d probably have to figure out exactly where in your draft you’re talking about)

- cospatial with dust?  
  “For a given excitation temperature, an upper limit on the total gas mass can be inferred from the upper limit on integrated ﬂux. We ﬁnd the upper limit on the total gas mass to range between xxx and xxx for excitation temperatures between 10 K and 250 K.”  
    
  **Cail:** For the calculation of the gas mass upper limit, Kevin wonders if the gas was assumed to be cospatial with dust. My guess is that no such assumption was made. I guess the range of excitation temperatures considered implies some range of distances from the star? Unless I’m missing something (quite possible) I would propose to leave things as is.
 *Note [page 6]:* Yes, that’s fine.  The gaswas only assumed to be within the 8x8” box. 

## Analysis

### model formalism

- scale height definition  
  “The disk vertical structure is set by the aspect ratio h = H(r)/r. At a given radius r, the vertical density proﬁle is assumed to be Gaussian with a standard deviation equal to the scale height H(r).”  
    
  **MW:** Can you call it “a CONSTANT aspect ratio”? As is I thought that H(r) is an arbitrary function of r whereas in fact it is linear because h is constant.   
    
  **SA:** I was a bit confused about the formulation of the model, particularly with respect to the vertical distribution.  You define h = H(r) / r, but maybe it would pay to explicitly state that h is independent of r.  Presuming that's correct, why would h be independent of r?  Maybe there's a theoretical reason, but its worth describing that.  
    
  **Cail:** This seems to have confused a whole bunch of people. I think it would be much clearer to say something like:  
    
  “At a given radius r, the vertical density proﬁle is assumed to be Gaussian with a standard deviation equal to the scale height H(r). The scale height is given by H(r) = h*r, where the aspect ratio h is a constant.” I know you said this is backwards, but it’s so much less confusing.. Also, I have no idea how to theoretically justify our choice of a linear scale height.   
    
  Unless you suggest otherwise, I’ll introduce the scale height/aspect ratio as a proposed above, and just say that the assuming h to be a linear function of r is standard for debris disk models.
 *Note [page 6]:* OK, fine.  And there’s no justification that I know of other than that it’s a common assumption and that we’re not justified in assuming anything more complicated given the resolution of our data compared with the size of the disk (I think we tested this at some point with the constant scale height assumption, right?)  

- Justification for Gaussian vertical structure  
  **DW:** A basic assumption that underlies the analysis is that a Gaussian profile in the vertical dimension describes the collection of dust orbits. This is not so obvious. Some justification from dynamics would be comforting. I know for our Kuiper Belt that the vertical distribution can be approximated very well by the sum of two Gaussians, see Brown 2001, AJ, 121, 2804.  
    
  **SA:** Speaking of which, is there a theoretical reason that the vertical distribution of particles would have a Gaussian density profile?  If so, a description of that would be helpful.  
    
  **Cail:** This definitely seems like something worth justifying. I’ll read/cite the paper David suggested and see if I can find a citation for extrasolar debris disks from there, but if you can think of any citations it would be appreciated!
 *Note [page 7]:* Yeah, David really helped you out there.  I’ve got nothing to add of the top of my head.  

- definition of disk center  
  **SA:** How did you decide on the disk center?  Why not fit the center position too?  
    
  **Cail:** Unless you suggest otherwise: I’ll add a sentence explaining that the disk center is set by the star, and that it wasn’t left as a free parameter because we already fit for the stellar location.
 *Note [page 7]:* Yup, sounds good.

### fitting / MCMC

- Ray tracing and optical depth?  
  “Previous studies of the scale height of debris disks have demonstrated a degeneracy between vertical structure, radial structure, and viewing geometry (e.g., Milli et al. 2014). In light of this, we adopt a modeling approach that combines appropriate ray tracing methods with a Markov Chain Monte Carlo (MCMC) algorithm in order to investigate the degree to which these known degeneracies impact our ability to measure the vertical structure of the disk.”  
    
  **MW:** Also, a small point is that you imply in 4 that the degeneracy leads you to adopt appropriate ray tracing methods. But here we are looking at optically thin thermal emission for which sophisticated ray tracing methods are not needed.  
    
  **Cail:** This seems like a place where we could be a little more precise, but I won’t do anything unless you suggest something (I’m over my head)!
 *Note [page 7]:* Yes, I see where Mark is confused (it’s his confusion, not yours).   I think you can just change your sentence to say “In light of this, we adopt a modeling approach that combines a simple ray-tracing code to properly project the radial and vertical flux distribution of the optically thin emission onto the sky plane with an MCMC fitting algorithm that allows us to explore the degree to which these known degeneracies…”  There are indeed radiative transfer codes that use MCMC methods just to do the sky-projection part, but that’s not what you’re doing, so I just tried to make your language more explicit about the purpose/scope of the ray-tracing vs. the MCMC.  

- how to explain prior bounds  
  “We assume uniform priors for all parameters; the dust mass was sampled in logarithmic space, formally equivalent to assuming a log uniform prior. No prior bounds were placed on the logarithm of the dust mass. Priors placed on the stellar ﬂux and disk inner radius, width, and aspect ratio ensured that these parameters were always greater than zero.”  
    
  **SA:** I hope you wouldn't get this, but I hear rumors of ApJ editors that are nitpicking on priors.  What you've done is fine (one could argue that choosing non-informative priors is not really Bayesian, but whatever), but there is a formally incorrect statement: you can't have a uniform prior with no bounds (its mathematically impossible), so just say you have some safely large bounds.  
    
  **Cail:** This is tricky, because when I say ‘unbounded’ I mean that the any value between -inf and inf was allowed. Unless I hear otherwise, I’ll say something like ‘Bounds placed on <whatever parameter> were chosen to be large enough so as not to impose any prior assumptions” but these feels just a tiny bit disingenuous.. maybe we should be explicit and say that we left the bounds at +/- infinity, but that these should be thought of as ‘safely large bounds’ rather than unbounded.
 *Note [page 8]:* I think the language you suggest is fine.  We can deal with it if the editor asks.  

- How to better discuss the inclination priors?  
  %**EIC:** "Leaving the inclination unbounded does not affect the preferred inclination ..."  
  %      This sounds funny and confusing for a couple reasons.  
  %      -- Leaving the inclination unbounded would seem to be the best possible thing  
  %      to do because "unbounded" implies "unbiased"; we don't want to affect the outcome in some unfair way.  
  %      -- Technically the inclination is bounded --- it is confined to be between 0 and 90 deg.  
  %      I think this needs to be reworded.  
  Leaving the inclination unbounded does not affect the preferred inclination; the MCMC chain simply converges to two positions symmetric about ang{90}.  
    
  **SA:** On p11, I was confused by the statement about the inclination that "MCMC chain simply converges to two positions symmetric about 90 degrees."  What does that mean?  
    
  **Cail:** I wouldn’t even ask about this, but I’ve really been struggling to find a clear way of expressing this. I’d propose replacing the offending sentence with the following, but please let me know if you can think of any improvements!  
    
  “Because AU Mic is so close to edge-on, the preferred inclination falls very close to the 90º prior upper bound. To ensure that the proximity of the solution to the edge of parameter space does not affect the posterior distribution, we investigated the effect of allowing the inclination to vary above 90º. Doing so produced a inclination distribution with two symmetric modes on either side of 90º. When the resulting ’unbounded’ inclination distribution was reparameterized such that all values fell between 0º and 90º, the original ‘bounded’ inclination distribution was recovered, indicating that the placing a 90º prior upper bound has no effect on the inclination posterior.”
 *Note [page 8]:* Awesome.  I think this lays it out really clearly.  

- chi^2 value too high?  
  “This model formalism resulted in a best-ﬁt chi 2 value of 626163.598 (reduced chi 2 = 2.053)”  
    
  **JC:** This is pretty high value and formally the model is rejected and the uncertainties should be suspect. Presumably this is to the residuals around the star that you show later.  
    
  **Cail:** I don’t think there’s anything in particular we can do about this, right? I guess we could add a sentence explaining that the high chi^2 value is due to deviations of the surface density profile from a simple power law, but I won’t do so unless you say so.
 *Note [page 9]:* Well, John is right, and I think it’s worth saying something (but being nonspecific, because we really don’t know).  You can say something like, “The reduced chi^2 value of 2.053 indicates that our best-fit model is not a perfect fit to the data, perhaps due to the 3-sigma residuals discussed in Section [whatever].  The statistical weights are likely not responsible, due to our recalculation of the weights according to the observed variance within each grid point in the u-v plane; these were checked to ensure that the distribution was indeed Gaussian.  

- add definition of lnprob?  
  **DW:** Table 3, if the Likelihood is going to reported here, then I think its definition should be in the text in S4.1.  
    
  **Cail:** I’m happy to add the definition to the text, but I wonder if we shouldn’t just remove the ‘ln Likelihood’ row from the table… The only purpose it really serves is to provide a point of comparison between the two model formalisms listed in the table, and there are of course better metrics.. I’ll take it out unless i hear otherwise. 
 *Note [page 9]:* I think it’s probably better to leave it in and add the definition as David suggests.  

- corner plot labeling  
  %**EIC:** "Leaving the inclination unbounded does not affect the preferred inclination ..."  
  %      This sounds funny and confusing for a couple reasons.  
  %      -- Leaving the inclination unbounded would seem to be the best possible thing  
  %      to do because "unbounded" implies "unbiased"; we don't want to affect the outcome in some unfair way.  
  %      -- Technically the inclination is bounded --- it is confined to be between 0 and 90 deg.  
  %      I think this needs to be reworded.  
  Leaving the inclination unbounded does not affect the preferred inclination; the MCMC chain simply converges to two positions symmetric about ang{90}.
 *Note [page 9]:* Um, didn’t we address this above?  or did you mean to put a different quote here, since this bullet point is labeled “corner plot labeling” and the Eugene quote given doesn’t actually seem to refer to corner plot labeling?

## Discussion

### compare flux density instead of dust mass?

**DW:** The "dust mass" here traces back to the flux density and assumed opacity, so maybe it is really the flux density measurements that should be compared?.  
  
**Cail:** My only would objection would be the impact of observing wavelength? I guess that we could scale all the fluxes to the same wavelength if we assume blackbody emission, but this seems like a mess. I would propose to leave things as is.
 *Note [page 9]:* Ugh, I mostly agree with you.  You’ve had too much experience in how hard it is to compare flux values!  That said, the flux is the real observable, and you’ve also had plenty of experience with how annoying it is to try to take a quoted dust mass and figure out what the original flux was by having to back out what opacity and temperature was assumed.  So, just make sure that both are given in the paper somewhere?

### radial structure

- log-log SB profile?  
  **DW:** Figure 6, since all the fits are power-laws, would it make more sense for this log-linear plot into a log-log plot?  
 *Note [page 10]:* Yup, fine with me.  It’s a tradeoff.  David is right that since the fits are power laws a log-log plot would show them most cleanly, but you are right that a log-log plot would also overemphasize the regions of the disk that we aren’t really constraining.  You’re the first author and therefore the tie-breaker!
    
  **Cail:** I don’t really like this idea, especially as a logarithmic x-axis would smoosh together the separations we care most about (~40 au). I’d propose to do nothing.

- constrain ‘bump’ separation with SED?  
  **MW:** One comment for this section is to pose the question of whether we can use the SED to distinguish between the bump at 10au projection being caused by emission at 10au or 40au? That is, would there be too much mid-IR emission from a ring at 10au at that level?  
    
  **Cail:** This is a cool idea, but seems beyond the scope/purpose of this paper. I’d propose to do nothing.
 *Note [page 10]:* This is an interesting idea indeed, and presumably one that Mark can explore as he publishes his student’s work.  Rather than doing nothing, though, I’d suggest that you add a statement to this effect as sort of a hypothetical, to show that you / your coauthors have discussed the point.  Just mention that in principle the SED might be able to help disentangle the radius of the bump at 10au projected separation, but that using the temperature to infer the radius would be subject to assumptions about the dust grain size distribution (I think that’s true, although it’s late so I’m not sure — double-check me!).  

### vertical structure

- do these sentences sound out of place?  
  “Because the measurements quoted above are determined from the observed vertical thickness of AU Mic’s disk, they can be aﬀected by the radial structure and viewing geometry of the disk as well as scattering eﬀects. As such, parametric modeling provides a more reliable way to assess AU Mic’s vertical structure. Krist et al. (2005) report a FWHM…”  
    
  **Kevin:** The first two sentences seem a little out of place here. It sounds like you are setting up your parametric modeling, but then you go on to discuss other results from the literature. Are you referring to your own parametric modeling? Or are you referring to modeling done by Krist et al.?  
    
  **Cail:** I think these sentences sound fine and do a good job of distinguishing between apparent and modeled FWHMs. As such, I’d propose to do nothing. However, if Kevin was confused, maybe other people will too?
 *Note [page 10]:* I also think they’re fine, but I think I’m too close to the project to be a reliable witness.  You could probably tweak it a bit to satisfy Kevin and readers like him (he’s a sharp guy, after all!).  What about adding “…provides a more reliable way to assess AU Mic’s vertical structure, and must be taken into account in any comparison of our results with values in the literature.  For example…”

- compare to protoplanetary scale heights?  
  **SA:** It might be worthwhile to make a few flippant comments comparing the mm scale height for AU Mic to primordial disks.  For example, HL Tau seems to have an even lower aspect ratio (h ~ 0.01; see Pinte et al. 2016).  There's a few edge-on cases with similar, ~few percent, aspect ratios.  I certainly recognize the physical differences in these scenarios, but its interesting to think about this in an evolutionary context.  
    
  **Cail:** This is a cool idea! I’ll start with the HL Tau reference, but can you think of any other citations/disks?
 *Note [page 10]:* Hm, interesting.  I don’t know of any off the top of my head (Sean knows the protoplanetary literature much better than I do these days), but I can give you some search ideas.  The butterfly star might be one.  There’s a recent Guilloteau ALMA study that I think looked at some edge-on disk or other (can’t remember which).  Any ALMA papers by Gaspard Duchene or Karl Stapelfeldt would also be a good place to look, since they care a lot about vertical structure of edge-on sources.  Otherwise, I think it’s up to you to do some searching.  Don’t go wild — as Sean acknowledges, this is more a point of curiosity than a fundamental point in the interpretation of our results.  

### uncertainty propagation with stellar mass?

**Cail:** DW wonders about propagating uncertainties on the stellar mass. This doesn’t seem very feasible, considering that the mass we use is an average of several values reported in the literature (see below). The only thing I can think of is estimating our own uncertainty, which seems bad.. I’d propose leaving things as is.  
  
From Schuppler, where I got 0.5 M_sun:  
“We assumed a stellar mass of 0.5 M_sun, motivated by the mean of the wide mass range given in the literature (0.3-0.6 M_sun , Plavchan et al. 2009; Houdebine & Doyle 1994).”  
 *Note [page 11]:* Yes, I agree.  Just leave as-is.



## Conclusions

### Plavchan’s Planet

- things aren’t so obvious to Eugene  
  %**EIC:** begin  
  %**EIC:** maybe "obvious" to you, but not obvious to me. How is a ring at 10 AU "explained" by a planet?  
  %       I deleted the speculation, but left intact the reference to Peter Plavchan.  
  %If a secondary ring of dust does exist in the AU Mic debris disk, planets are of course an obvious explanation;   
  P.~Plavchan et al.~(in preparation) propose a Jovian-mass exoplanet candidate interior to SI{1}{au}, but AU Mic's stellar activity make it difficult to confirm the radial velocity detection.  
 *Note [page 11]:* This is a fair statement, I think. 
  %Even if the candidate is confirmed, it is unlikely that a planet so close to the star would be responsible for a ring at SI{10}{au}.  
  %**EIC:** Sorry, I disagree strongly with the last sentence. Perhaps you can explain to me why a confirmed  
  %       detection of a planet (particularly one interior to 1 AU) would explain AU Mic's fast-moving features?  
  %       I am happy to discuss why I do not think such an explanation is viable, if you would like (skype/phone/email).  
  %That being said, a confirmed detection could provide a promising origin for AU Mic's fast-moving features.  
  %**EIC:** end  
    
  **Cail:** This one made me chuckle. To Eugene’s first objection, I thought it was common knowledge that planets can sculpt rings through pressure ‘bumps,’ etc. I guess I should be a little more conservative and say that planets are often invoked to explain rings, rather than saying that a planet provides an obvious explanation. I should also probably include a citation.. maybe your debris disk review paper, or Mark Wyatt’s 1999 thesis for a more theoretical approach. Can you think of a more relevant paper?  
 *Note [page 11]:* Pressure bumps don’t work where there’s no gas.  It’s more of a stirring argument.  Yes, it’s pretty unobjectionable to say that planets are often invoked to explain rings — I don’t think even Eugene would argue with that, especially since he has done it himself!
 *Note [page 12]:* Only the recent paper by Lee & Chiang that basically seeks to explain all debris disk morphology with planets… but is it kosher to cite Eugene’s paper to explain a point that he disagrees with? :-) I think I’d just include the vague, unobjectionable statement and leave out the citation, honestly.  
    
  To Eugene’s second point: In Plavchan’s email to me about the planet, he said: “Such a planet would be an excellent explanation for the moving structures observed further out in the disk with direct imaging. While that's a very compelling story to tell, I haven't convinced myself yet that the planet is genuine.” Going back to look at the Sezestre (2017), the best-fit planet separation is ~8 au, but separations <= 5 au are still acceptable. Thus it doesn’t seem so crazy to link Plavchan’s planet to the moving features. I suspect Eugene is disagreeing so strongly because he doesn’t believe the conclusions of Sezestre?  
    
  I’d propose rewriting the sentences Eugene doesn’t like in a more conservative manner (as outlined above), carefully outlining the logic and doing a better job of citing things. I’m not particularly interested in getting involved in an email discussion with Eugene, so ideally Eugene can just add his comments to the next draft.. but maybe I should email him about it directly, since he seems so opinionated? I won’t do this last part unless you think it’s necessary. 
 *Note [page 12]:* Yes, I think it’s important to be fair to papers that are already refereed and published, even if Eugene disagrees with them (I mean, Eugene is a smart guy, so if he disagrees with them there’s probably a good reason, but it is still responsible to cite the relevant literature even if we do it in a circumspect way).  I think it’s fine to say something like “The < 1au separation of the potential planet detection by Plavchan et al. is far interior to the preferred 8au separation that Sezestre et al. (2017) derive to explain the origin of the fast-moving features, but separations < 5 au are not ruled out by their analysis.”  

- move to discussion of ring?  
  **Cail:** Someone (I think Kevin?) suggested moving the discussion of Plavchan’s planet to, well, the discussion. Unless you suggest otherwise, this is what I’ll do!
 *Note [page 12]:* Sure!  Sounds like a good idea. 
